#!/usr/bin/env bash

# SEQC2 targeted cancer panel - variant analysis pipeline
#
# sample    : SampleB_BRP1_ST28_LIB3 (Burning Rock OncoScreen Plus, 520 genes)
# reference : UCSC hg19, primary chromosomes only (chr1-22, X, Y, M)
# targets   : regions/BRP_hg19.bed (SEQC2, figshare)
#
# usage     : run from the project folder (~/seqc2_panel):
#               conda activate panel
#               bash scripts/pipeline.sh

# error handling:
#   -e          stop the script if any command fails
#   -u          stop the script if an undefined variable is used
#   -o pipefail also catch failures in the middle of piped ("|") commands

set -euo pipefail

#settings: to analyse another sample, change the variables

SAMPLE="SampleB_BRP1_ST28_LIB3" 
THREADS=8
BED="regions/BRP_hg19.bed"
REF="ref/hg19_primary.fa"
R1="data/${SAMPLE}_R1.fastq.gz"
R2="data/${SAMPLE}_R2.fastq.gz"
OUT="results"
FASTQC_DIR="${OUT}/fastqc"
FASTP_DIR="${OUT}/fastp"
ALIGN_DIR="${OUT}/alignment"
COV_DIR="${OUT}/coverage"
VAR_DIR="${OUT}/variants"
ANN_DIR="${OUT}/annotation"
mkdir -p "$FASTQC_DIR" "$FASTP_DIR" "$ALIGN_DIR" "$COV_DIR" "$VAR_DIR" "$ANN_DIR"

echo "variant analysis pipeline starting..."

echo "checking required tools..."

for tool in fastqc fastp bwa samtools bcftools mosdepth snpEff; do
    if command -v "$tool" > /dev/null 2>&1; then
        echo " - $tool: found"
    else
        echo "ERROR: $tool not found. Did you run 'conda activate panel'?"
        exit 1
    fi
done

echo "All required tools found."

for f in "$R1" "$R2" "$BED" "$REF" "${REF}.fai" "${REF}.bwt"; do
    [[ -s "$f" ]] || { echo "ERROR: missing file: $f"; exit 1; }
done
echo "All input files found."

#quality control of raw reads 
echo "starting quality control of raw reads..."
if [[ -s "${FASTQC_DIR}/${SAMPLE}_R1_fastqc.html" && -s "${FASTQC_DIR}/${SAMPLE}_R2_fastqc.html" ]]; then
    echo "fastQC reports already exist, skipping."
else
    fastqc --threads "$THREADS" --outdir "$FASTQC_DIR" "$R1" "$R2"
    echo "quality control finished. see ${FASTQC_DIR}/ for the reports."
fi

# Read trimming
echo "starting read trimming..."
TRIM_R1="${FASTP_DIR}/${SAMPLE}_trimmed_R1.fastq.gz"
TRIM_R2="${FASTP_DIR}/${SAMPLE}_trimmed_R2.fastq.gz"

if [[ -s "$TRIM_R1" && -s "$TRIM_R2" ]]; then
    echo "trimmed reads already exist, skipping."
else
    fastp \
        -i "$R1" -I "$R2" \
        -o "$TRIM_R1" -O "$TRIM_R2" \
        --trim_poly_g \
        --detect_adapter_for_pe \
        --thread "$THREADS" \
        --html "${FASTP_DIR}/${SAMPLE}_fastp.html" \
        --json "${FASTP_DIR}/${SAMPLE}_fastp.json"
    echo "fastp finished. Report: ${FASTP_DIR}/${SAMPLE}_fastp.html"
fi

echo "starting alignment..."
SORTED_BAM="${ALIGN_DIR}/${SAMPLE}.sorted.bam"
if [[ -s "$SORTED_BAM" ]]; then
    echo "Sorted BAM already exists, skipping."
else
    bwa mem -t "$THREADS" \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
        "$REF" "$TRIM_R1" "$TRIM_R2" \
      | samtools fixmate -m - - \
      | samtools sort -@ "$THREADS" -o "$SORTED_BAM" -
    samtools index "$SORTED_BAM"
    echo "Alignment finished: $SORTED_BAM"
fi

# Duplicate marking
# Why: PCR copies of the same DNA fragment are not independent evidence.
#      Marking them prevents a single PCR error from looking like a
#      well-supported variant. Reads are flagged, not removed.

echo "Starting duplicate marking..."
MARKDUP_BAM="${ALIGN_DIR}/${SAMPLE}.markdup.bam"

if [[ -s "$MARKDUP_BAM" ]]; then
    echo "Duplicate-marked BAM already exists, skipping."
else
    samtools markdup -@ "$THREADS" \
        -f "${ALIGN_DIR}/${SAMPLE}.markdup_stats.txt" \
        "$SORTED_BAM" "$MARKDUP_BAM"
    samtools index "$MARKDUP_BAM"
    samtools flagstat -@ "$THREADS" "$MARKDUP_BAM" > "${ALIGN_DIR}/${SAMPLE}.flagstat.txt"
    echo "Duplicate marking finished: $MARKDUP_BAM"
fi

echo "Starting coverage analysis..."
COV_PREFIX="${COV_DIR}/${SAMPLE}"

if [[ -s "${COV_PREFIX}.mosdepth.summary.txt" ]]; then
    echo "Coverage results already exist, skipping."
else
    mosdepth --threads 4 --by "$BED" \
        --thresholds 100,250,500,1000 \
        -Q 20 --no-per-base \
        "$COV_PREFIX" "$MARKDUP_BAM"

    # On-target rate: unique mapped reads inside the BED / all unique mapped reads
    # -F 0xD04 excludes unmapped, secondary, duplicate and supplementary reads
    TOTAL_READS=$(samtools view -c -@ "$THREADS" -F 0xD04 "$MARKDUP_BAM")
    ONTARGET_READS=$(samtools view -c -@ "$THREADS" -F 0xD04 -L "$BED" "$MARKDUP_BAM")
    echo -e "total_reads\ton_target_reads\n${TOTAL_READS}\t${ONTARGET_READS}" \
        > "${COV_DIR}/${SAMPLE}.on_target.tsv"

    echo "Coverage analysis finished: ${COV_DIR}/"
fi
#Variant Calling
echo "Starting variant calling..."
RAW_VCF="${VAR_DIR}/${SAMPLE}.raw.vcf.gz"

if [[ -s "$RAW_VCF" ]]; then
    echo "Raw VCF already exists, skipping."
else
    bcftools mpileup --threads "$THREADS" \
        -f "$REF" -R "$BED" \
        -d 10000 -L 10000 -q 20 -Q 20 \
        -a FORMAT/AD,FORMAT/DP \
        -Ou "$MARKDUP_BAM" \
      | bcftools call --threads "$THREADS" -m -v -Oz -o "$RAW_VCF"
    bcftools index "$RAW_VCF"
    bcftools stats "$RAW_VCF" > "${VAR_DIR}/${SAMPLE}.raw.stats.txt"
    echo "Variant calling finished: $RAW_VCF"
fi


#Normalizing VCF
# Why: -m -any splits multiallelic sites into one line per alternative
#      allele, so each allele can be filtered on its own (e.g. stutter
#      alleles at homopolymers). -f REF left-aligns indels, so the same
#      indel is always written the same way.

echo "Starting normalization..."
NORM_VCF="${VAR_DIR}/${SAMPLE}.normalized.vcf.gz"

if [[ -s "$NORM_VCF" ]]; then
    echo "Normalized VCF already exists, skipping."
else
    bcftools norm -f "$REF" -m -any -Oz -o "$NORM_VCF" "$RAW_VCF"
    bcftools index "$NORM_VCF"
    bcftools stats "$NORM_VCF" > "${VAR_DIR}/${SAMPLE}.normalized.stats.txt"
    echo "Normalizing VCF finished! $NORM_VCF"
fi 
#Filtering
#Why: Not every variant is real. Low-quality calls, few reads positions and alleles supported by only small fraction of reads 
#     are usally sequencşing or PCR errors.
#  +fill-tags -t FORMAT/VAF   adds the variant allele fraction
#                              (alt reads / all reads) to every variant
#   QUAL<30         call confidence too low (error chance above 1 in 1000)
#   FORMAT/DP<20    too few reads to decide the genotype
#   FORMAT/VAF<0.2  germline variants are expected near 0.5 (het) or 1.0 (hom);
#                   much lower fractions are likely artefacts
#   -s LowQual      variants are not removed, they are labelled LowQual in the
#                   FILTER column; the rest are labelled PASS

echo "Starting filtering..."
FILT_VCF="${VAR_DIR}/${SAMPLE}.filtered.vcf.gz"

if [[ -s "$FILT_VCF" ]]; then
   echo "Filtered VCF already exists, skipping..."
else
   bcftools +fill-tags "$NORM_VCF" -Ou -- -t FORMAT/VAF | bcftools filter -s LowQual -e 'QUAL<30 || FORMAT/DP<20 || FORMAT/VAF<0.2' -Oz -o "$FILT_VCF" -
   bcftools index "$FILT_VCF"
   bcftools stats "$FILT_VCF" > "${VAR_DIR}/${SAMPLE}.filtered.stats.txt"
   bcftools stats -f PASS "$FILT_VCF" > "${VAR_DIR}/${SAMPLE}.pass.stats.txt"
   echo "Filtering VCF finished!"
fi

# Annotation
# Why: Find which gene each variant falls in and how it changes the protein
#      (e.g. missense, frameshift, synonymous).
#   hg19        SnpEff database for the genome build we used (UCSC, "chr"
#               names, matches our reference)
#   -Xmx8g      SnpEff runs on Java; give Java 8 GB of memory
#   -stats      write SnpEff's HTML summary report
# Limitation: SnpEff only predicts the effect of a variant. Whether a
#             variant is linked to disease needs databases such as ClinVar.

echo "Starting VCF Annotation..."
ANN_VCF="${ANN_DIR}/${SAMPLE}.annotated.vcf.gz"

if [[ -s "$ANN_VCF" ]]; then
   echo "Annotated VCF already exists, skipping..."
else
   snpEff -Xmx8g ann -v -stats "${ANN_DIR}/${SAMPLE}.snpEff.html" hg19 "$FILT_VCF" | bcftools view -Oz -o "$ANN_VCF"
   bcftools index "$ANN_VCF"
   echo "Annotation is complete!"
fi

# ---------------------------------------------------------------------
# STEP 10: Export tables for the R report
# Why: R reads plain tables more easily than VCF files.
#   variants.tsv : one row per variant (position, filter, genotype,
#                  depth, VAF and the SnpEff annotation)
#   context.tsv  : reference sequence from 10 bp before to 20 bp after
#                  each variant, used in the report to flag variants in
#                  homopolymers and short repeats
# ---------------------------------------------------------------------
echo "10: Export tables for the report"

VARIANTS_TSV="${VAR_DIR}/${SAMPLE}.variants.tsv"
CONTEXT_TSV="${VAR_DIR}/${SAMPLE}.context.tsv"

if [[ -s "$VARIANTS_TSV" && -s "$CONTEXT_TSV" ]]; then
    echo "Report tables already exist, skipping."
else
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%QUAL\t[%GT]\t[%DP]\t[%VAF]\t%INFO/ANN\n' \
        "$ANN_VCF" > "$VARIANTS_TSV"

    # Indels are left-aligned, so a repeat usually extends to the right
    bcftools query -f '%CHROM\t%POS\n' "$ANN_VCF" \
        | awk '{print $1 ":" $2-10 "-" $2+20}' \
        > "${VAR_DIR}/${SAMPLE}.context_regions.txt"

    # paste - - joins each FASTA header with its sequence on one line
    samtools faidx "$REF" -r "${VAR_DIR}/${SAMPLE}.context_regions.txt" \
        | paste - - > "$CONTEXT_TSV"

    echo "Tables written: $VARIANTS_TSV and $CONTEXT_TSV"
fi
