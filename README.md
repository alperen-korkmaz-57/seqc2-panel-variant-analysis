# Targeted cancer panel variant analysis (SEQC2 Sample B)

A learning project: a reproducible germline variant-calling pipeline for targeted panel sequencing data, built on a public reference sample, with an R/Quarto report.

**Report:** https://alperen-korkmaz-57.github.io/seqc2-panel-variant-analysis/

## Overview

| | |
|---|---|
| Data | SEQC2 Oncopanel Sequencing Working Group, NCBI SRA SRP295113 (BioProject PRJNA677997) |
| Sample | Sample B, a normal male reference cell line (Agilent Human Reference DNA) |
| Panel | Burning Rock OncoScreen Plus: 520 cancer-related genes, about 1.6 Mb of target regions |
| Library | `SampleB_BRP1_ST28_100ng_LIB3` (one of four technical replicates) |
| Sequencing | Illumina NovaSeq 6000, 2 × 150 bp paired-end |
| Reference | UCSC hg19, primary chromosomes only |

I used a public reference sample on purpose. Patient sequencing data cannot be shared, and open data lets anyone rerun the analysis and check the results.

## Pipeline

All steps are in [`scripts/pipeline.sh`](scripts/pipeline.sh). Each step is skipped if its output already exists.

| Step | Tool | Purpose |
|---|---|---|
| 1 | FastQC | Quality control of raw reads |
| 2 | fastp | Trim poly-G tails (a NovaSeq artefact) and adapters |
| 3 | bwa mem, samtools | Align to hg19, add read groups, sort and index |
| 4 | samtools markdup | Mark PCR duplicates |
| 5 | mosdepth | Depth over target regions; on-target rate |
| 6 | bcftools mpileup / call | Germline variant calling within the target regions |
| 7 | bcftools norm | Split multi-allelic sites, left-align indels |
| 8 | bcftools filter | Soft filter: QUAL < 30, depth < 20 or VAF < 0.2 → `LowQual` |
| 9 | SnpEff | Functional annotation |
| 10 | bcftools query, samtools faidx | Export tables for the R report |

The report ([`report/report.qmd`](report/report.qmd)) reads the pipeline outputs and builds all tables and figures from them.

## Main results

| Metric | Value |
|---|---|
| Read pairs (raw / after trimming) | 39.4 M / 38.9 M (98.7%) |
| Mapped reads | 99.96% |
| Duplicate rate | 57.1% (estimated library size 19.2 M fragments) |
| On-target rate (unique reads) | 77.3% |
| Mean depth over targets | 1,493x |
| Target bases at ≥ 100x | 99.67% |
| Variants (PASS) | 1,246 (1,121 SNVs, 125 indels) |
| Ti/Tv (PASS) | 2.32 |

## What I learned

1. **Default settings can fail silently.** At about 1,500x, `bcftools mpileup` with default settings called only 2 indels, because it looks at no more than 250 reads per position when calling indels. Raising the limit (`-L 10000`) gave 150 indels, with no change in the SNV calls. I found the problem only because the indel fraction was biologically implausible.
2. **"No reads" and "reads with no confident position" are different.** CRLF2 showed 0x coverage. It lies in PAR1, a region identical on chrX and chrY, so its reads receive mapping quality 0 and are removed by the MAPQ ≥ 20 filter. Without the filter, one CRLF2 region had 4,032 reads.
3. **Short reads struggle in repeats.** 69% of indels, but only 1.6% of SNVs, lie in homopolymers or dinucleotide repeats. Several of these sites showed three alleles at once. Most of the seven HIGH-impact calls also looked doubtful on closer inspection. I report them as limitations, not as findings.

## Limitations

Germline calling only; CRLF2 and HIST2H3C cannot be called; indels in repeats are unreliable; only one library was analysed; no comparison with the SEQC2 truth set yet. See section 6 of the report for details.

## Repository structure

```
├── scripts/pipeline.sh     # the full pipeline (bash)
├── report/report.qmd       # R/Quarto report
├── docs/index.html         # rendered report (GitHub Pages)
├── results/                # small pipeline outputs used by the report
├── environment.yml         # conda environment for the pipeline
└── README.md
```

Large files (FASTQ, reference genome, BAM) are not included in the repository.

## How to run

**1. Environment** (Linux or WSL2, conda/mamba):

```bash
mamba env create -f environment.yml
conda activate panel
mkdir -p data ref regions
```

**2. Input data**

- FASTQ files for library `SampleB_BRP1_ST28_100ng_LIB3` from ENA (project PRJNA677997), saved as `data/SampleB_BRP1_ST28_LIB3_R1.fastq.gz` and `data/SampleB_BRP1_ST28_LIB3_R2.fastq.gz`. Check the md5 sums against ENA.
- Panel BED file `BRP_hg19.bed` from figshare ([doi:10.6084/m9.figshare.19128005](https://doi.org/10.6084/m9.figshare.19128005)), saved in `regions/`.
- Reference genome:

```bash
cd ref
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
samtools faidx hg19.fa
samtools faidx hg19.fa chr{1..22} chrX chrY chrM > hg19_primary.fa
samtools faidx hg19_primary.fa
bwa index hg19_primary.fa          # about one hour
snpEff download hg19
cd ..
```

**3. Pipeline**

```bash
bash scripts/pipeline.sh 2>&1 | tee results/pipeline.log
```

**4. Report** (R ≥ 4.3 with Quarto; packages: tidyverse, jsonlite, knitr, DT). Open `report/report.qmd` in RStudio and click Render, or run:

```bash
quarto render report/report.qmd
```

## Use of AI assistance

I built this project step by step with guidance from an AI assistant (Claude). I wrote or edited the code myself under its review, ran every step, and checked the outputs against expected values. The explanatory text in the report was written together with the assistant. I am still learning many of the underlying methods; this project documents that learning process.

## Data source

SEQC2 Oncopanel Sequencing Working Group data set: Gong B, et al. *Scientific Data* 9, 288 (2022). https://doi.org/10.1038/s41597-022-01359-6

## About

I am a final-year biology student at Gazi University learning bioinformatics. Feedback is welcome.
