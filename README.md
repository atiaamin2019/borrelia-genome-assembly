
# Borrelia Genome Assembly Pipeline

This repository implements a comprehensive telomere-to-telomere hybrid genome assembly pipeline for *Borrelia burgdorferi* using both Illumina short reads and Oxford Nanopore long reads.

## Pipeline Overview
1. Quality control of Illumina & ONT reads
2. Assembly of linear chromosome using Trycycler
3. Polishing with Medaka and Polypolish
4. Telomeric trimming of linear chromosome
5. Plasmid assembly using Plassembler
6. Telomeric extension of linear plasmids
7. Annotation using Bakta

## Installation
```bash
conda env create -f environment.yml
conda activate borrelia_assembly
```

## Usage
Follow the script order in the `scripts/` directory for each stage of the pipeline.

## Directory Structure
See `docs/protocol_summary.md` for step-by-step descriptions.

## Step-by-Step Pipeline
## 1. Quality Control of Reads
1.1 Illumina Paired-End Reads (using fastp)

```bash
bash scripts/01_qc/fastp_illumina.sh
```
Output:
	•	bb_illumina_R1_P.fastq.gz, bb_illumina_R2_P.fastq.gz

1.2 Nanopore Reads (using filtlong + NanoPlot)
```bash
bash scripts/01_qc/filtlong_nanopore.sh
bash scripts/01_qc/nanoplot_summary.sh
```
Output:
	•	bb_ont_filtered.fastq
	•	QC plots under bb_ont_filtered_QC/


## Acknowledgments
Telomere-to-telomere assembly detects genomic diversity in Canadian strains of Borrelia burgdorferi (Amin et al. 2025)

