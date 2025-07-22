
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
bb_illumina_R1_P.fastq.gz
bb_illumina_R2_P.fastq.gz

1.2 Nanopore Reads (using filtlong + NanoPlot)
```bash
bash scripts/01_qc/filtlong_nanopore.sh
bash scripts/01_qc/nanoplot_summary.sh
```
Output:
bb_ont_filtered.fastq
QC plots under bb_ont_filtered_QC/

## 2. Linear Chromosome Assembly
2.1 Subsample Reads (Trycycler)
```
bash scripts/02_assembly/trycycler_subsample.sh
```
2.2 Assemble with Flye, Canu, and Raven
```
bash scripts/02_assembly/trycycler_assembly.sh
```
2.3 Cluster, Align, Partition, Consensus
```
bash scripts/02_assembly/trycycler_clustering.sh
```
Output:
7_final_consensus.fasta in trycycler/cluster_001/

2.4 Polishing (Medaka + Polypolish)
```
bash scripts/02_assembly/polishing.sh
```
Final Output:
trycycler_medaka_polypolish.fasta

## 3. Telomeric Trimming of Linear Chromosome
3.1 Prepare trimming_position.txt
Manually inspect telomere ends and create:
```
bb1 11463 921975
bb3 19033 921811
...
```
3.2 Run trimming script
```
python scripts/03_trimming/trimming_based_on_inversion_position.py
```
Output:
trimmed_chromosomes.fasta

## 4. Plasmid Assembly (Plassembler)
4.1 Run Plassembler
```
bash scripts/04_plasmid_assembly/plassembler.sh
```
Output:
plassembler_output/plasmids.fasta

4.2 Split circular vs linear based on summary.tsv
Manually separate into:
	•	circular_plasmid_contigs.fasta
	•	linear_plasmid_contigs.fasta
4.3 Circularize plasmids (optional manual fix)
```
python scripts/04_plasmid_assembly/simple_circularise.py circular_plasmid_contigs.fasta circularized_contigs.fasta -r 10 -min 5000
```
## 5. Telomeric Extension for Linear Plasmids
5.1 Extract telomeric reads
```
bash scripts/05_telomere_extension/extract_telomeric_reads.sh
```
5.2 Assemble telomeric ends (Flye)
```
bash scripts/05_telomere_extension/flye_telomere_assembly.sh
```
5.3 Merge telomeric ends with main contig
Make sure:
```
contig_x_L_end.fasta
contig_x_R_end.fasta
bbX_contig_x_lpYY.fasta
```
are all in their contig_* folder.
Then run:
```
bash scripts/05_telomere_extension/batch_merge_contigs.sh
```
5.4 Trim merged linear plasmids
Update trimming_position.txt, then:
```
python scripts/03_trimming/trimming_based_on_inversion_position.py
```
## 6. Genome Annotation (Bakta)
6.1 Download Bakta DB
```
bash scripts/06_annotation/bakta_annotation.sh
```
Output files:
annotation.gbff, annotation.tsv, annotation.json in results/annotation/

📁 Output Summary
	•	trimmed_chromosomes.fasta: Final linear chromosome
	•	circularized_contigs.fasta: Circular plasmids
	•	trimmed_linear_plasmids.fasta: Final linear plasmids
	•	bb_chromosome_plasmids.fa: Combined for annotation
	•	annotation/: Annotated genome files

## Acknowledgments
Telomere-to-telomere assembly detects genomic diversity in Canadian strains of Borrelia burgdorferi
Atia Amin, Ana Victoria Ibarra Meneses, Simon Gagnon, Georgi Merhi, Martin Olivier, Momar Ndao, Mathieu Blanchette, Christopher Fernandez Prada, David Langlais
bioRxiv 2025.07.18.665644; doi: https://doi.org/10.1101/2025.07.18.665644

