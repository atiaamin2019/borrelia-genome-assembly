# Protocol Summary: Telomere-to-Telomere Assembly of *Borrelia burgdorferi*

This document summarizes the genome assembly pipeline as outlined in the paper:
*A Comprehensive Protocol for Telomere-to-telomere Genome Assembly of a Spirochete Bacteria Borrelia burgdorferi.*

## Overview

The pipeline integrates long-read Oxford Nanopore and short-read Illumina sequencing to:
- Assemble the linear chromosome with telomeric resolution
- Assemble and curate both circular and linear plasmids
- Extend and trim telomeric regions for accuracy
- Perform genome annotation with `Bakta`

## Pipeline Stages

### 1. Quality Control
- **Illumina reads**: trimmed using `fastp`
- **Nanopore reads**: filtered using `filtlong` and analyzed via `NanoPlot`

### 2. Linear Chromosome Assembly
- Assembled using `Trycycler` from ONT reads
- Polished with `Medaka` (ONT) and `Polypolish` (Illumina)

### 3. Telomeric Trimming
- Hairpin telomere ends are identified using conserved sequences
- Trimmed using a BLAST-guided Python script

### 4. Plasmid Assembly
- Conducted using `Plassembler`
- Circular contigs are trimmed with `simple_circularise.py`
- Linear plasmids are extended at telomeric ends via `Flye`

### 5. Telomeric Extension
- ONT reads mapped to plasmid ends are extracted
- Left and right end contigs are assembled and merged into final plasmid contigs

### 6. Annotation
- Final genome annotated with `Bakta` using downloaded genus-specific database

## Outputs
- Telomere-resolved chromosome FASTA
- Circular and linear plasmid FASTAs
- Annotated genome in GenBank, TSV, and JSON formats

See individual scripts in `scripts/` for detailed steps.
