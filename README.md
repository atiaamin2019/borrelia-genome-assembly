<<<<<<< HEAD
# Borrelia Genome Assembly Pipeline

This repository implements a comprehensive telomere-to-telomere hybrid genome assembly pipeline for *Borrelia burgdorferi* using both Illumina short reads and Oxford Nanopore long reads.

## Pipeline Overview
1. Quality control of Illumina & ONT reads
2. Assembly of linear chromosome using Trycycler
3. Polishing with Medaka and Polypolish
4. Telomeric trimming of linear replicons
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

## Acknowledgments
Based on the protocol: *A Comprehensive Protocol for Telomere-to-telomere Genome Assembly of a Spirochete Bacteria Borrelia burgdorferi* (Amin et al.)
=======
# borrelia-genome-assembly
Telomere-to-telomere genome assembly pipeline for Borrelia burgdorferi
>>>>>>> 28a9642930f9f7b7088c906894b60a8ec09bd244
