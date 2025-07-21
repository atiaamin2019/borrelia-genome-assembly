#!/bin/bash
threads=16
mkdir assemblies

flye --nano-hq read_subsets/sample_01.fastq --threads $threads --out-dir assembly_01 && \
    cp assembly_01/assembly.fasta assemblies/assembly_01.fasta

canu -p assembly_02 -d canu_assembly_02 genomeSize=1m -nanopore read_subsets/sample_02.fastq \
     useGrid=false maxThreads=$threads && \
     cp canu_assembly_02/assembly_02.contigs.fasta assemblies/assembly_02.fasta

raven --threads $threads --disable-checkpoints read_subsets/sample_03.fastq > assemblies/assembly_03.fasta
