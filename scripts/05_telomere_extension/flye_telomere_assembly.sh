#!/bin/bash
mkdir -p filtered_reads flye_output

# Filter reads longer than 6kb
for file in telomeric_end_reads/*.fasta; do
    base=$(basename $file)
    seqtk seq -L 6000 $file > filtered_reads/filtered_$base
done

# Run Flye for each
for fasta in filtered_reads/*.fasta; do
    name=$(basename $fasta .fasta)
    mkdir -p flye_output/$name
    flye -t 6 --nano-raw $fasta --genome-size 0.05m --out-dir flye_output/$name
done
