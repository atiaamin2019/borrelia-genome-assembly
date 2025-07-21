#!/bin/bash
for dir in ./contig*/; do
    if [ -d "$dir" ]; then
        large_contig=$(find "$dir" -maxdepth 1 -type f -name "bb*.fasta" | head -n 1)
        small_contig_1=$(find "$dir" -maxdepth 1 -type f -name "*_L_end.fasta" | head -n 1)
        small_contig_2=$(find "$dir" -maxdepth 1 -type f -name "*_R_end.fasta" | head -n 1)

        output_dir="${dir}/merging_of_contigs"
        python merging_contigs_blast.py "$large_contig" "$small_contig_1" "$small_contig_2" "$output_dir"
    fi
done
