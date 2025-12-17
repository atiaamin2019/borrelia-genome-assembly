#!/bin/bash
# Align ONT reads to plasmid contigs
minimap2 -a linear_plasmid_contigs.fasta bb_ont_filtered.fastq > ONT_mapped_to_linear_plasmid.sam

# Convert to sorted BAM and index
samtools view -bS ONT_mapped_to_linear_plasmid.sam \
| samtools view -h -F 0x900 - \
| samtools sort -o ONT_mapped_to_linear_plasmid.sorted.bam - \
&& samtools index ONT_mapped_to_linear_plasmid.sorted.bam

# Filter reads based on BED regions
samtools view -b ONT_mapped_to_linear_plasmid.sorted.bam | bedtools intersect -a - -b regions.bed > ends_mapped_reads.bam
samtools index ends_mapped_reads.bam

# Extract telomeric reads
samtools view -b ends_mapped_reads.bam lp17:0-5000 | samtools fasta - > telomeric_end_reads/lp17_left_end_mapped_ONT_reads.fasta
samtools view -b ends_mapped_reads.bam lp17:16000-21000 | samtools fasta - > telomeric_end_reads/lp17_right_end_mapped_ONT_reads.fasta
