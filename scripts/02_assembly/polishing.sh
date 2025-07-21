#!/bin/bash
threads=16

# Medaka
medaka_consensus -i trycycler/cluster_001/4_reads.fastq -d trycycler/cluster_001/7_final_consensus.fasta \
                 -o trycycler/cluster_001/medaka -m r104_e81_sup_g610

mv trycycler/cluster_001/medaka/consensus.fasta trycycler/cluster_001/trycycler_medaka.fasta

# Polypolish
bwa index trycycler_medaka.fasta
bwa mem -t $threads -a trycycler_medaka.fasta bb_illumina_R1_P.fastq.gz > alignments_1.sam
bwa mem -t $threads -a trycycler_medaka.fasta bb_illumina_R2_P.fastq.gz > alignments_2.sam

polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish polish trycycler_medaka.fasta filtered_1.sam filtered_2.sam > trycycler_medaka_polypolish.fasta
