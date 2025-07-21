#!/bin/bash
fastp --in1 bb_illumina_R1.fastq.gz --in2 bb_illumina_R2.fastq.gz \
      --out1 bb_illumina_R1_P.fastq.gz --out2 bb_illumina_R2_P.fastq.gz \
      --unpaired1 bb_illumina_R1_U.fastq.gz --unpaired2 bb_illumina_R2_U.fastq.gz
