#!/bin/bash

# Running plassembler
plassembler \
--longreads bb_ont_filtered.fastq \
--shortreads1 bb_illumina_R1_P.fastq \
--shortreads2 bb_illumina_R2_P.fastq \
--output plassembler_output \
--threads 8
