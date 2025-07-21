#!/bin/bash
trycycler cluster --assemblies assemblies/*.fasta --reads bb_ont_filtered.fastq --out_dir trycycler
trycycler msa --cluster_dir trycycler/cluster_001
trycycler partition --reads bb_ont_filtered.fastq --cluster_dir trycycler/cluster_001
trycycler consensus --cluster_dir trycycler/cluster_001
