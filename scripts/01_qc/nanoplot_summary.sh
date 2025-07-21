#!/bin/bash
NanoPlot --fastq bb_ont_filtered.fastq -o ./bb_ont_filtered_QC --tsv_stats -c magenta -cm plasma \
         -f tiff --plots kde hex dot --N50 --title bb_ONT --dpi 300 --verbose
