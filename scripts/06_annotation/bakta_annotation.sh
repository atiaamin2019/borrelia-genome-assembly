#!/bin/bash
bakta_db download --output config/bakta_db --type light

bakta -v --db config/bakta_db --output results/annotation --genus Borreliella \
      --species burgdorferi --force --compliant --threads 16 bb_chromosome_plasmids.fa
