#!/bin/bash - 

grep < ARGannot_r3.fasta \> | sed 's/ /;/' | \
    tr ';' '\t' | awk -v OFS='\t' '{ print $5, $4}' | \
    tee >( sort -k1,2 > antibiotics_code.v3.raw.tsv) \
    >( cut -f1 | sort | uniq > antibiotics_categories.txt) \
    > /dev/null

