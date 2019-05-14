#!/bin/bash

kmer_finder(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p kmerfinder\_$run_name RESULTS/

  DB="/mnt/disk1/bin/kmerfinder_DB/bacteria.organisms.ATGAC"

  for i in $PWD/ASSEMBLY/*assembly.fa
  do
    genome_name=`basename $i | cut -d\- -f1`
    findTemplate -i $i -t $DB -x ATGAC -w -o kmerfinder\_$run_name/$genome_name.kmer
  done

  for i in kmerfinder\_$run_name/*kmer
  do
    echo $i | cut -d\/ -f2 | cut -d\_ -f1
    cat $i | tail -n+2 | head -4 | awk -F'\t' -v OFS='\t' '{ print $1, $2, $5 }' | \
        sed 's/ //g; s/_/#/; s/\t/\&\t/; s/_.*&//; s/#/ /'
  done > RESULTS/kmerfinder\_$run_name\_all.txt
}
