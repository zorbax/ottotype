#!/bin/bash

assembly_spades() {
  mkdir -p ASSEMBLY
  memory=`awk '{ printf "%.2f", $2/1024/1024 ; exit}' /proc/meminfo | cut -d\. -f1`
  for i in $(ls *trim.fastq.gz | cut -d\_ -f1,2 | sort | uniq)
  do
    spades.py --pe1-1 $i\_R1.trim.fastq.gz --pe1-2 $i\_R2.trim.fastq.gz \
              --pe1-s 1U2U/$i.1U.trim.fastq.gz --pe1-s 1U2U/$i.2U.trim.fastq.gz \
              -o $i\_spades -t $(nproc) -m $memory &>/dev/null
  # --careful --only-assembler using sga corrected reads
    find $i\_spades -maxdepth 2 -type f -name 'scaffolds.fasta' -exec cp {} $PWD/$i.tmp \;
    cat $i.tmp | sed ':a;N;/^>/M!s/\n//;ta;P;D' | \
       awk '/^>/ { getline seq } length(seq) >500 { print $0 "\n" seq }' > ASSEMBLY/$i-spades-assembly.fa
    rm -rf $i\_spades $i.tmp
  done
}
