#!/bin/bash

assembly_spades() {

mkdir -p ASSEMBLY
memory=$(awk '{ printf "%.2f", $2/1024/1024 ; exit}' /proc/meminfo | cut -d\. -f1)

for r1 in TRIMMING/*R1.trim.fastq.gz
do
  r2="${r1/R1/R2}"
  name="${r1##*/}"; name="${name%%_R1*}"
  spades.py --pe1-1 $r1 --pe1-2 $r2 \
            --pe1-s TRIMMING/1U2U/${name}.1U.trim.fastq.gz \
            --pe1-s TRIMMING/1U2U/${name}.2U.trim.fastq.gz \
            -o ${name}_spades -t $(nproc) -m $memory &>/dev/null

  if [ $? != 0 ]; then
    cp ${name}/contig.fa ${name}-spades-assembly.fa
  else
    cp ${name}/scaffolds.fasta ${name}-spades-assembly.fa
  fi

  find ${name}_spades -maxdepth 2 -type f -name 'scaffolds.fasta' -exec cp {} ${name}.tmp \;

  cat ${name}.tmp | sed ':a;N;/^>/M!s/\n//;ta;P;D' | \
       awk '/^>/ { getline seq } length(seq) >500 { print $0 "\n" seq }' \
       > ASSEMBLY/${name}-spades-assembly.fa

  rm -rf ${name}_spades ${name}.tmp
done
}
