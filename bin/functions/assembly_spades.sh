#!/bin/bash

assembly_spades() {

  mkdir -p ASSEMBLY/RAW
  memory=$(awk '{ printf "%.2f", $2/1024/1024 ; exit}' /proc/meminfo | cut -d\. -f1)

  for r1 in TRIMMING/*R1.trim.fastq.gz
  do
    r2="${r1/R1/R2}"
    name="${r1##*/}"; name="${name%%_R1*}"
    spades.py --pe1-1 $r1 --pe1-2 $r2 \
              --pe1-s TRIMMING/1U2U/${name}.1U.trim.fastq.gz \
              --pe1-s TRIMMING/1U2U/${name}.2U.trim.fastq.gz \
              -o ${name}_spades -t $(nproc) -m $memory &>/dev/null

    if [ -s "${name}_spades/scaffolds.fasta" ]; then
      cp ${name}_spades/scaffolds.fasta ${name}_spades/${name}-spades-assembly.fa
    else
      cp ${name}_spades/contigs.fasta ${name}_spades/${name}-spades-assembly.fa
    fi

    assembly="${name}_spades/${name}-spades-assembly.fa"
    cat ${assembly} | sed ':a;N;/^>/M!s/\n//;ta;P;D' | \
        awk '/^>/ { getline seq } length(seq) >500 { print $0 "\n" seq }' \
        > ASSEMBLY/${name}-spades-assembly.fa

    rm -rf ${name}_spades
  done
}
