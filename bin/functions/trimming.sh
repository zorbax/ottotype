#!/bin/bash

trimming() {

  mkdir -p TRIMMING/1U2U
  for r1 in *R1.fastq.gz
  do
    r2=${r1/R1/R2}
    name=${r1%%_R1*}
    trimmomatic PE -phred33 -threads $(nproc) $r1 $r2 \
              TRIMMING/${name}_R1.trim.fastq.gz TRIMMING/1U2U/${name}.1U.trim.fastq.gz \
              TRIMMING/${name}_R2.trim.fastq.gz TRIMMING/1U2U/${name}.2U.trim.fastq.gz \
              SLIDINGWINDOW:4:20 MINLEN:125 &> ${name}.trim.log

    if [ $? -eq 0 ]; then
      rm *trim.log
    fi
  done
}
