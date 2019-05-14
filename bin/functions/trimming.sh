#!/bin/bash

trimming() {

  mkdir -p $PWD/TRIMMING/1U2U
  for i in $(ls *fastq.gz | cut -d\_ -f1,2 | sort | uniq )
  do
    trimmomatic PE -phred33 -threads $(nproc) $i\_R1.fastq.gz $i\_R2.fastq.gz \
       $i\_R1.trim.fastq.gz $i.1U.trim.fastq.gz \
       $i\_R2.trim.fastq.gz $i.2U.trim.fastq.gz \
       SLIDINGWINDOW:4:20 MINLEN:70 &> $i.trim.log
    mv *U.trim.fastq.gz TRIMMING/1U2U && mv *trim*gz TRIMMING

    if [ $? -eq 0 ]; then
        rm $i.trim.log
    fi
  done
}
