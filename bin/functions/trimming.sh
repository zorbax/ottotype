#!/bin/baslocal h

trimming() {

  mkdir -p TRIMMING/1U2U
  local adapters="/mnt/disk1/bin/adapters/NexteraPE-PE.fa"
  for r1 in *R1.fastq.gz
  do
    local r2=${r1/R1/R2}
    local name=${r1%%_R1*}
    trimmomatic PE -phred33 -threads "$(nproc)" $r1 $r2 \
              TRIMMING/${name}_R1.trim.fastq.gz TRIMMING/1U2U/${name}.1U.trim.fastq.gz \
              TRIMMING/${name}_R2.trim.fastq.gz TRIMMING/1U2U/${name}.2U.trim.fastq.gz \
              ILLUMINACLIP:${adapters}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:70
  done
}
