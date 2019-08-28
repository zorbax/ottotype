#!/bin/bash

lis_type(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p {SRST2Ec_$run_name,OUTPUT,RESULTS}

  repo="https://raw.githubusercontent.com/CNRDOGM"
  docker_cmd="docker run --rm -it -v $(pwd):/data -u $(id -u):$(id -g) -w /data"

  wget -q -nc $repo/srst2/master/data/EcOH.fasta

  $docker_cmd srst2 getmlst.py --species "Listeria monocytogenes" &> /dev/null
  
  $docker_cmd srst2 srst2 --log --output /data/SRST2_Lys --input_pe *fastq.gz \
      --forward R1 --reverse R2 --mlst_db Listeria_monocytogenes.fasta \
      --mlst_definitions lmonocytogenes.txt \
      --mlst_delimiter '_' --threads $(nproc) &> /dev/null
