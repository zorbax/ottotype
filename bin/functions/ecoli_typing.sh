#!/bin/bash

ecoli_type(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p {SRST2Ec_${run_name},OUTPUT,RESULTS}

  repo="https://raw.githubusercontent.com/CNRDOGM"
  docker_cmd="docker run --rm -it -v $(pwd):/data -u $(id -u):$(id -g) -w /data"

  wget -q -nc $repo/srst2/master/data/EcOH.fasta

  $docker_cmd srst2 srst2 --log --output /data/SRST2_EcOH \
      --input_pe *fastq.gz --forward R1 --reverse R2 \
      --gene_db EcOH.fasta --threads $(nproc) &> /dev/null

   wget -q -nc $repo/srst2/master/data/{LEE_mlst.fasta,LEE_profiles.txt}
   $docker_cmd srst2 srst2 --log --output /data/SRST2_LEE \
       --input_pe *fastq.gz --forward R1 --reverse R2 \
       --mlst_db LEE_mlst.fasta --mlst_definitions LEE_profiles.txt \
       --threads $(nproc) &> /dev/null

  wget -q -nc "$repo/srst2/master/data/ARGannot_r3.fasta" -O ARGannot.fasta
  $docker_cmd srst2 srst2 --log --output /data/SRST2_ARG \
       --input_pe *fastq.gz --forward R1 --reverse R2 \
       --gene_db ARGannot.fasta --threads $(nproc) &> /dev/null

  $docker_cmd srst2 getmlst.py --species "Escherichia coli#1" &> /dev/null && \
  $docker_cmd srst2 srst2 --log --output /data/SRST2_Ec1 \
       --input_pe *fastq.gz --forward R1 --reverse R2 \
       --mlst_db Escherichia_coli#1.fasta --mlst_definitions ecoli.txt \
       --mlst_delimiter '_' --threads $(nproc) &> /dev/null

  $docker_cmd srst2 getmlst.py --species "Escherichia coli#2" &> /dev/null && \
  $docker_cmd srst2 srst2 --log --output /data/SRST2_Ec2 \
       --input_pe *fastq.gz --forward R1 --reverse R2 \
       --mlst_db Escherichia_coli#2.fasta --mlst_definitions ecoli_2.txt \
       --mlst_delimiter '_' --threads $(nproc) &> /dev/null

  find . -maxdepth 1 -iname "*results.txt" -type f -exec mv {} SRST2Ec_${run_name} \;
  rm -f EcOH* LEE* ARG* *tfa SRST2_*.{bam,pileup} Escherichia* ecoli* *log

  cat SRST2Ec_${run_name}/SRST2_EcOH__fullgenes__EcOH__results.txt | tail -n+2 | \
     sort -nk 1 | sed -E 's/_S[0-9]{1,}_//; s/_//' | awk -v OFS='\t' '{ print $1, $3, $NF }' | \
     sed -E 's/\t[A-Z0-9]{1,}.[0-9]{1};/\t/; s/\t/\tGen: /; s/;/: /' | \
     awk -F'\t' -v OFS='\t' '{x=$1; $1=""; a[x]=a[x]$0 } END { for (x in a) print x,a[x] }' | \
     sort -nk 1 | sed -E 's/\t{1,}Gen/\nGen/g; s/\t/; /g; s/polyermase/polymerase/' \
     > RESULTS/srst2Ec_${run_name}_EcOH.tsv

  cat SRST2Ec_${run_name}/SRST2_ARG__genes__ARGannot__results.txt | sed 's/[?*]//g' | \
     perl -pe 's/\t\-/#/g; s/\#{1,}//g;s/\_[0-9]{1,}//g'| perl -pe 's/_S[0-9]{1,}_//g;
     s/_//; s/f{3,}/\tFAILED/' | sed "s/''/-/" | tail -n+2 > RESULTS/srst2Ec_${run_name}_argannot.tsv

  translate.py RESULTS/srst2Ec_${run_name}_argannot.tsv RESULTS/antibioticsEc_${run_name}.tsv
  translate.py RESULTS/srst2Ec_${run_name}_argannot.tsv RESULTS/antibioticsEc_${run_name}_freq.tsv -freq

  cat SRST2Ec_${run_name}/SRST2_Ec1__mlst__Escherichia_coli#1__results.txt | \
      cut -d$'\t' -f 1-9 | sed 's/[?*]//g; /^$/d' | perl -pe 's/_S[0-9]{1,}_//g;
      s/_//' | sort | sed 's/ST/ST1/' > SRST2Ec_${run_name}/srst2Ec_${run_name}_achtman#1.tsv

  cat SRST2Ec_${run_name}/SRST2_Ec2__mlst__Escherichia_coli#2__results.txt | \
      cut -d$'\t' -f 1-9 | sed 's/[?*]//g; /^$/d' | perl -pe 's/_S[0-9]{1,}_//g;
      s/_//' | sort | sed 's/ST/ST2/' > SRST2Ec_${run_name}/srst2Ec_${run_name}_achtman#2.tsv

  paste SRST2Ec_${run_name}/srst2Ec_${run_name}_achtman#1.tsv \
      SRST2Ec_${run_name}/srst2Ec_${run_name}_achtman#2.tsv | \
      awk -v OFS='\t' '$1 == $10 { $10="" ; print}' \
      > RESULTS/srst12Ec_${run_name}_achtman12.tsv

  cat SRST2Ec_${run_name}/SRST2_LEE__mlst__LEE_mlst__results.txt | \
      cut -d$'\t' -f 1-9 | sed 's/[?*]//g; /^$/d' | perl -pe 's/_S[0-9]{1,}_//g;
      s/_//' | tail -n+2 | sort > RESULTS/srst2Ec_${run_name}_LEE.tsv

  mv SRST2Ec_${run_name} OUTPUT
}
