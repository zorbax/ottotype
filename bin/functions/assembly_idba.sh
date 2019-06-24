#!/bin/bash

assembly_idba(){

  for r1 in TRIMMING/*R1.trim.fastq.gz
  do
    mkdir -p ASSEMBLY

    r2="${r1/R1/R2}"
    name="${r1##*/}"; name="${name%%_R1*}"

    echo "${name}"
    sga preprocess -q 25 -f 20 --pe-mode=1 ${r1} ${r2} \
        > ${name}_12.t.pp.fq 2> /dev/null
    median=$(cat ${name}_12.t.pp.fq | awk 'NR%4==2 { print length }' | head -20000 | \
            Rscript -e 'summary (as.numeric (readLines ("stdin")))' | tail -n+2 | \
            awk '{ print $3 }'  | cut -d\. -f1)

    if [[ $median -le 200 ]]; then
      sga index -t $(nproc) -a ropebwt ${name}_12.t.pp.fq \
          > ${name}.index.out 2> ${name}.index.err
      sga correct -t $(nproc) ${name}_12.t.pp.fq -o ${name}_12.t.pp.ec.fq \
          > ${name}.correct.out 2> ${name}.correct.err
    else
      sga index -t $(nproc) -a sais ${name}_12.t.pp.fq \
          > ${name}.index.out 2> ${name}.index.err
      sga correct -t $(nproc) ${name}_12.t.pp.fq -o ${name}_12.t.pp.ec.fq \
          > ${name}.correct.out 2> ${name}.correct.err
    fi

    sed -n '1~4s/^@/>/p;2~4p' ${name}_12.t.pp.ec.fq > ${name}_12.t.pp.ec.fa
    medianec=$(cat ${name}_12.t.pp.ec.fa | awk '$0 !~ /^>/ { print length }' | \
              head -20000 | Rscript -e 'summary (as.numeric (readLines ("stdin")))' | \
              tail -n+2 | awk '{ print $3}'  | cut -d\. -f1)

    if [[ $medianec -ge 128 ]]; then
      idba_ud500 -r ${name}_12.t.pp.ec.fa -o ${name} --mink 35 --maxk 249 \
          --num_threads $(nproc) --min_pairs 2 > /dev/null
    else
      idba_ud -r ${name}_12.t.pp.ec.fa -o ${name} --mink 35 --maxk 124 \
          --num_threads $(nproc) --min_pairs 2 > /dev/null
    fi

    if [[ $? != 0 ]]; then
       cp ${name}/contig.fa ${name}-idba-assembly.fa
    else
       cp ${name}/scaffold.fa ${name}-idba-assembly.fa
    fi

    rm -rf *pp.{fq,bwt,sai,rbwt} *out *err *rsai *ec.{fq,fa} ${name}

    cat ${name}-idba-assembly.fa | sed ':a;N;/^>/M!s/\n//;ta;P;D' | \
         awk '/^>/ { getline seq } length(seq) >500 { print $0 "\n" seq }' \
         > ASSEMBLY/${name}-idba-assembly.fa
    rm ${name}-idba-assembly.fa
  done
}
