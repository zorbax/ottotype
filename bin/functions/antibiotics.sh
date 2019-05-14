#!/bin/bash

antibiotics(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p OUTPUT RESULTS

  wget -q -nc https://raw.githubusercontent.com/CNRDOGM/bin/master/srst2/data/ARGannot.fasta
  docker run --rm -it -v $(pwd):/data -w /data srst2 srst2 --log --output /data/SRST2 --input_pe *fastq.gz \
                        --forward R1 --reverse R2 --gene_db ARGannot.fasta --threads $(nproc) &> /dev/null
  cat SRST2__genes__ARGannot__results.txt | tail -n +2 | sed 's/[?*]//g' | sed -E 's/_S[0-9]{1,}_//' | \
          perl -pe "s/\t\-/#/g; s/\#{1,}//g; s/\_[0-9]{1,}//g; s/\'\'/\-/g" | tail -n+2 | sort \
          > $PWD/RESULTS/srst2_$run_name\_argannot.tsv
  mv SRST2_*genes__ARGannot__results.txt $PWD/OUTPUT

  translate.py $PWD/RESULTS/srst2_$run_name\_argannot.tsv $PWD/RESULTS/antibiotics_$run_name.tsv
  translate.py $PWD/RESULTS/srst2_$run_name\_argannot.tsv $PWD/RESULTS/antibiotics_$run_name\_freq.tsv -freq

  docker run --rm -it -v $(pwd):/data ariba ariba getref card card &> /dev/null && \
  docker run --rm -it -v $(pwd):/data ariba ariba prepareref -f card.fa -m card.tsv card.prepareref &> /dev/null

  mkdir -p ariba\_$run_name
  for i in $(ls *gz | grep -v trim | cut -d\_ -f1,2 | sort | uniq)
  do
    docker run --rm -it -v $(pwd):/data -w /data ariba ariba run /data/card.prepareref \
                             $i\_R1.fastq.gz $i\_R2.fastq.gz $i\_card &> /dev/null && \
    mv $i\_card ariba\_$run_name
  done

  rm -rf card.* *ARGannot*{bam,pileup,bt2,fasta,fai}
}
