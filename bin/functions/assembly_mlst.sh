#!/bin/bash

assembly_mlst(){
  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p mlst\_$run_name OUTPUT
  cd $PWD/ASSEMBLY
  for i in *assembly.fa
  do
    mlst --csv $i --threads $(nproc) >> ../mlst\_$run_name/assembly_mlst\_$run_name.csv 2> /dev/null
  done
  cd ..
  cat mlst\_$run_name/assembly_mlst\_$run_name.csv | sed 's/-assembly.fa//; s/,/\t/g' | \
        sort  > OUTPUT/assembly_mlst\_$run_name.tsv
  cp OUTPUT/assembly_mlst\_$run_name.tsv RESULTS/assembly_mlst\_$run_name.tsv
}
