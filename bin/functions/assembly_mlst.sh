#!/bin/bash

assembly_mlst(){
  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p ASSEMBLY_MLST_${run_name} OUTPUT RESULTS
  cd ASSEMBLY

  for i in *assembly.fa
  do
    mlst --csv $i --threads $(nproc) \
        >> ../ASSEMBLY_MLST_${run_name}/assembly_mlst_${run_name}.csv 2> /dev/null
  done

  cd ..

  cat ASSEMBLY_MLST_${run_name}/assembly_mlst_${run_name}.csv | \
      sed 's/-assembly.fa//; s/,/\t/g' | sort | \
      tee OUTPUT/assembly_mlst_${run_name}.tsv \
      RESULTS/assembly_mlst_${run_name}.tsv > /dev/null

  cat OUTPUT/assembly_mlst_${run_name}.tsv | sed 's/^/\n>/; s/\t/\tST /2; s/\t/\n/g' \
      > RESULTS/assembly_mlst_${run_name}.list
}
