#!/bin/bash

assembly_mlst(){
    local run_name
    run_name=$(basename "$(pwd)" | cut -d '_' -f1)
    mkdir -p ASSEMBLY_MLST_${run_name} OUTPUT RESULTS
    cd ASSEMBLY || exit

    for i in *assembly.fa
    do
        mlst --csv $i --threads "$(nproc)" \
            >> ../ASSEMBLY_MLST_${run_name}/assembly_mlst_${run_name}.csv 2> /dev/null
    done

    cd ..

    sed 's/-assembly.fa//; s/,/\t/g'  ASSEMBLY_MLST_${run_name}/assembly_mlst_${run_name}.csv | \
        sort | tee OUTPUT/assembly_mlst_${run_name}.tsv \
        RESULTS/assembly_mlst_${run_name}.tsv > /dev/null

   sed 's/^/\n>/; s/\t/\tST /2; s/\t/\n/g' OUTPUT/assembly_mlst_${run_name}.tsv \
        > RESULTS/assembly_mlst_${run_name}.list
}
