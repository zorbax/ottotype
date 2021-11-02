#!/bin/bash

screen_assembly_tax(){
    local dbNCBI dbRDP dbSILVA
    dbNCBI="$HOME/bin/16S/NCBI.gz"
    dbRDP="$HOME/bin/16S/RDP.gz"
    dbSILVA="$HOME/bin/16S/SILVA.gz"

    for i in ASSEMBLY/*assembly.fa
    do
        DB="$dbNCBI $dbRDP $dbSILVA"
        for k in "${DB[@]}"
        do
            local acc
            acc=$(minimap2 -t "$(nproc)" -x sr ${k} $i 2> /dev/null | awk -F "\t" '$12>0 {print}' | \
                        cut -f 6 | sort | uniq -c | sort -nr | head -n 1 | sed 's/^ *//' | cut -d ' ' -f 2)
            local name
            name=$(basename ${k} .gz)
            if [ -z "$acc" ]; then
                echo "No match found!"
            else
                echo "##########$name##########"
                echo "Reads: $( basename $i -assembly.fa)" | tee -a results_assembly_$name.txt
                echo "Top hit: $acc" | tee -a results_assembly_$name.txt
                local hit binomial
                hit=$(gzip -c -d -f ${k} | grep -m 1 -F $acc)
                echo "Description: $hit" | tee -a results_assembly_$name.txt
                binomial=$(echo $hit | cut -d ' ' -f 2,3)
                echo "Species: $binomial" | tee -a results_assembly_$name.txt
            fi
        done
    done
}
