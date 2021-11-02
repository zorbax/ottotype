#!/bin/bash

kraken_tax(){
    local run_name
    local YGGDRASIL
    run_name=$(basename "$(pwd)" | cut -d\_ -f1)
    mkdir -p KRAKEN2_${run_name} OUTPUT RESULTS

    YGGDRASIL="/mnt/disk2/bin/Kraken2/yggdrasil"

    for r1 in *R1.fastq.gz
    do
        local r2 name
        r2="${r1/R1/R2}"
        name="${r1%%_R1*}"

        echo "$r1"  # --use-mpa-style
        docker run --rm -it -v $YGGDRASIL:/database -v "$(pwd)":/data \
            -u "$(id -u)":"$(id -g)" -w /data kraken2 kraken2 --paired \
            --gzip-compressed --threads "$(nproc)" --db /database \
            --report KRAKEN2_${run_name}/${name}.kraken2-report.tsv \
            ${r1} ${r2} > /dev/null 2> KRAKEN2_${run_name}/${name}.kraken2.log

        tail -2 KRAKEN2_${run_name}/${name}.kraken2.log | paste - - | \
            sed "s/^  /${name}\t/" | sponge KRAKEN2_${run_name}/${name}.kraken2.log
        # rs and sponge > sudo apt-get install moreutils rs
        local tax
        tax=$(awk -F'\t' '{if($1>3) print }' KRAKEN2_${run_name}/${name}.kraken2-report.tsv | \
            grep -P '\t[S]\t'| awk -F'\t' '{ print $6, "#"$1}' | \
            sed -E 's/^ {1,}//; s/[ ]{1,}/_/g' | rs -TeC | sed 's/_#_/ /g')

        echo -e "${name}\t$tax" >> OUTPUT/kraken2_${run_name}.tax.tsv
    done

    cp OUTPUT/kraken2_${run_name}.tax.tsv RESULTS/kraken2_${run_name}.tax.tsv
}
