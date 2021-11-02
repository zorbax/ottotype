#!/bin/bash

antibiotics(){

    local run_name
    run_name=$(basename "$(pwd)" | cut -d\_ -f1)
    mkdir -p OUTPUT RESULTS ANTIBIOTICS_${run_name}

    local repo docker_cmd
    repo="https://raw.githubusercontent.com/CNRDOGM"
    wget -q -nc ${repo}/srst2/master/data/ARGannot_r3.fasta
    docker_cmd="docker run --rm -it -v $(pwd):/data -w /data"

    ${docker_cmd} srst2 srst2 --log --output /data/SRST2 --input_pe ./*fastq.gz \
                  --forward R1 --reverse R2 --gene_db ARGannot_r3.fasta \
                  --threads "$(nproc)" &> /dev/null

    find . -maxdepth 1 -name "SRST2_*" -type f -not -path "ANTIBIOTICS_$run_name/*" \
            -exec mv {} ANTIBIOTICS_${run_name}/ \;

    tail -n +2 ANTIBIOTICS_${run_name}/SRST2__genes__ARGannot_r3__results.txt | \
            sed 's/[?*]//g' | sed -E 's/_S[0-9]{1,}_//' | \
            perl -pe "s/\t\-/#/g; s/\#{1,}//g; s/\_[0-9]{1,}//g; s/\'\'/\-/g;
            s/f{3,}//" | sort > RESULTS/antibiotics_${run_name}_argannot.tsv

    translate.py RESULTS/antibiotics_${run_name}_argannot.tsv \
                              RESULTS/antibiotics_r3_${run_name}.tsv
    translate.py RESULTS/antibiotics_${run_name}_argannot.tsv \
                              RESULTS/antibiotics_r3_${run_name}_freq.tsv -freq

    ${docker_cmd} ariba ariba getref card card &> /dev/null && \
    ${docker_cmd} ariba ariba prepareref -f card.fa -m card.tsv card.prepareref &> /dev/null
    #megares
    for r1 in *R1.fastq.gz
    do
        local r2 name
        r2="${r1/R1/R2}"
        name="${r1%%_R1*}"
        ${docker_cmd} ariba ariba run /data/card.prepareref \
                                ${r1} ${r2} ${name}_card &> /dev/null && \
        mv ${name}_card ANTIBIOTICS_${run_name}
    done

    rm -rf card.* ./*ARGannot*{bt2,fasta,fai} SRST2.log
    ${docker_cmd} ariba ariba summary --preset all out.summary \
                    "$(find ARIBA_PIPELINE/ -type f -name report.tsv -printf "%p ")"
}
