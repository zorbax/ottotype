#!/bin/bash

kmer_finder(){
    local run_name KmerFinder_DB kraw kgenome
    run_name=$(basename ""$(pwd)"" | cut -d '_' -f1)
    mkdir -p KMERFINDER_${run_name}/{raw,genome}_${run_name} RESULTS/

    KmerFinder_DB=/mnt/disk1/bin/KmerFinder_DB
    kraw="KMERFINDER_$run_name/raw_$run_name"
    kgenome="KMERFINDER_$run_name/genome_$run_name"

    for r1 in *R1.fastq.gz
    do
        local name
        name="${r1%%_R1*}"
        docker run --rm -it -v $KmerFinder_DB:/database -v "$(pwd)":/workdir \
                    -u "$(id -u)":"$(id -g)" -w /data kmerfinder -i /workdir/${r1} \
                    -o /workdir/$kraw/${name} -db /database/bacteria.ATG \
                    -tax /database/bacteria.name -x &> /dev/null
    done

    for i in ASSEMBLY/*assembly.fa
    do
        local genome_name
        genome_name="${i##*/}"; name="${name%%-*}"
        echo $genome_name
        docker run --rm -it -v $KmerFinder_DB:/database -v "$(pwd)":/workdir \
                    -u "$(id -u)":"$(id -g)" -w /data kmerfinder -i /workdir/$i \
                    -o /workdir/$kgenome/$genome_name -db /database/bacteria.ATG \
                    -tax /database/bacteria.name -x &> /dev/null
    done

    for i in *R1.fastq.gz
    do
        name="${i%%_R1*}"
        echo "# RAW_${name}"
        tail -n+2  $kraw/${name}/results.txt | awk -F'\t' -v OFS='\t' '{ print $NF, $3, $13}'
        echo
        #echo "# GENOME_${name}"
        #cat $kgenome/${name}/results.txt | tail -n+2 | awk -F'\t' -v OFS='\t' '{ print $NF, $3, $13}'
        echo
    done > RESULTS/kmerfinder_${run_name}.txt
}
