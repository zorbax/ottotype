#!/bin/bash

#echo -e "\n#Executing" ${FUNCNAME[0]} "\n"
# TODO: update to SeqSero2
# shellcheck disable=SC2012
run_salmonella() {
    local run_name
    run_name=$(basename "$(pwd)" | cut -d '_' -f1)

    docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v
    local docker_cmd
    docker_cmd="docker run --rm -it -v $(pwd):/data -w /data"
    mkdir -p SEQSERO_${run_name} OUTPUT RESULTS

    if [ -z $file ]; then
        local file
        file="../SCREENING/salm_id.txt"
    fi

    local line
    cut -f1 < $file | while read -r line
    do
        if [[ "find . -type l -name "${line}*"" ]]; then
            local R1 R2
            R1=$(ls -lt ${line}* | awk '{ print $NF }' | awk '($1 ~ /R1/) { print $1 }')
            R2=$(ls -lt ${line}* | awk '{ print $NF }' | awk '($1 ~ /R2/) { print $1 }')
            $docker_cmd seqsero SeqSero.py -m 2 -i $R1 $R2 && \
            mv SeqSero_result* SEQSERO_${run_name}
        else
            $docker_cmd seqsero SeqSero.py -m 2 -i ${line}_R1.fastq.gz ${line}_R2.fastq.gz \
                    && mv SeqSero_result* SEQSERO_${run_name}
        fi
    done

    find SEQSERO_${run_name} -type f -name '*_result.txt' \
            -exec sh -c 'cat $1 > OUTPUT/seqsero_$1_serotype.txt' _ {} \;

    grep ':' OUTPUT/seqsero_${run_name}_serotype.txt | grep -v "^\*\|^Sdf" | \
            grep "serotype\|R1" | sed '/^$/d; s/\:.*Ty/: Ty/; s/_R1.*$//
            s/^.*(s):/Salmonella enterica subsp. enterica serovar/'| tr '\t' ' ' | \
            tr -d '*' | paste - - | cut -d ':' -f1 --complement | sed 's/^ //' | sort \
            > RESULTS/seqsero_${run_name}_st01_name.tsv

    grep ':' OUTPUT/seqsero_${run_name}_serotype.txt | grep -v "^\*\|^Sdf" | \
            sed 's/H[12].*(//; s/_R1.*$//; s/^.*(s)/Serotipo/;
            s/[)*]//; s/^.*ofile/Perfil antigénico/; s/^O.*:/Antígeno O:/' | \
            tr '\t' ' ' | paste - - - - - - | cut -d ':' -f1 --complement | \
            sed 's/^ //' | tr '\t' '\n' > RESULTS/seqsero_${run_name}_st02_antigen.txt
            # N/A values?

    echo "SeqSero: DONE"

    docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v

    mkdir -p SRST2_${run_name}

    $docker_cmd srst2 getmlst.py --species "Salmonella"

    if [ -f "ARGannot.fasta" ]; then
        "ARGannot is already downloaded"
    else
        local repo
        repo="https://raw.githubusercontent.com/CNRDOGM"
        wget -q -nc "$repo/srst2/master/data/ARGannot_r3.fasta" -O ARGannot.fasta
    fi

    $docker_cmd srst2 srst2 --log --output /data/SRST2 --input_pe ./*fastq.gz \
            --forward R1 --reverse R2 --mlst_db Salmonella_enterica.fasta \
            --mlst_definitions senterica.txt --mlst_delimiter '_' \
            --gene_db ARGannot.fasta --threads "$(nproc)"

    find . -maxdepth 1 -name "SRST2_*" -type f -not -path "SRST2_$run_name/*" \
            -exec mv {} SRST2_${run_name}/ \;

    cp SRST2_${run_name}/SRST2__compiledResults.txt OUTPUT/srst2_${run_name}_compiledResults.txt

    rm -f senterica.txt ./*tfa Salmonella* mlst* ARG* ./*bt2 ./*fai SRST2.log

    tail -n+2 OUTPUT/srst2_${run_name}_compiledResults.txt | cut -d$'\t' -f 1-9 | \
            sed 's/[?*]//g; /^$/d' | perl -pe 's/_S[0-9]{1,}_//g' | sort \
            > RESULTS/srst2_${run_name}_achtman.tsv

    tail -n+2 OUTPUT/srst2_${run_name}_compiledResults.txt | \
            awk -v OFS='\t' -v f=2 -v t=13 '{
                for( i=1;i<=NF;i++ )
                if( i>=f&&i<=t )
                    continue;
                else
                    printf( "%s%s",$i,( i!=NF )?OFS:ORS )
            }' | sed 's/[?*]//g' | perl -pe 's/\t\-/#/g; s/\#{1,}//g;
            s/\_[0-9]{1,}//g' | perl -pe 's/_S[0-9]{1,}_//g; s/f{3,}/\tNA/' | \
            sed "s/''/-/" | sort > RESULTS/srst2_${run_name}_argannot.tsv

    translate.py RESULTS/srst2_${run_name}_argannot.tsv RESULTS/antibiotics_$run_name.tsv
    translate.py RESULTS/srst2_${run_name}_argannot.tsv RESULTS/antibiotics_$run_name.freq -freq
    local line2
    while read -r line2
    do
        local id st
        id=$(echo $line2 | cut -d ' ' -f1)
        st=$(echo $line2 | awk '{ print $2 }')

        if [ ! -f "$HOME/bin/Strain_Senterica.tsv" ]; then
            wget -q -nc -O $HOME/bin/Strain_Senterica.tsv https://git.io/fhpbf
        fi

        grep -w ^$st $HOME/bin/Strain_Senterica.tsv | \
            awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
            awk 'BEGIN {FS="\t"} $10!="" {print}' \
            >> RESULTS/srst2_${run_name}_enterobase.tsv

        grep -w ^$st $HOME/bin/Strain_Senterica.tsv | \
        awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
                awk 'BEGIN {FS="\t"} $10=="" {print}' >> RESULTS/null.tsv
    done < RESULTS/srst2_${run_name}_achtman.tsv

    awk '{ print $1, $2, $11, $10 }' RESULTS/srst2_${run_name}_enterobase.tsv | \
            sed 's/ /\nST:#/; s/ /\neBG:#/; s/ /\nSerotipo:#/; s/#/ /g' \
            > RESULTS/srst2_${run_name}_st02_enterobase.tsv

    local id1 id2
    id1=$(cut -f1 RESULTS/srst2_${run_name}_enterobase.tsv | sort | uniq)
    id2=$(cut -f1 RESULTS/null.tsv | sort | uniq)
    diff <(echo "$id1") <(echo "$id2") | grep "^\>" \
              > RESULTS/srst2_${run_name}_NF_enterobase.tsv

    find RESULTS/ -type f -name "null.tsv" -delete -or -size 0 -delete

    echo "SRST2: DONE"

    docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v
    mkdir -p ARIBA_${run_name}

    $docker_cmd ariba ariba pubmlstget "Salmonella enterica" Salmonella && \
    $docker_cmd ariba ariba getref card card && \
    $docker_cmd ariba ariba prepareref -f card.fa -m card.tsv card.prepareref
    local line3
    cut -f1 < $file | while read -r line3
    do
        $docker_cmd ariba ariba run /data/Salmonella/ref_db \
                                ${line3}_R1.fastq.gz ${line3}_R2.fastq.gz ${line3}_ariba && \
        mv ${line3}_ariba ARIBA_${run_name}

        $docker_cmd ariba ariba run /data/card.prepareref \
                                ${line3}_R1.fastq.gz ${line3}_R2.fastq.gz ${line3}_card && \
        mv ${line3}_card ARIBA_${run_name}
    done

    rm -rf Salmonella card.*
    cd ARIBA_${run_name} || exit

    for i in *ariba
    do
        cd $i || exit
        local id
        id=$(echo $i | cut -d '_' -f1)
        if [ -e "mlst_report.tsv" ]; then
            tail -n+2 mlst_report.tsv | \
                awk -v var="$id" -v OFS='\t' '{ print var, $0}' | sed 's/[*?]//g' | sort
        fi
        cd ..
    done > ../OUTPUT/ariba_${run_name}_achtman.tsv && cd ..

    cp OUTPUT/ariba_${run_name}_achtman.tsv RESULTS/ariba_${run_name}_achtman.tsv

    while read -r line
    do
        id=$(echo $line | cut -d ' ' -f1)
        st=$(echo $line | awk '{ print $2 }')

        if [ ! -f "$HOME/bin/Strain_Senterica.tsv" ]; then
            wget -q -nc -O $HOME/bin/Strain_Senterica.tsv https://git.io/fhpbf
        fi

       grep -w ^$st $HOME/bin/Strain_Senterica.tsv | \
            awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
            awk 'BEGIN {FS="\t"} $10!="" {print}' \
            >> RESULTS/ariba_${run_name}_enterobase.tsv

        grep -w ^$st $HOME/bin/Strain_Senterica.tsv | \
        awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
                awk 'BEGIN {FS="\t"} $10=="" {print}' >> RESULTS/null.tsv
    done < OUTPUT/ariba_${run_name}_achtman.tsv

    awk '{ print $1, $2, $11, $10 }' RESULTS/ariba_${run_name}_enterobase.tsv | \
            sed 's/ /\nST:#/; s/ /\neBG:#/; s/ /\nSerotipo:#/; s/#/ /g' \
            > RESULTS/ariba_${run_name}_st02_enterobase.tsv

    id1=$(cut -f1 RESULTS/ariba_${run_name}_enterobase.tsv | sort | uniq)
    id2=$(cut -f1 RESULTS/null.tsv | sort | uniq)
    diff <(echo "$id1") <(echo "$id2") | grep "^\>" \
              > RESULTS/ariba_${run_name}_NF_enterobase.tsv 2>/dev/null

    find RESULTS/ -type f -name "null.tsv" -delete -or -size 0 -delete
    touch -c RESULTS
    docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v

    echo "ARIBA: DONE"
}
