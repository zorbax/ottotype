#!/bin/bash

#echo -e "\n#Executing" ${FUNCNAME[0]} "\n"

run_salmonella() {
  run_name=$(basename `pwd` | cut -d\_ -f1)

  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v &> /dev/null
  docker_cmd="docker run --rm -it -v $(pwd):/data -w /data"
  mkdir -p SEQSERO_${run_name} OUTPUT RESULTS

  if [ -z $file ]; then
    file="SCREENING/salm_id.txt"
  fi

  for i in $(cat $file | cut -f1)
  do
    if [[ "find . -type l -name "$i*"" ]]; then
      R1=$(ls -lt $i* | awk '{ print $NF }' | awk '($1 ~ /R1/) { print $1 }')
      R2=$(ls -lt $i* | awk '{ print $NF }' | awk '($1 ~ /R2/) { print $1 }')
      $docker_cmd seqsero SeqSero.py -m 2 -i $R1 $R2 &> /dev/null && \
      mv SeqSero_result* SEQSERO_${run_name}
    else
      $docker_cmd seqsero SeqSero.py -m 2 -i ${i}_R1.fastq.gz ${i}_R2.fastq.gz \
          &> /dev/null && mv SeqSero_result* SEQSERO_${run_name}
    fi
  done

  find SEQSERO_${run_name} -type f -name '*_result.txt' -exec cat {} \
                 > OUTPUT/seqsero_${run_name}_serotype.txt \;

  cat OUTPUT/seqsero_${run_name}_serotype.txt | grep \: | grep -v "^\*\|^Sdf" | \
      grep "serotype\|R1" | sed '/^$/d; s/\:.*Ty/: Ty/; s/_R1.*$//
      s/^.*(s):/Salmonella enterica subsp. enterica serovar/'| tr '\t' ' ' | \
      tr -d '*' | paste - - | cut -d\: -f1 --complement | sed 's/^ //' | sort \
      > RESULTS/seqsero_${run_name}_st01_name.tsv

  cat OUTPUT/seqsero_${run_name}_serotype.txt | grep \: | grep -v "^\*\|^Sdf" | \
      sed 's/H[12].*(//; s/_R1.*$//; s/^.*(s)/Serotipo/;
      s/[)*]//; s/^.*ofile/Perfil antigénico/; s/^O.*:/Antígeno O:/' | \
      tr '\t' ' ' | paste - - - - - - | cut -d\: -f1 --complement | \
      sed 's/^ //' | tr '\t' '\n' > RESULTS/seqsero_${run_name}_st02_antigen.txt
      # N/A values?

  echo "SeqSero: DONE"

  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v &> /dev/null

  mkdir -p SRST2_${run_name}

  $docker_cmd srst2 getmlst.py --species "Salmonella" &> /dev/null

  if [ -f "ARGannot.fasta" ]; then
    "ARGannot is already downloaded"
  else
    repo="https://raw.githubusercontent.com/CNRDOGM"
    wget -q -nc "$repo/srst2/master/data/ARGannot_r3.fasta" -O ARGannot.fasta
  fi

  $docker_cmd srst2 srst2 --log --output /data/SRST2 --input_pe *fastq.gz \
      --forward R1 --reverse R2 --mlst_db Salmonella_enterica.fasta \
      --mlst_definitions senterica.txt --mlst_delimiter '_' \
      --gene_db ARGannot.fasta --threads $(nproc) &> /dev/null

  find . -maxdepth 1 -name "SRST2_*" -type f -not -path "SRST2_$run_name/*" \
      -exec mv {} SRST2_${run_name}/ \;

  cp SRST2_${run_name}/SRST2__compiledResults.txt OUTPUT/srst2_${run_name}_compiledResults.txt

  rm -f senterica.txt *tfa Salmonella* mlst* ARG* *bt2 *fai SRST2.log

  cat OUTPUT/srst2_${run_name}_compiledResults.txt | tail -n+2 | cut -d$'\t' -f 1-9 | \
      sed 's/[?*]//g; /^$/d' | perl -pe 's/_S[0-9]{1,}_//g' | sort \
      > RESULTS/srst2_${run_name}_achtman.tsv

  cat OUTPUT/srst2_${run_name}_compiledResults.txt | tail -n+2 | \
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

  while read line
  do
    id=$(echo $line | cut -d ' ' -f1)
    st=$(echo $line | awk '{ print $2 }')

    if [ ! -f "$HOME/bin/Strain_Senterica.tsv" ]; then
      wget -q -nc -O $HOME/bin/Strain_Senterica.tsv https://git.io/fhpbf
    fi

    cat $HOME/bin/Strain_Senterica.tsv | grep -w ^$st | \
      awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
      awk 'BEGIN {FS="\t"} $10!="" {print}' \
      >> RESULTS/srst2_${run_name}_enterobase.tsv

    cat $HOME/bin/Strain_Senterica.tsv | grep -w ^$st | \
    awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
        awk 'BEGIN {FS="\t"} $10=="" {print}' >> RESULTS/null.tsv
  done < RESULTS/srst2_${run_name}_achtman.tsv

  cat RESULTS/srst2_${run_name}_enterobase.tsv | awk '{ print $1, $2, $11, $10 }' | \
      sed 's/ /\nST:#/; s/ /\neBG:#/; s/ /\nSerotipo:#/; s/#/ /g' \
      > RESULTS/srst2_${run_name}_st02_enterobase.tsv

  id1=$(cat RESULTS/srst2_${run_name}_enterobase.tsv | cut -f1 | sort | uniq)
  id2=$(cat RESULTS/null.tsv | cut -f1 | sort | uniq)
  diff <(echo "$id1") <(echo "$id2") | grep "^\>" \
       > RESULTS/srst2_${run_name}_NF_enterobase.tsv 2>/dev/null

  find RESULTS/ -type f -name "null.tsv" -delete -or -size 0 -delete

  echo "SRST2: DONE"

  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v &> /dev/null
  mkdir -p ARIBA_${run_name}

  $docker_cmd ariba ariba pubmlstget "Salmonella enterica" Salmonella && \
  $docker_cmd ariba ariba getref card card && \
  $docker_cmd ariba ariba prepareref -f card.fa -m card.tsv card.prepareref

  for i in $(cat $file | cut -f1)
  do
    $docker_cmd ariba ariba run /data/Salmonella/ref_db \
                ${i}_R1.fastq.gz ${i}_R2.fastq.gz ${i}_ariba && \
    mv ${i}_ariba ARIBA_${run_name}

    $docker_cmd ariba ariba run /data/card.prepareref \
                ${i}_R1.fastq.gz ${i}_R2.fastq.gz ${i}_card && \
    mv ${i}_card ARIBA_${run_name}
  done

  rm -rf Salmonella card.*
  cd ARIBA_${run_name}

  for i in *ariba
  do
    cd $i
    id=$(echo $i | cut -d\_ -f1)
    if [ -e "mlst_report.tsv" ]; then
      cat mlst_report.tsv | tail -n+2 | awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
      sed 's/[*?]//g' | sort
    fi
    cd ..
  done > ../OUTPUT/ariba_${run_name}_achtman.tsv && cd ..

  cp OUTPUT/ariba_${run_name}_achtman.tsv RESULTS/ariba_${run_name}_achtman.tsv

  while read line
  do
    id=$(echo $line | cut -d ' ' -f1)
    st=$(echo $line | awk '{ print $2 }')

    if [ ! -f "$HOME/bin/Strain_Senterica.tsv" ]; then
      wget -q -nc -O $HOME/bin/Strain_Senterica.tsv https://git.io/fhpbf
    fi

    cat $HOME/bin/Strain_Senterica.tsv | grep -w ^$st | \
      awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
      awk 'BEGIN {FS="\t"} $10!="" {print}' \
      >> RESULTS/ariba_${run_name}_enterobase.tsv

    cat $HOME/bin/Strain_Senterica.tsv | grep -w ^$st | \
    awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
        awk 'BEGIN {FS="\t"} $10=="" {print}' >> RESULTS/null.tsv
  done < OUTPUT/ariba_${run_name}_achtman.tsv

  cat RESULTS/ariba_${run_name}_enterobase.tsv | awk '{ print $1, $2, $11, $10 }' | \
      sed 's/ /\nST:#/; s/ /\neBG:#/; s/ /\nSerotipo:#/; s/#/ /g' \
      > RESULTS/ariba_${run_name}_st02_enterobase.tsv

  id1=$(cat RESULTS/ariba_${run_name}_enterobase.tsv | cut -f1 | sort | uniq)
  id2=$(cat RESULTS/null.tsv | cut -f1 | sort | uniq)
  diff <(echo "$id1") <(echo "$id2") | grep "^\>" \
       > RESULTS/ariba_${run_name}_NF_enterobase.tsv 2>/dev/null

  find RESULTS/ -type f -name "null.tsv" -delete -or -size 0 -delete
  touch -c RESULTS
  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v

  echo "ARIBA: DONE"
}
