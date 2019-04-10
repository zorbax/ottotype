#!/bin/bash

####   _____                 _   _
####  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___
####  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
####  |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
####  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
####

T=$(date +%s)

end(){

  echo -e "\n\n========"
  echo "  DONE"
  echo -e "========\n\n"
  R=$(($(date +%s)-$T))
  D=$((R/60/60/24))
  H=$((R/60/60%24))
  M=$((R/60%60))
  S=$((R%60))

  printf 'CMD-> %s\n' "$0"
  printf 'RUNTIME-> '
  (( $D > 0 )) && printf '%d d ' $D
  (( $H > 0 )) && printf '%d h ' $H
  (( $M > 0 )) && printf '%d m ' $M
  (( $D > 0 || $H > 0 || $M > 0 )) && printf 'and '
  printf '%d s\n' $S
}

tm() {

  local T=$(date +%s)
  local command="$@"
  $@
  rc=$?
  local R=$(($(date +%s)-$T))
  local D=$((R/60/60/24))
  local H=$((R/60/60%24))
  local M=$((R/60%60))
  local S=$((R%60))

  printf 'CMD-> %s\n' "$command"
  printf 'RUNTIME-> '
  (( $D > 0 )) && printf '%d d ' $D
  (( $H > 0 )) && printf '%d h ' $H
  (( $M > 0 )) && printf '%d m ' $M
  (( $D > 0 || $H > 0 || $M > 0 )) && printf 'and '
  printf '%d s\n' $S
  return $rc
}

error(){
  local parent_lineno="$1"
  local script="$2"
  local message="$3"
  local code="${4:-1}"

  if [[ -n "$message" ]] ; then
    echo -e "\n---------------------------------------\n"
    echo -e "ERROR in $script on or near line ${parent_lineno}; exiting with status ${code}"
    echo -e "MESSAGE:\n"
    echo -e "$message"
    echo -e "\n---------------------------------------\n"
  else
    echo -e "\n---------------------------------------\n"
    echo -e "ERROR in $script on or near line ${parent_lineno}; exiting with status ${code}"
    echo -e "\n---------------------------------------\n"
  fi

  exit "${code}"
}

mkdir -p $PWD/log
log_file=log/ottotype.log
if [ -f $log_file ];then
  rm -f $log_file
fi

cat << EOF > $log_file
====================================
Serotyping pipeline from SSB-CNRDOGM
====================================
EOF

echo -e "\nLOG FILE OTTOTYPE" >> $log_file
echo $(date +"%Y-%m-%d %H:%M") >> $log_file

check_path(){

  if [ ! -d "$HOME/bin" ]; then
    mkdir -p $HOME/bin
    cat <<EOF >> $HOME/.bashrc

if [ -d "\$HOME/bin" ] ; then
  export PATH="\$HOME/bin:\$PATH"
fi

EOF
  else
    perl -e 'exit(!(grep(m{^$ENV{HOME}/bin$},split(":", $ENV{PATH}))) > 0)'
    if [ $? != 0 ]; then
      cat <<EOF >> $HOME/.bashrc

if [ -d "\$HOME/bin" ] ; then
  export PATH="\$HOME/bin:\$PATH"
fi

EOF
    fi
  fi

. ~/.bashrc
}

check_connection(){
  echo -e "Network Status\n" &>> $log_file
  adress=$(ip r | grep default | cut -d ' ' -f 3)
  up=$(ping -q -w 1 -c 1 "$adress" > /dev/null && echo ok || echo error)
  if [ "$up" == "ok" ]; then
    echo -e "SESSION_TYPE=local\n" &>> $log_file
  elif [ -n "$SSH_CLIENT" ] || [ -n "$SSH_TTY" ]; then
    echo -e "SESSION_TYPE=remote/ssh\n" &>> $log_file
  else
    echo -e "Network Error, please check your internet connection\n" &>> $log_file
    exit 1
  fi
}

check_dependencies(){
  counter=0
  printf '\n%s\t%20s\n' "DEPENDENCY" "STATUS"
  printf '%s\t%20s\n'   "----------" "------"

  for command in "$@"
  do
    length_command=$(echo $command | wc -m)
    distance_table=$((30 - $length_command))
    distance_expression=$(echo "%${distance_table}s")
    printf '%s' $command
    if ! [ -x "$(which $command 2> /dev/null)" ]; then
      printf $distance_expression
      printf "NOT INSTALLED\n"
      let counter++
    else
      printf $distance_expression
      printf "INSTALLED\n"
    fi
  done

  if [ $counter -gt 0 ]; then
    printf "ERROR: $counter missing dependencies, aborting execution\n" &>> $log_file
    exit 1
  fi
}

check_dockers(){
  counter=0
  printf '\n%s\t%28s\n' "DOCKER" "STATUS"
  printf '%s\t%28s\n' "------" "------"

  for docker in "$@"
  do
    length_docker=$(echo $docker | wc -m)
    distance_table=$((30 - $length_docker))
    distance_expression=$(echo "%${distance_table}s")
    printf '%s' $docker

    if [[ "$(docker images -q $docker 2> /dev/null)" == "" ]]; then
      printf $distance_expression
      printf "NOT INSTALLED\n"
      let counter++
    else
      printf $distance_expression
      printf "INSTALLED\n"
    fi
  done

  if [ $counter -gt 0 ]; then
    printf "ERROR: $counter missing dockers, aborting execution\n" &>> $log_file
    exit 1
  fi
}

check_databases(){

  dbNCBI="$HOME/bin/16S/NCBI.gz"
  dbRDP="$HOME/bin/16S/RDP.gz"
  dbSILVA="$HOME/bin/16S/SILVA.gz"
  dbKraken="/mnt/disk2/bin/Kraken2/yggdrasil"
  #Old version of Kmerfinder database
  # dbKmerfinder_prefix="/mnt/disk1/bin/kmerfinder_DB/bacteria.organisms.ATGAC"
  # dbKmerfinder_list=(.desc.p .p .len.p .ulen.p)
  dbKmerfinder_prefix="/mnt/disk1/bin/KmerFinder_DB/bacteria.ATG."
  dbKmerfinder_list=(seq.b name length.b comp.b)
  dbKmerfinder="${dbKmerfinder_list[@]/#/$dbKmerfinder_prefix}"
  databases="$dbNCBI $dbRDP $dbSILVA $dbKraken $dbKmerfinder"

  counter=0
  printf '\n%s\t%60s\n' "DATABASE" "STATUS"
  printf '%s\t%60s\n'   "--------" "------"

  for db in ${databases}
  do
    length_db=$(echo ${db} | wc -m)
    distance_table=$((70 - $length_db))
    distance_expression=$(echo "%${distance_table}s")
    printf '%s' $db

    if [ ! -e "$db" ] && [ ! -f  "$db" ]; then
      printf $distance_expression
      printf "NOT INSTALLED\n"
      let counter++
    elif [ ! -e "$db" ] && [ ! -d "$db" ]; then
      printf $distance_expression
      printf "NOT INSTALLED\n"
      let counter++
    else
      printf $distance_expression
      printf "INSTALLED\n"
    fi
  done

  if [ $counter -gt 0 ]; then
    printf "ERROR: $counter missing DB, aborting execution\n" &>> $log_file
    echo
    printf "INSTALL DATABASES: cd ottotype/ && bash databases.sh\n" &>> $log_file
    exit 1
  fi
}

checklist(){
  check_path &>> $log_file # || error ${LINENO} $(basename $0)
  check_connection
  check_dependencies salmonella.py minimap2 translate.py trimmomatic \
      sga Rscript idba_ud500 idba_ud spades.py bwa \
      samtools mlst kraken2 bioawk &>> $log_file # || error ${LINENO} $(basename $0)
  check_dockers seqsero srst2 ariba &>> $log_file  # || error ${LINENO} $(basename $0)
  check_databases &>> $log_file  # || error ${LINENO} $(basename $0)
}

clean() {
  find -maxdepth 1 -name "*fastq.gz" -type f -or -type l | rename 's/_L001//; s/_001//;
                          s/_1.fastq/_R1.fastq/ ; s/_2.fastq/_R2.fastq/'
  for i in *fastq.gz
  do
    if  [[ $i == SRR* ]] || [[ $i == ERR* ]] || [[ $i == DRR* ]]; then
      rename "s/_R([12])/_S01_R$1/" $i
    fi
  done
  echo -e "\n\n# The filenames were renamed with the ${FUNCNAME[0]} function" &>> $log_file
}

small_samples() {

  mkdir small_samples && touch small_samples/small_samples_size.txt

  for i in *fastq.gz
  do
    if [[ -L "$i" && -f "$i" ]]; then
      if [[ "$(( $( stat -Lc '%s' $i ) / 1024 / 1024 ))" -le "45" ]]; then
        small_links=`ls -l $i | awk '{ print $NF }'`
        ls -lh $small_links | awk '{ print $NF": " $5}' | \
        sed 's/_L001//; s/_001// ; s/.*\///; s/_R[12].fastq.gz//' \
        >> small_samples/small_samples_size.txt
      fi
    fi
  done

  find . -maxdepth 1 -size -45M -name "*fastq.gz" -exec ls -lh {} \; | \
        awk '{ print $NF": " $5}' | cut -d\/ -f2 | \
        sed 's/_R[12].fastq.gz//' >> small_samples/small_samples_size.txt

  if [ ! -s "small_samples/small_samples_size.txt" ]; then
    rm -rf small_samples
    echo "# Not found small samples in dataset." #&>> $log_file
  else
    cat small_samples/small_samples_size.txt | cut -d\: -f1 | sort | uniq \
                                       > small_samples/small_samples_ids.txt
    n_samples=`cat small_samples/small_samples_ids.txt | wc -l`
    echo -e "\n# This $n_samples sample(s) are too small:\n" #&>> $log_file
    echo -e "ID\tReads" #&>> $log_file

    for i in $(cat small_samples/small_samples_size.txt | cut -d\: -f1 | sort | uniq)
    do
      echo -e "${i}\t$(zcat ${i}*R1.fastq.gz | awk 'NR%4==1' | wc -l )"
      mv -i ${i}_R*.fastq.gz $PWD/small_samples
    done #&>> $log_file
  fi
}

screen_tax() {
  salmonella.py -d . -e _R1.fastq.gz | \
    sed 's/_R[12].fastq.gz//' | tail -n+2 | \
    awk -F'\t' -v OFS='\t' '{
      if($6 == "100")
        { print $1, $6 > "salm_id.txt"}
      else
        { print $1, $6 > "nosalm_id.txt"}
      }'

  dbNCBI="$HOME/bin/16S/NCBI.gz"

  for i in *R1.fastq.gz
  do
    acc=$(minimap2 -t $(nproc) -x sr $dbNCBI $i 2> /dev/null | \
          awk -F "\t" '$12 > 0 { print }' | cut -f 6 | sort | uniq -c | \
          sort -nr | head -n 1 | sed 's/^ *//' | cut -d ' ' -f 2)

    if [ -z "$acc" ]; then
      echo "No match found!"
    else
      echo "Reads: $( echo $i | cut -d\_ -f1,2)"
      echo "Top hit: $acc"
      hit=$(gzip -c -d -f $dbNCBI | grep -m 1 -F $acc)
      echo "Description: $hit"
      echo "Species: $(echo $hit | cut -d ' ' -f2,3)"
    fi
  done > screen_tax_raw_ncbi.txt

  cat screen_tax_raw_ncbi.txt | grep "Reads\|Species" | sed 's/ /#/' | cut -d\# -f2 | \
      grep -B1 Salmonella | sed 's/--//; /^$/d' | paste - - | sort -k1 > salm_id_ncbi.txt 2>/dev/null
  cat screen_tax_raw_ncbi.txt | grep "Reads\|Species" | sed 's/ /#/' | cut -d\# -f2 | \
      sed 's/--//; /^$/d' | paste - - | grep -v Salmonella | sed '/^$/d' | \
      sort -k1 > nosalm_id_ncbi.txt 2>/dev/null

  if [[ -s "salm_id.txt" && -s "salm_id_ncbi.txt" ]]; then
    grep -vwif <(sort -k1 salm_id.txt | cut -f1) <(cut -f1 salm_id_ncbi.txt) \
    > salm-like.txt 2>/dev/null
  else
    if [ -s "salm_id.txt" ]; then
      echo $(sort -k1 salm_id.txt | cut -f1) | tr ' ' '\n' > salm-like.txt 2>/dev/null
    else
      echo $(cut -f1 salm_id_ncbi.txt) | tr ' ' '\n' > salm-like.txt 2>/dev/null
    fi
  fi

  for i in *txt
  do
    if [ ! -s "$i" ]; then
      rm -f $i
    fi
  done
}

run_salmonella() {
  run_name=$(basename `pwd` | cut -d\_ -f1)

  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v &> /dev/null
  mkdir -p SEQSERO_$run_name OUTPUT RESULTS

  if [ -z $file ]; then
    file="salm_id.txt"
  fi

  for i in $(cat $file | cut -f1)
  do
    if [[ "find . -type l -name "$i*"" ]]; then
      R1=`ls -lt $i* | awk '{ print $NF }' | awk '($1 ~ /R1/) { print $1 }'`
      R2=`ls -lt $i* | awk '{ print $NF }' | awk '($1 ~ /R2/) { print $1 }'`
      docker run --rm -it -v $(pwd):/data -w /data seqsero SeqSero.py -m 2 -i $R1 $R2 &> /dev/null && \
      mv SeqSero_result* SEQSERO_$run_name
    else
      docker run --rm -it -v $(pwd):/data -w /data seqsero SeqSero.py -m 2 \
             -i ${i}_R1.fastq.gz ${i}_R2.fastq.gz &> /dev/null && mv SeqSero_result* SEQSERO_$run_name
    fi
  done

  find SEQSERO_$run_name -type f -name '*_result.txt' -exec cat {} \
                 > OUTPUT/seqsero_${run_name}\_serotype.txt \;

  cat OUTPUT/seqsero_${run_name}\_serotype.txt | grep \: | grep -v "^\*\|^Sdf" | \
      grep "serotype\|R1" | sed '/^$/d; s/\:.*Ty/: Ty/; s/_R1.*$//
      s/^.*(s):/Salmonella enterica subsp. enterica serovar/'| tr '\t' ' ' | \
      tr -d '*' | paste - - | cut -d\: -f1 --complement | sed 's/^ //' | sort \
      > RESULTS/seqsero_$run_name\_st01_name.tsv

  cat OUTPUT/seqsero_${run_name}\_serotype.txt | grep \: | grep -v "^\*\|^Sdf" | \
      sed 's/H[12].*(//; s/_R1.*$//; s/^.*(s)/Serotipo/;
      s/[)*]//; s/^.*ofile/Perfil antigénico/; s/^O.*:/Antígeno O:/' | \
      tr '\t' ' ' | paste - - - - - - | cut -d\: -f1 --complement | \
      sed 's/^ //' | tr '\t' '\n' > RESULTS/seqsero_$run_name\_st02_antigen.txt
      # N/A values?

  echo "SeqSero: DONE"

  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v &> /dev/null

  mkdir -p SRST2_$run_name

  docker run --rm -it -v $(pwd):/data -w /data srst2 getmlst.py --species "Salmonella" &> /dev/null

  if [ -f "ARGannot.fasta" ]; then
    "ARGannot is already downloaded"
  else
    repo="https://raw.githubusercontent.com/CNRDOGM"
    wget -q -nc "$repo/srst2/master/data/ARGannot_r3.fasta" -O ARGannot.fasta
  fi

  docker run --rm -it -v $(pwd):/data -w /data srst2 srst2 --log --output /data/SRST2 \
         --input_pe *fastq.gz --forward R1 --reverse R2 --mlst_db Salmonella_enterica.fasta \
         --mlst_definitions senterica.txt --mlst_delimiter '_' --gene_db ARGannot.fasta \
         --threads $(nproc) &> /dev/null

  find . -maxdepth 1 -name "SRST2_*" -type f -not -path "SRST2_$run_name/*" -exec mv {} SRST2_$run_name/ \;

  cp SRST2_$run_name/SRST2__compiledResults.txt OUTPUT/srst2_$run_name\_compiledResults.txt

  rm -f senterica.txt *tfa Salmonella* mlst* ARG* *bt2 *fai SRST2.log

  cat OUTPUT/srst2_$run_name\_compiledResults.txt | tail -n+2 | cut -d$'\t' -f 1-9 | \
      sed 's/[?*]//g; /^$/d' | perl -pe 's/_S[0-9]{1,}_//g' | sort \
      > RESULTS/srst2_$run_name\_achtman.tsv

  cat OUTPUT/srst2_$run_name\_compiledResults.txt | tail -n+2 | \
      awk -v OFS='\t' -v f=2 -v t=13 '{
        for( i=1;i<=NF;i++ )
        if( i>=f&&i<=t )
          continue;
        else
          printf( "%s%s",$i,( i!=NF )?OFS:ORS )
      }' | sed 's/[?*]//g' | perl -pe 's/\t\-/#/g; s/\#{1,}//g;
      s/\_[0-9]{1,}//g' | perl -pe 's/_S[0-9]{1,}_//g; s/f{3,}/\tNA/' | \
      sed "s/''/-/" | sort > RESULTS/srst2_$run_name\_argannot.tsv

  translate.py RESULTS/srst2_$run_name\_argannot.tsv RESULTS/antibiotics_$run_name.tsv
  translate.py RESULTS/srst2_$run_name\_argannot.tsv RESULTS/antibiotics_$run_name.freq -freq

  while read line
  do
    id=`echo $line | cut -d ' ' -f1`
    st=`echo $line | awk '{ print $2 }'`

    if [ ! -f "$HOME/bin/Strain_Senterica.tsv" ]; then
      wget -q -nc -O $HOME/bin/Strain_Senterica.tsv https://git.io/fhpbf
    fi

    cat $HOME/bin/Strain_Senterica.tsv | grep -w ^$st | \
      awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
      awk 'BEGIN {FS="\t"} $10!="" {print}' \
      >> RESULTS/srst2_${run_name}\_enterobase.tsv

    cat $HOME/bin/Strain_Senterica.tsv | grep -w ^$st | \
    awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
        awk 'BEGIN {FS="\t"} $10=="" {print}' >> RESULTS/null.tsv
  done < RESULTS/srst2_$run_name\_achtman.tsv

  cat RESULTS/srst2_${run_name}\_enterobase.tsv | awk '{ print $1, $2, $11, $10 }' | \
      sed 's/ /\nST:#/; s/ /\neBG:#/; s/ /\nSerotipo:#/; s/#/ /g' \
      > RESULTS/srst2_${run_name}\_st02_enterobase.tsv

  id1=`cat RESULTS/srst2_${run_name}\_enterobase.tsv | cut -f1 | sort | uniq`
  id2=`cat RESULTS/null.tsv | cut -f1 | sort | uniq`
  diff <(echo "$id1") <(echo "$id2") | grep "^\>" \
       > RESULTS/srst2_${run_name}\_NF_enterobase.tsv 2>/dev/null

  find RESULTS/ -type f -name "null.tsv" -delete -or -size 0 -delete

  echo "SRST2: DONE"

  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v
  mkdir ARIBA_$run_name

  docker run --rm -it -v $(pwd):/data ariba ariba pubmlstget "Salmonella enterica" Salmonella &> /dev/null && \
  docker run --rm -it -v $(pwd):/data ariba ariba getref card card &> /dev/null && \
  docker run --rm -it -v $(pwd):/data ariba ariba prepareref -f card.fa -m card.tsv card.prepareref &> /dev/null

  if [ -z $file ]; then
    file="salm_id.txt"
    for i in $(cat $file | cut -f1)
    do
      docker run --rm -it -v $(pwd):/data -w /data ariba ariba run /data/Salmonella/ref_db \
                                ${i}_R1.fastq.gz ${i}_R2.fastq.gz ${i}_ariba &> /dev/null && \
      mv ${i}_ariba ARIBA_${run_name}

      docker run --rm -it -v $(pwd):/data -w /data ariba ariba run /data/card.prepareref \
                               ${i}_R1.fastq.gz ${i}_R2.fastq.gz ${i}_card &> /dev/null && \
      mv ${i}_card ARIBA_${run_name}
    done
  fi

  rm -rf Salmonella card.*

  cd ARIBA_$run_name
  for i in *ariba
  do
    cd $i
    id=`echo $i | cut -d\_ -f1`
    if [ -e "mlst_report.tsv" ]; then
      cat mlst_report.tsv | tail -n+2 | awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
      sed 's/[*?]//g' | sort
      cd ..
    else
      cd ..
    fi
  done > ../OUTPUT/ariba_$run_name\_achtman.tsv && cd ..

  cp OUTPUT/ariba_$run_name\_achtman.tsv RESULTS/ariba_$run_name\_achtman.tsv

  while read line
  do
    id=`echo $line | cut -d ' ' -f1`
    st=`echo $line | awk '{ print $2 }'`

    if [ ! -f "$HOME/bin/Strain_Senterica.tsv" ]; then
      wget -q -nc -O $HOME/bin/Strain_Senterica.tsv https://git.io/fhpbf
    fi

    cat $HOME/bin/Strain_Senterica.tsv | grep -w ^$st | \
      awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
      awk 'BEGIN {FS="\t"} $10!="" {print}' \
      >> RESULTS/ariba_$run_name\_enterobase.tsv

    cat $HOME/bin/Strain_Senterica.tsv | grep -w ^$st | \
    awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
        awk 'BEGIN {FS="\t"} $10=="" {print}' >> RESULTS/null.tsv
  done < OUTPUT/ariba_$run_name\_achtman.tsv

  cat RESULTS/ariba_$run_name\_enterobase.tsv | awk '{ print $1, $2, $11, $10 }' | \
      sed 's/ /\nST:#/; s/ /\neBG:#/; s/ /\nSerotipo:#/; s/#/ /g' \
      > RESULTS/ariba_$run_name\_st02_enterobase.tsv

  id1=`cat RESULTS/ariba_${run_name}\_enterobase.tsv | cut -f1 | sort | uniq`
  id2=`cat RESULTS/null.tsv | cut -f1 | sort | uniq`
  diff <(echo "$id1") <(echo "$id2") | grep "^\>" \
       > RESULTS/ariba_$run_name\_NF_enterobase.tsv 2>/dev/null

  find RESULTS/ -type f -name "null.tsv" -delete -or -size 0 -delete
  touch RESULTS
  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v

  echo "ARIBA: DONE"

}

trimming() {

  mkdir -p $PWD/TRIMMING/1U2U
  for i in $(ls *fastq.gz | cut -d\_ -f1,2 | sort | uniq )
  do
    trimmomatic PE -phred33 -threads $(nproc) ${i}_R1.fastq.gz ${i}_R2.fastq.gz \
              TRIMMING/${i}_R1.trim.fastq.gz TRIMMING/1U2U/$i.1U.trim.fastq.gz \
              TRIMMING/${i}_R2.trim.fastq.gz TRIMMING/1U2U/$i.2U.trim.fastq.gz \
              SLIDINGWINDOW:4:20 MINLEN:70 &> $i.trim.log

    if [ $? -eq 0 ]; then
        rm $i.trim.log
    fi
  done
}

assembly() {

  for i in $(ls $PWD/TRIMMING/*gz | cut -d\_ -f1,2 | sort | uniq)
  do
    mkdir -p $PWD/ASSEMBLY

    name=$(basename $i)
    echo "${name}"
    sga preprocess -q 25 -f 20 --pe-mode=1  ${i}_R1.trim.fastq.gz ${i}_R2.trim.fastq.gz > ${name}_12.t.pp.fq 2> /dev/null
    median=`cat ${name}_12.t.pp.fq | awk 'NR%4==2 { print length }' | head -20000 | \
            Rscript -e 'summary (as.numeric (readLines ("stdin")))' | tail -n+2 | \
            awk '{ print $3 }'  | cut -d\. -f1`

    if [ $median -le 200 ]; then
      #mem=`cat /proc/meminfo | grep MemTotal | awk '{ print $2}'`
      #kb=$(echo "($mem*0.5)" | bc | cut -d\. -f1)
      sga index -t $(nproc) -a ropebwt ${name}_12.t.pp.fq > ${name}.index.out 2> ${name}.index.err
      sga correct -t $(nproc) ${name}_12.t.pp.fq -o ${name}_12.t.pp.ec.fq > ${name}.correct.out 2> ${name}.correct.err
    else
      sga index -t $(nproc) -a sais ${name}_12.t.pp.fq > ${name}.index.out 2> ${name}.index.err
      sga correct -t $(nproc) ${name}_12.t.pp.fq -o ${name}_12.t.pp.ec.fq > ${name}.correct.out 2> ${name}.correct.err
    fi

    sed -n '1~4s/^@/>/p;2~4p' ${name}_12.t.pp.ec.fq > ${name}_12.t.pp.ec.fa
    medianec=`cat ${name}_12.t.pp.ec.fa | awk '$0 !~ /^>/ { print length }' | head -20000 | \
              Rscript -e 'summary (as.numeric (readLines ("stdin")))' | tail -n+2 | \
              awk '{ print $3}'  | cut -d\. -f1`

    if [[ $medianec -ge 128 ]]; then
      idba_ud500 -r ${name}_12.t.pp.ec.fa -o ${name} --mink 35 --maxk 249 --num_threads $(nproc) --min_pairs 2 > /dev/null
    else
      idba_ud -r ${name}_12.t.pp.ec.fa -o ${name} --mink 35 --maxk 124 --num_threads $(nproc) --min_pairs 2 > /dev/null
    fi

    if [ $? != 0 ]; then
       cp ${name}/contig.fa ${name}-idba-assembly.fa
    else
       cp ${name}/scaffold.fa ${name}-idba-assembly.fa
    fi

    rm -rf *pp.{fq,bwt,sai,rbwt} *out *err *rsai *ec.{fq,fa} ${name}

    cat ${name}-idba-assembly.fa | sed ':a;N;/^>/M!s/\n//;ta;P;D' | \
         awk '/^>/ { getline seq } length(seq) >500 { print $0 "\n" seq }' > $PWD/ASSEMBLY/${name}-idba-assembly.fa
    rm ${name}-idba-assembly.fa
  done
}

assembly_spades() {
  mkdir -p ASSEMBLY
  memory=`awk '{ printf "%.2f", $2/1024/1024 ; exit}' /proc/meminfo | cut -d\. -f1`
  for i in $(ls TRIMMING/*gz | cut -d\_ -f1,2 | sort | uniq)
  do
    name=`echo $i | cut -d\/ -f2`
    spades.py --pe1-1 TRIMMING/${name}_R1.trim.fastq.gz --pe1-2 TRIMMING/${name}_R2.trim.fastq.gz \
        --pe1-s TRIMMING/1U2U/${name}.1U.trim.fastq.gz --pe1-s TRIMMING/1U2U/${name}.2U.trim.fastq.gz \
        -o ${name}_spades -t $(nproc) -m $memory &>/dev/null
  # --careful --only-assembler using sga corrected reads
  find ${name}_spades -maxdepth 2 -type f -name 'scaffolds.fasta' -exec cp {} ${name}.tmp \;
  cat ${name}.tmp | sed ':a;N;/^>/M!s/\n//;ta;P;D' | \
       awk '/^>/ { getline seq } length(seq) >500 { print $0 "\n" seq }' > ASSEMBLY/${name}-spades-assembly.fa
  rm -rf ${name}_spades ${name}.tmp
done
}

assembly_stats_cov() {

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p OUTPUT $PWD/ASSEMBLY/Stats RESULTS
  echo -e "Assembly\tNumber_of_contigs\tCoverage\tAssembly_size\tLargest_contig\tN50\tN90" > assembly_$run_name.stats

  for i in $PWD/ASSEMBLY/*assembly.fa
  do
    assembly=ASSEMBLY/$(basename "$i" | cut -d\- -f3 --complement)
    name=`echo $assembly | cut -d\/ -f2`
    reads=TRIMMING/$(basename $i | cut -d\- -f1)

    echo -e "Assembly:\t${name}" | tee ${name}.stats
    bwa index -a bwtsw $i -p ${name} &> /dev/null
    bwa mem -t $(nproc) ${name} $reads\_R1.trim.fastq.gz $reads\_R2.trim.fastq.gz \
    2> /dev/null | samtools view -Sb -F4 -@ $(nproc) - | \
    samtools sort -@ $(nproc) - -o ${name}.mapped.sorted.bam 2>/dev/null
    samtools index -@ $(nproc) ${name}.mapped.sorted.bam
    rm ${name}.{bwt,pac,ann,amb,sa}

    cov=$(samtools mpileup ${name}.mapped.sorted.bam 2> /dev/null | \
      awk '{ count++ ; sum += $4 } END { printf "%s\t%.2f\n", "Coverage:", sum/count }')
    cover=`echo $cov | awk '{ print $2 }'`

    cat $i | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }'| \
                                 sed '/^>/ d'| awk '{ print length($0) }' | sort -gr > ${name}_contig_lengths.stat

    contigs_number=$(cat ${name}_contig_lengths.stat | wc -l)
    assembly_size=$(paste -sd+ ${name}_contig_lengths.stat | bc)
    awk 'BEGIN {sum=0} {sum= sum+$0; print sum}' ${name}_contig_lengths.stat > ${name}_contig_lengths_cum.stat

    awk -v var=$assembly_size 'BEGIN {FS=OFS=","} {for (i=1; i<=NF; i++) $i/=var;}1' \
                                    ${name}_contig_lengths_cum.stat > ${name}_cum_perc.stat

    paste ${name}_contig_lengths.stat ${name}_cum_perc.stat > ${name}_matrix.stat

    N50=$(awk '$2 >= 0.50' ${name}_matrix.stat | head -1 | awk '{ print $1 }')
    N90=$(awk '$2 >= 0.90' ${name}_matrix.stat | head -1 | awk '{ print $1 }')
    large_contig=$(head -1 ${name}_contig_lengths.stat)

    rm *stat
    echo -e "Number of contigs:\t$contigs_number\n$cov" >> ${name}.stats
    echo -e "Assembly size:\t$assembly_size" >> ${name}.stats
    echo -e "Largest contig:\t$large_contig" >> ${name}.stats
    echo -e "N50:\t$N50\nN90:\t$N90" >> ${name}.stats
    mv ${name}.stats ASSEMBLY/Stats
    echo -e "${name}\t$contigs_number\t$cover\t$assembly_size\t$large_contig\t$N50\t$N90" >> assembly_$run_name.stats
    find . -maxdepth 1 -type f \( -name "*.bam" -o -name "*.bai" \) -exec rm {} \;
  done
  mv assembly_$run_name.stats $PWD/OUTPUT && cp $PWD/OUTPUT/assembly_$run_name.stats $PWD/RESULTS
}

assembly_mlst(){
  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p ASSEMBLY_MLST_${run_name} OUTPUT RESULTS
  cd ASSEMBLY
  for i in *assembly.fa
  do
    mlst --csv $i --threads $(nproc) >> ../ASSEMBLY_MLST_${run_name}/assembly_mlst_${run_name}.csv 2> /dev/null
  done
  cd ..
  cat ASSEMBLY_MLST_${run_name}/assembly_mlst_${run_name}.csv | sed 's/-assembly.fa//; s/,/\t/g' | \
        sort | tee OUTPUT/assembly_mlst_${run_name}.tsv RESULTS/assembly_mlst_${run_name}.tsv > /dev/null

  cat OUTPUT/assembly_mlst_${run_name}.tsv | sed 's/^/\n>/; s/\t/\tST /2; s/\t/\n/g' \
      > RESULTS/assembly_mlst_${run_name}.list
}

kmer_finder(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p KMERFINDER_${run_name}/{raw,genome}_${run_name} RESULTS/

  KmerFinder_DB=/mnt/disk1/bin/KmerFinder_DB
  kraw="KMERFINDER_$run_name/raw_$run_name"
  kgenome="KMERFINDER_$run_name/genome_$run_name"

  for i in $(ls *fastq.gz | cut -d\_ -f1,2 | sort | uniq)
  do
    docker run --rm -it -v $KmerFinder_DB:/database -v $(pwd):/workdir \
          -u $(id -u):$(id -g) -w /data kmerfinder -i /workdir/${i}_R1.fastq.gz \
          -o /workdir/$kraw/$i -db /database/bacteria.ATG \
          -tax /database/bacteria.name -x &> /dev/null
  done

  for i in ASSEMBLY/*assembly.fa
  do
    genome_name=`basename $i | cut -d\- -f1 | cut -d\/ -f2`
    docker run --rm -it -v $KmerFinder_DB:/database -v $(pwd):/workdir \
          -u $(id -u):$(id -g) -w /data kmerfinder -i /workdir/$i \
          -o /workdir/$kgenome/$genome_name -db /database/bacteria.ATG \
          -tax /database/bacteria.name -x &> /dev/null
  done

  for i in $(ls *fastq.gz | cut -d\_ -f1,2 | sort | uniq)
  do
    name=`echo $i | cut -d\_ -f1`
    echo "# RAW_${name}"
    cat $kraw/$i/results.txt | tail -n+2 | awk -F'\t' -v OFS='\t' '{ print $NF, $3, $13}'
    echo
    echo "# GENOME_${name}"
    cat $kgenome/$i/results.txt | tail -n+2 | awk -F'\t' -v OFS='\t' '{ print $NF, $3, $13}'
    echo
  done > RESULTS/kmerfinder_${run_name}.txt
}

kraken_tax(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p KRAKEN2_${run_name} OUTPUT RESULTS


  YGGDRASIL=/mnt/disk2/bin/Kraken2/yggdrasil
: <<'END'
  docker run --rm -it -v $YGGDRASIL:/database -v $(pwd):/workdir \
        -u $(id -u):$(id -g) -w /data kmerfinder -i /workdir/${i}_R1.fastq.gz \
        -o /workdir/$kraw/$i -db /database/bacteria.ATG \
        -tax /database/bacteria.name -x &> /dev/null

  docker run --rm -it -v $(pwd):/data -u $(id -u):$(id -g) -w /data kraken2 kraken2 --help
END

  for i in $(ls *fastq.gz | cut -d\_ -f1,2 | sort | uniq)
  do
    echo "$i"  # --use-mpa-style
    kraken2 --paired --gzip-compressed --threads $(nproc) \
      --db $YGGDRASIL --report KRAKEN2_${run_name}/$i.kraken2-report.tsv \
      ${i}_R1.fastq.gz ${i}_R2.fastq.gz > /dev/null 2> KRAKEN2_${run_name}/$i.kraken2.log

    cat KRAKEN2_${run_name}/$i.kraken2.log | tail -2 | paste - - | sed "s/^  /$i\t/" | \
                                               sponge KRAKEN2_${run_name}/$i.kraken2.log
    # sponge > sudo apt-get install moreutils
    # Test output format
    tax=`cat KRAKEN2_${run_name}/$i.kraken2-report.tsv | awk -F'\t' '{if($1>5) print }' | \
         grep -P '\t[DPCOFGS]\t' | sed '0,/D/s//K/' | awk -F'\t' '$4=tolower($4){ print $4"_", $6}' | \
         sed -E 's/[ ]{1,}/_/g' | tr  "\n" ";" | sed 's/;$/\n/'`
    echo -e "$i\t$tax" > KRAKEN2_${run_name}/$i.tax.tsv

  done
 #Test
  cat KRAKEN2_${run_name}/*tax.tsv | awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1'| \
      sed -e 's/k__/#/; s/s__/#/; s/\(#\).*\(#\)//' | sed -E 's/_S[0-9]{1,}//' | \
      awk -F'\t' -v OFS='\t' '{gsub(";s_"," |",$2);gsub("_"," ",$2)}1' | \
      perl -pe 'if(/\#/){s/\ /\_/g}' > OUTPUT/kraken_${run_name}.tax.tsv
  #sed 's/;p  /#/; s/|/#|/; s/\(#\).*\(#\)/ /'
  cp OUTPUT/kraken_${run_name}.tax.tsv RESULTS/kraken_${run_name}.tax.tsv
}

antibiotics(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p OUTPUT RESULTS ANTIBIOTICS_$run_name

  wget -q -nc https://raw.githubusercontent.com/CNRDOGM/srst2/master/data/ARGannot_r3.fasta
  docker run --rm -it -v $(pwd):/data -w /data srst2 srst2 --log --output /data/SRST2 --input_pe *fastq.gz \
                        --forward R1 --reverse R2 --gene_db ARGannot_r3.fasta --threads $(nproc) &> /dev/null

  find . -maxdepth 1 -name "SRST2_*" -type f -not -path "ANTIBIOTICS_$run_name/*" -exec mv {} ANTIBIOTICS_$run_name/ \;

  cat ANTIBIOTICS_$run_name/SRST2__genes__ARGannot_r3__results.txt | tail -n +2 | \
      sed 's/[?*]//g' | sed -E 's/_S[0-9]{1,}_//' | perl -pe "s/\t\-/#/g; s/\#{1,}//g;
      s/\_[0-9]{1,}//g; s/\'\'/\-/g; s/f{3,}//" | sort > RESULTS/antibiotics_$run_name\_argannot.tsv

  translate.py RESULTS/antibiotics_$run_name\_argannot.tsv RESULTS/antibiotics_r3_$run_name.tsv
  translate.py RESULTS/antibiotics_$run_name\_argannot.tsv RESULTS/antibiotics_r3_$run_name\_freq.tsv -freq

  docker run --rm -it -v $(pwd):/data ariba ariba getref card card &> /dev/null && \
  docker run --rm -it -v $(pwd):/data ariba ariba prepareref -f card.fa -m card.tsv card.prepareref &> /dev/null

  for i in $(ls *gz | cut -d\_ -f1,2 | sort | uniq)
  do
    docker run --rm -it -v $(pwd):/data -w /data ariba ariba run /data/card.prepareref \
                             ${i}_R1.fastq.gz ${i}_R2.fastq.gz ${i}_card &> /dev/null && \
    mv ${i}_card ANTIBIOTICS_$run_name
  done

  rm -rf card.* *ARGannot*{bt2,fasta,fai} SRST2.log
  docker run --rm -it -v $(pwd):/data -w /data ariba ariba summary --preset all out.summary \
                    $(find ARIBA_PIPELINE/ -type f -name report.tsv -printf "%p ")
  # Filter out.summary from CARD
}

plasmids(){
  run_name=$(basename `pwd` | cut -d\_ -f1)
  memory=`awk '{ printf "%.2f", $2/1024 ; exit}' /proc/meminfo | cut -d\. -f1`

  ln /mnt/disk1/bin/plasmidid_db/plasmid.salmonella.nr100.fna .

  for i in $(ls *gz | grep -v trim | cut -d\_ -f1,2 | sort | uniq)
  do
    cp ASSEMBLY/$i-idba-assembly.fa $i.fna
    docker run --rm -it -v $(pwd):/data -w /data \
          buisciii/plasmidid plasmidID.sh \
          -1 ${i}_R1.fastq.gz -2 ${i}_R2.fastq.gz -T $(nproc)\
          -d plasmid.salmonella.nr100.fna -M $memory\
          -c $i.fna --no-trim -s $i -g plasmids_${run_name}
    rm $i.fna
  done
  rm plasmid.complete.nr100.fna
}



####   __  __       _
####  |  \/  | __ _(_)_ __
####  | |\/| |/ _` | | '_ \
####  | |  | | (_| | | | | |
####  |_|  |_|\__,_|_|_| |_|
####

docker images --no-trunc | grep '<none>' | awk '{ print $3 }' | xargs -r docker rmi
run_name=$(basename `pwd` | cut -d\_ -f1)

echo "Clean"
clean
echo "Checklist"
checklist
echo "Small samples"
small_samples
echo "Taxa screening"
screen_tax

: <<'END'
echo "SALMONELLA"
if [ -s "salm_id.txt" ]; then
  mkdir -p SALMONELLA

  while read -r fastq
  do
    find . -name "$fastq*fastq.gz" -type f -not -path "$PWD/SALMONELLA/*" -print0 | \
           xargs -0 mv -t "$PWD/SALMONELLA" 2>/dev/null
  done < <(cat salm_id.txt | cut -f1)

  cd SALMONELLA

  file="../salm_id.txt"
  echo "run_salmonella"
  run_salmonella #&>> $log_file || error ${LINENO} $(basename $0)
  echo "Trimming"
  trimming
  echo "Assembly"
  assembly
  echo "Assembly stats"
  assembly_stats_cov
  echo "Plasmids"
  plasmids
  cd ..
fi

echo "SALM-LIKE"

if [ -s "salm-like.txt" ]; then
  mkdir -p SALM-LIKE

  while read -r fastq
  do
    find . -name "$fastq*fastq.gz" -type f -not -path "$PWD/SALM-LIKE/*" -print0 | \
           xargs -0 mv -t "$PWD/SALM-LIKE" 2>/dev/null
  done < <(cat salm-like.txt | cut -f1)

  cd SALM-LIKE

  file="../salm-like.txt"
  echo "run_salmonella"
  run_salmonella
  echo "trimm"
  trimming
  echo "assembly"
  assembly
  echo "stats"
  assembly_stats_cov
  echo "Assembly mlst"
  assembly_mlst
  echo "kmerfinder"
  kmer_finder
  echo "Kraken"
  kraken_tax
  cd ..
fi

echo "OTHERS"

if [ -s "nosalm_id_ncbi.txt" ]; then
  file="nosalm_id_ncbi.txt"
else
  file="nosalm_id.txt"
fi
  mkdir -p OTHERS

  while read -r fastq
  do
    find . -name "$fastq*fastq.gz" -type f -not -path "$PWD/OTHERS/*" -print0 | \
           xargs -0 mv -t "$PWD/OTHERS" 2>/dev/null
  done < <(cat $file | cut -f1)

  cd OTHERS

  echo "Trimming"
  trimming
  echo "Assembly"
  assembly
  echo "Assembly stats"
  assembly_stats_cov
  echo "Assembly mlst"
  assembly_mlst
  echo "Kmerfinder"
  kmer_finder
  echo "Kraken"
  kraken_tax
  echo "Antibiotics"
  antibiotics
  cd ..

fi
END

run_name=$(basename `pwd` | cut -d\_ -f1)
path_results=`find . -type d -name RESULTS`
mkdir RESULTS_${run_name}

for i in $path_results
do
  cp $i/* RESULTS_${run_name}
done

end &>> $log_file
