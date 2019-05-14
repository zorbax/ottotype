#!/bin/bash

####   _____                 _   _
####  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___
####  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
####  |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
####  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
####

source bin/functions/log.sh
source bin/functions/check.sh
source bin/functions/clean.sh
source bin/functions/small_samples.sh
source bin/functions/screen_raw_tax.sh
source bin/functions/salmonella_typing.sh
source bin/functions/ecoli_typing.sh
source bin/functions/trimming.sh
source bin/functions/assembly_idba.sh
source bin/functions/assembly_spades.sh


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

mkdir -p log
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



kraken_tax(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p KRAKEN2_${run_name} OUTPUT RESULTS

  YGGDRASIL="/mnt/disk2/bin/Kraken2/yggdrasil"

  for r1 in *R1.fastq.gz
  do
    r2="${r1/R1/R2}"
    name="${r1%%_R1*}"

    echo "$i"  # --use-mpa-style
    docker run --rm -it -v $YGGDRASIL:/database -v $(pwd):/data \
         -u $(id -u):$(id -g) -w /data kraken2 kraken2 --paired \
         --gzip-compressed --threads $(nproc) --db /database \
         --report KRAKEN2_${run_name}/${name}.kraken2-report.tsv \
         ${r1} ${r2} > /dev/null 2> KRAKEN2_${run_name}/${name}.kraken2.log

    cat KRAKEN2_${run_name}/${name}.kraken2.log | tail -2 | paste - - | \
        sed "s/^  /${name}\t/" | sponge KRAKEN2_${run_name}/${name}.kraken2.log
    # rs and sponge > sudo apt-get install moreutils rs
    tax=$(cat KRAKEN2_${run_name}/${name}.kraken2-report.tsv | \
          awk -F'\t' '{if($1>5) print }' | grep -P '\t[S]\t'| \
          awk -F'\t' '{ print $6, "#"$1}' | sed -E 's/^ {1,}//; s/[ ]{1,}/_/g' | \
          rs -TeC | sed 's/_#_/ /g')

    echo -e "${name}\t$tax" >> OUTPUT/kraken2_${run_name}.tax.tsv
  done

  cp OUTPUT/kraken2_${run_name}.tax.tsv RESULTS/kraken2_${run_name}.tax.tsv
}

antibiotics(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p OUTPUT RESULTS ANTIBIOTICS_$run_name

  repo="https://raw.githubusercontent.com/CNRDOGM"
  wget -q -nc $repo/srst2/master/data/ARGannot_r3.fasta
  docker_cmd="docker run --rm -it -v $(pwd):/data -w /data"

  $docker_cmd srst2 srst2 --log --output /data/SRST2 --input_pe *fastq.gz \
     --forward R1 --reverse R2 --gene_db ARGannot_r3.fasta --threads $(nproc) &> /dev/null

  find . -maxdepth 1 -name "SRST2_*" -type f -not -path "ANTIBIOTICS_$run_name/*" \
      -exec mv {} ANTIBIOTICS_$run_name/ \;

  cat ANTIBIOTICS_$run_name/SRST2__genes__ARGannot_r3__results.txt | \
      tail -n +2 | sed 's/[?*]//g' | sed -E 's/_S[0-9]{1,}_//' | \
      perl -pe "s/\t\-/#/g; s/\#{1,}//g; s/\_[0-9]{1,}//g; s/\'\'/\-/g;
      s/f{3,}//" | sort > RESULTS/antibiotics_${run_name}_argannot.tsv

  translate.py RESULTS/antibiotics_${run_name}_argannot.tsv RESULTS/antibiotics_r3_${run_name}.tsv
  translate.py RESULTS/antibiotics_${run_name}_argannot.tsv RESULTS/antibiotics_r3_${run_name}_freq.tsv -freq

  $docker_cmd ariba ariba getref card card &> /dev/null && \
  $docker_cmd ariba ariba prepareref -f card.fa -m card.tsv card.prepareref &> /dev/null

  for r1 in *R1.fastq.gz
  do
    r2="${r1/R1/R2}"
    name="${r1%%_R1*}"
    $docker_cmd ariba ariba run /data/card.prepareref \
                ${r1} ${r2} ${name}_card &> /dev/null && \
    mv ${name}_card ANTIBIOTICS_$run_name
  done

  rm -rf card.* *ARGannot*{bt2,fasta,fai} SRST2.log
  $docker_cmd ariba ariba summary --preset all out.summary \
          $(find ARIBA_PIPELINE/ -type f -name report.tsv -printf "%p ")
}

plasmids(){

  plasmid_db="/mnt/disk1/bin/plasmidid_db/plasmid.complete.nr100.fna"
  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir PLASMIDS_${run_name}

  for i in *R1.fastq.gz
  do
    hit=$(minimap2 -t $(nproc) -x sr $plasmid_db $i -t $(nproc) 2> /dev/null | \
            awk -F "\t" '$12 > 0 { print }' | cut -f 6 | sort | uniq -c | \
            sort -nr | head -1 | sed 's/^ *//')

    plasmid_acc=$(echo $hit | cut -d ' ' -f2)
    sample_name=$( echo $i | cut -d\_ -f1,2)

    if [ -z "$plasmid_acc" ]; then
      echo -e "$sample_name\tNF"
    else
      description=$(cat $plasmid_db | grep -m 1 -F $plasmid_acc | \
                    cut -d ' ' -f1 --complement | tr -d '>')
      reads=$(echo $hit | cut -d ' ' -f1)
      echo -e "$sample_name\t$reads\t$plasmid_acc\t$description"
    fi
  done > PLASMIDS_${run_name}/plasmid_candidates_${run_name}.tsv

  memory=`awk '{ printf "%.2f", $2/1024 ; exit}' /proc/meminfo | cut -d\. -f1`
  file=PLASMIDS_${run_name}/plasmid_id_NF_${run_name}.tsv

  ln -f /mnt/disk1/bin/plasmidid_db/plasmid.complete.nr100.fna .

  for r1 in *R1.fastq.gz
  do
    r2="${r1/R1/R2}"
    name="${r1%%_R1*}"
    cp ASSEMBLY/$name-idba-assembly.fa ${name}.fna
    grep -q "$name" PLASMIDS_${run_name}/plasmid_candidates_${run_name}.tsv

    if [ $? -eq 0 ]; then
      docker run --rm -it -v $(pwd):/data -w /data \
            -u $(id -u):$(id -g) plasmidid plasmidID.sh \
            -1 ${r1} -2 ${r2} -T $(nproc) \
            -d /data/plasmid.complete.nr100.fna -M $memory \
            -c ${name}.fna --no-trim -s ${name} -g PLASMIDS_${run_name}
      rm ${name}.fna

      out="PLASMIDS_${run_name}/${name}/images/${name}_summary.png"
      if [ ! -f "$out" ]; then
        rm -rf PLASMIDS_${run_name}/${name}
        echo -e "${name}\tNF" >> $file
      fi
    else
      :
    fi
  done

  if [ ! -s "$file" ]; then
    rm $file
  fi
  rm plasmid.complete.nr100.fna

}



#lis_type(){}

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

echo "SALMONELLA"
if [ -s "salm_id.txt" ]; then
  mkdir -p SALMONELLA

  while read -r fastq
  do
    find . -name "$fastq*fastq.gz" -type f -not -path "SALMONELLA/*" \
        -print0 | xargs -0 mv -t "SALMONELLA/" 2>/dev/null
  done < <(cat salm_id.txt | cut -f1)

  cd SALMONELLA

  file="../salm_id.txt"
  echo "run_salmonella"
  run_salmonella #&>> $log_file || error ${LINENO} $(basename $0)
  echo "Trimming"
  trimming
  echo "Assembly"
  assembly_idba
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
    find . -name "$fastq*fastq.gz" -type f -not -path "SALM-LIKE/*" \
        -print0 | xargs -0 mv -t "SALM-LIKE/" 2>/dev/null
  done < <(cat salm-like.txt | cut -f1)

  cd SALM-LIKE

  file="../salm-like.txt"
  echo "run_salmonella"
  run_salmonella
  echo "trimm"
  trimming
  echo "assembly"
  assembly_idba
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
  find . -name "$fastq*fastq.gz" -type f -not -path "OTHERS/*" \
      -print0 | xargs -0 mv -t "OTHERS/" 2>/dev/null
done < <(cat $file | cut -f1)

cd OTHERS

echo "Trimming"
trimming
echo "Assembly"
assembly_idba
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

run_name=$(basename `pwd` | cut -d\_ -f1)
path_results=`find . -type d -name "RESULTS"`
mkdir RESULTS_${run_name}

for i in $path_results
do
  cp $i/* RESULTS_${run_name}
done

end &>> $log_file
