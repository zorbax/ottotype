#!/bin/bash

find -maxdepth 1 -name "*fastq.gz" -type f -or -type l | \
     rename 's/_L001//; s/_001//; s/_1.fastq/_S01_R1.fastq/;
     s/_2.fastq/_S01_R2.fastq/'

trimming() {

  mkdir -p $PWD/TRIMMING/1U2U
  for r1 in *R1.fastq.gz
  do
    r2=${r1/R1/R2}
    name=${r1%%_R1*}
    trimmomatic PE -phred33 -threads $(nproc) $r1 $r2 \
              TRIMMING/${name}_R1.trim.fastq.gz TRIMMING/1U2U/${name}.1U.trim.fastq.gz \
              TRIMMING/${name}_R2.trim.fastq.gz TRIMMING/1U2U/${name}.2U.trim.fastq.gz \
              SLIDINGWINDOW:4:20 MINLEN:70 &> ${name}.trim.log
    rm ${name}.trim.log
  done
}

assembly_spades(){
  mkdir -p ASSEMBLY
  memory=$(awk '{ printf "%.2f", $2/1024/1024 ; exit}' /proc/meminfo | cut -d\. -f1)

  for r1 in *TRIMMING/*R1.trim.fastq.gz
  do
    r2="${r1/R1/R2}"
    name="${r1##*/}"; name="${name%%_R1*}"
    spades.py --pe1-1 $r1 --pe1-2 $r2 \
              --pe1-s TRIMMING/1U2U/${name}.1U.trim.fastq.gz \
              --pe1-s TRIMMING/1U2U/${name}.2U.trim.fastq.gz \
              -o ${name}_spades -t $(nproc) -m $memory &>/dev/null

  find ${name}_spades -maxdepth 2 -type f -name 'scaffolds.fasta' -exec cp {} ${name}.tmp \;

  cat ${name}.tmp | sed ':a;N;/^>/M!s/\n//;ta;P;D' | \
       awk '/^>/ { getline seq } length(seq) >500 { print $0 "\n" seq }' \
       > ASSEMBLY/${name}-spades-assembly.fa

  rm -rf ${name}_spades ${name}.tmp
done
}

if [ ! -d "ASSEMBLY" ]; then
  trimming && assembly_spades
else
  samples=$(ls *fastq.gz | cut -d\_ -f1,2 | sort | uniq)
  assemblies=$(ls ASSEMBLY/*fa | cut -d\- -f1 | cut -d\/ -f2 | sort  | uniq)
  difference=$(diff <(echo "$samples") <(echo "$assemblies"))

  if [[ ! -z $difference ]]; then

    mkdir -p TRIMMING/1U2U
    for i in $(echo $difference | cut -d ' ' -f3)
    do
      trimmomatic PE -phred33 -threads $(nproc) ${i}_R1.fastq.gz ${i}_R2.fastq.gz \
                TRIMMING/${i}_R1.trim.fastq.gz TRIMMING/1U2U/${i}.1U.trim.fastq.gz \
                TRIMMING/${i}_R2.trim.fastq.gz TRIMMING/1U2U/${i}.2U.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:70 &> ${i}.trim.log

      rm ${i}.trim.log

      memory=$(awk '{ printf "%.2f", $2/1024/1024 ; exit}' /proc/meminfo | cut -d\. -f1)

      spades.py --pe1-1 TRIMMING/${i}_R1.trim.fastq.gz \
                --pe1-2 TRIMMING/${i}_R2.trim.fastq.gz \
                --pe1-s TRIMMING/1U2U/${i}.1U.trim.fastq.gz \
                --pe1-s TRIMMING/1U2U/${i}.2U.trim.fastq.gz \
                -o ${i}_spades -t $(nproc) -m $memory &>/dev/null

      find ${i}_spades -maxdepth 1 -type f -name 'scaffolds.fasta' -exec cp {} ${i}.tmp \;

      cat ${i}.tmp | sed ':a;N;/^>/M!s/\n//;ta;P;D' | \
          awk '/^>/ { getline seq } length(seq) >500 { print $0 "\n" seq }' \
          > ASSEMBLY/${i}-spades-assembly.fa

      rm -rf ${i}_spades ${i}.tmp
    done
  fi
fi

plasmid_db="/mnt/disk1/bin/plasmidid_db/plasmid.complete.nr100.fna"
run_name=$(basename `pwd` | cut -d\_ -f1)
mkdir -p PLASMIDS_${run_name}

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

ln -f /mnt/disk1/bin/plasmidid_db/plasmid.complete.nr100.fna .

for r1 in *R1.fastq.gz
do
  r2="${r1/R1/R2}"
  name="${r1%%_R1*}"
  memory=$(awk '{ printf "%.2f", $2/1024 ; exit}' /proc/meminfo | cut -d\. -f1)
  cp ASSEMBLY/$name-spades-assembly.fa ${name}.fna
  grep -q "$name" PLASMIDS_${run_name}/plasmid_candidates_${run_name}.tsv

  if [ $? -eq 0 ]; then
    docker run --rm -it -v $(pwd):/data -w /data \
          -u $(id -u):$(id -g) plasmidid plasmidID.sh \
          -1 ${r1} -2 ${r2} -T $(nproc) \
          -d /data/plasmid.complete.nr100.fna -s ${name} \
          -M $memory -c ${name}.fna -g PLASMIDS_${run_name} \
          --no-trim &> /dev/null

    rm ${name}.fna

    out="PLASMIDS_${run_name}/${name}/images/${name}_summary.png"
    file=PLASMIDS_${run_name}/plasmid_id_NF_${run_name}.tsv

    if [ ! -f "$out" ]; then
      rm -rf PLASMIDS_${run_name}/${name}
      echo -e "${name}\tNF" >> $file
    fi
  else
    :
  fi
done

if [ ! -s "$file" ]; then
  rm -f $file 2> /dev/null
fi

rm plasmid.complete.nr100.fna
