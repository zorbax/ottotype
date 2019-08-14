#!/bin/bash

plasmids(){

  plasmid_db="/mnt/disk1/bin/plasmidid_db/plasmid.complete.nr100.fna"
  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p PLASMIDS_${run_name} RESULTS

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

  memory=$(awk '{ printf "%.2f", $2/1024 ; exit}' /proc/meminfo | cut -d\. -f1)
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
      result=$(cat PLASMIDS_${run_name}/${name}/data/${name}.karyotype_summary.txt | \
               sed 's/chr - //' | cut -d\. -f1 | tr '\n' ' ' | sed 's/ $//')
      echo -e "$name\t$result" >> RESULTS/plasmids_${run_name}.tsv

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
