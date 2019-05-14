#!/bin/bash

screen_tax() {
  salmonella.py -d . -e _R1.fastq.gz | grep -v file | sed 's/_R[12].fastq.gz//' | \
      awk -F'\t' -v OFS='\t' '{if($6 == "100") { print $1, $6 > "salm_id.txt"} else { print $1, $6 > "nosalm_id.txt"}}'

  for i in *R1.fastq.gz
  do
    acc=$(minimap2 -t $(nproc) -x sr $HOME/bin/16S/NCBI.gz $i 2> /dev/null | awk -F "\t" '$12>0{print}' | \
          cut -f 6 | sort | uniq -c | sort -nr | head -n 1 | sed 's/^ *//' | cut -d ' ' -f 2)

    if [ -z "$acc" ]; then
      echo "No match found!"
    else
      echo "Reads: $( echo $i | cut -d\_ -f1,2)"
      echo "Top hit: $acc"
      hit=$(gzip -c -d -f $HOME/bin/16S/NCBI.gz | grep -m 1 -F $acc)
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
