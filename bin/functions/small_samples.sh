#!/bin/bash

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
    echo "# Not found small samples in dataset." &>> $log_file
  else
    cat small_samples/small_samples_size.txt | cut -d\: -f1 | sort | uniq \
                                       > small_samples/small_samples_ids.txt
    n_samples=`cat small_samples/small_samples_ids.txt | wc -l`
    echo -e "\n# This $n_samples sample(s) are too small:\n" &>> $log_file
    echo -e "ID\tReads" &>> $log_file

    for i in $(cat small_samples/small_samples_size.txt | cut -d\: -f1 | sort | uniq)
    do
      echo -e "${i}\t$(zcat ${i}*R1.fastq.gz | awk 'NR%4==1' | wc -l )"
      mv -i ${i}_R*.fastq.gz small_samples/
    done &>> $log_file
  fi
}
