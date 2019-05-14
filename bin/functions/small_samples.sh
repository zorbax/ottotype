#!/bin/bash

small_samples() {

  mkdir $PWD/small_samples

  for i in *fastq.gz
  do
    if [[ -L "$i" && -f "$i" ]]; then
      if [[ "$(( $( stat -Lc '%s' $i ) / 1024 / 1024 ))" -le "4" ]]; then
        small_links=`ls -l $i | awk '{ print $NF }'`
        ls -lh $small_links | awk '{ print $NF": " $5}' | \
        sed 's/_L001//; s/_001// ; s/.*\///; s/_R[12].fastq.gz//' \
        >> $PWD/small_samples/small_samples_size.txt
      fi
    else
      find . -maxdepth 1 -size -4M -name "*fastq.gz" -exec ls -lh {} \; | \
           awk '{ print $NF": " $5}' | cut -d\/ -f2 | \
           sed 's/_R[12].fastq.gz//' >> $PWD/small_samples/small_samples_size.txt
    fi
  done

  x=`ls -l $PWD/small_samples/small_samples_size.txt | awk '{ print $5 }'`
  if [ $x == 0 ]; then
    rm -rf $PWD/small_samples
    echo "# Not found small samples in dataset."
  else
    small="$PWD/small_samples/small_samples_size.txt"
    id_samples=`cat $small | awk '{ print $1 }' | cut -d\_ -f1,2 | cut -d\: -f1 | sort | uniq`
    echo $id_samples | cut -d\: -f1 > $PWD/small_samples/small_samples_ids.txt
    n_samples=`cat $id_samples | wc -l`
    echo -e "\n# This $n_samples are too small:\n# $(echo $id_samples)"

    for i in $id_samples
    do
      echo "${i}\t$(zcat ${i}*R1.fastq.gz | awk 'NR%4==1' | wc -l )"
    done

    for i in $id_samples
    do
      mv $i\_R*.fastq.gz $PWD/small_samples
    done
  fi
}
