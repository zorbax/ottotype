#!/bin/bash

kraken_tax(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p kraken2\_$run_name OUTPUT RESULTS

  for i in $(ls *fastq.gz | grep -v trim | cut -d\_ -f1,2 | sort | uniq)
  do
    kraken2 --paired --gzip-compressed --threads $(nproc) \
      --db $YGGDRASIL --report kraken2\_$run_name/$i.kraken2-report.tsv \
      $i\_R1.fastq.gz $i\_R2.fastq.gz > /dev/null 2> kraken2\_$run_name/$i.kraken2.log

    # sponge > sudo apt-get install moreutils
    tax=`cat kraken2\_$run_name/$i.kraken2-report.tsv | awk -F'\t' '{if($1>5) print }' | \
         grep -P '\t[DPCOFGS]\t' | sed 's/D/K/' | awk -F'\t' '$4=tolower($4){ print $4"_", $6}' | \
         sed -E 's/[ ]{1,}/_/g' | tr  "\n" ";" | sed 's/;$/\n/'`
    echo -e "$i\t$tax" > kraken2\_$run_name/$i.tax.tsv
    cat kraken2\_$run_name/$i.kraken2.log | tail -2 | paste - - | sed "s/^  /$i\t/" | sponge kraken2\_$run_name/$i.kraken2.log
  done

  cat kraken2\_$run_name/*tax.tsv | awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1'| \
      sed -e 's/k__/#/; s/s__/#/; s/\(#\).*\(#\)//' | sed -E 's/_S[0-9]{1,}//' | \
      awk -F'\t' -v OFS='\t' '{gsub(";s_"," |",$2);gsub("_"," ",$2)}1' | \
      perl -pe 'if(/\#/){s/\ /\_/g}' | sed 's/;p  /#/; s/|/#|/; s/\(#\).*\(#\)/ /' \
      > OUTPUT/kraken\_$run_name.tax.tsv
  cp OUTPUT/kraken\_$run_name.tax.tsv RESULTS/kraken\_$run_name.tax.tsv
}
