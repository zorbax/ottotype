#!/bin/bash
 
#PNUSAS_ID/SAMN_ID > SRA ID > FASTQ
 
list=`cat $1`
 
echo -e "SampleName\tRun\tBioProject\tSample\tBioSample\tScientificName" > pnusas2sra.tsv
 
for i in ${list}
do
  table=`esearch -db sra -query ${i} | efetch -format runinfo | \
  awk -F',' -v OFS='\t' '{ print $30, $1, $22, $25, $26, $29}' | \
  tail -n+2`

  echo -e "$table\n" | sort | uniq | sed '/^$/d' >> pnusas2sra.tsv
  sra=`echo "$table" | grep -v Run | awk '{ print $2 }' | sed '/^$/d' | sort | uniq`
  echo $sra >> tmp.txt
  for k in ${sra}
  do
    echo -e "${i}\t${k}"
    fastq-dump --split-files --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --gzip ${k}
  done
done
