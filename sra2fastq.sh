#!/bin/bash
 
#PNUSAS_ID/SAMN_ID > SRA ID > FASTQ
 
list=`cat $1`
 
echo -e "SampleName\tRun\tBioProject\tSample\tBioSample\tScientificName" > pnusas2sra.tsv
 
for i in ${list}
do
  table=`esearch -db sra -query ${i} | efetch -format runinfo | \
  awk -F',' -v OFS='\t' '{ print $30, $1, $22, $25, $26, $29}' | tail -n+2`

  echo $table >> pnusas2sra.tsv
  sra=`echo $table | awk '{ print $2 }'`
  fastq-dump --split-files --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --gzip $sra
done
