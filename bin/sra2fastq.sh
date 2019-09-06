#!/bin/bash

#PNUSAS_ID/SAMN_ID > SRA ID > FASTQ

if [ -f "$1" ]; then
  list=$(cat $1)
else
  list=$@
fi

echo -e "SampleName\tRun\tBioExperiment\tBioProject\t\
         Sample\tBioSample\tLibraryStrategy\tLibrarySelection\t\
         LibrarySource\tScientificName" > pnusas2sra.tsv

for i in ${list}
do
  table=$(esearch -db sra -query ${i} | efetch -format runinfo | \
          awk -F',' -v OFS='\t' '{ print $30, $1, $11, $22, $25,
              $26, $13, $14, $15, $29}' | tail -n+2)

  echo -e "$table\n" | sort | uniq | sed '/^$/d' >> pnusas2sra.tsv
  sra=$(echo "$table" | grep -v Run | awk '{ print $2 }' | \
        sed '/^$/d' | sort | uniq)
  echo $sra >> sra.txt
  for k in ${sra}
  do
    echo -e "${i}\t${k}"
    fastq-dump --split-files --defline-seq '@$ac.$si.$sg/$ri' \
               --defline-qual '+' --gzip ${k}
    # fasterq-dump ${k} -t /dev/shm -e $(nproc)
    rename 's/_1.fastq/_S01_R1_L001.fastq/ ; s/_2.fastq/_S01_R2_L001.fastq/' *gz
  done
done

rm -rf $HOME/ncbi
