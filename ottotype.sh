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
source bin/functions/assembly_stats.sh
source bin/functions/assembly_mlst.sh
source bin/functions/antibiotics.sh
source bin/functions/kraken_tax.sh
source bin/functions/plasmids.sh

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
