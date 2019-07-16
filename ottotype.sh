#!/bin/bash

####   _____                 _   _
####  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___
####  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
####  |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
####  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
####

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
for i in $dir/bin/functions/*sh
do
  source $i
done

T=$(date +%s)

docker images --no-trunc | grep '<none>' | awk '{ print $3 }' | xargs -r docker rmi
run_name=$(basename `pwd` | cut -d\_ -f1)

logfile
tmpmk

echo "Clean"
clean
echo "Checklist"
checklist
echo "Small samples"
small_samples
echo "Screen tax"
screen_tax

if [ -s "SCREENING/salm_id.txt" ]; then
  mkdir -p SALMONELLA

  while read -r fastq
  do
    find . -name "$fastq*fastq.gz" -type f -not -path "SALMONELLA/*" \
        -print0 | xargs -0 mv -t "SALMONELLA/" 2>/dev/null
  done < <(cat SCREENING/salm_id.txt | cut -f1)

  cd SALMONELLA
  file="../SCREENING/salm_id.txt"
  echo "SALMONELLA"
  run_salmonella &>> $log_file || error ${LINENO} $(basename $0)
  echo "Trimming"
  trimming &>> $log_file || error ${LINENO} $(basename $0)
  #echo "Assembly"
  #assembly_idba &>> $log_file || error ${LINENO} $(basename $0)
  #echo "Stats"
  #assembly_stats_cov &>> $log_file || error ${LINENO} $(basename $0)
  #echo "Virulence"
  #plas_vir_vfdb &>> $log_file || error ${LINENO} $(basename $0)
  #echo "Plasmids"
  #plasmids &>> $log_file || error ${LINENO} $(basename $0)
  cd ..
fi

: <<'END'
echo "SALM-LIKE"

if [ -s "SCREENING/salm-like.txt" ]; then
  mkdir -p SALM-LIKE

  while read -r fastq
  do
    find . -name "$fastq*fastq.gz" -type f -not -path "SALM-LIKE/*" \
        -print0 | xargs -0 mv -t "SALM-LIKE/" 2>/dev/null
  done < <(cat SCREENING/salm-like.txt | cut -f1)

  cd SALM-LIKE
  file="../SCREENING/salm-like.txt"
  echo "SALM-LIKE"
  run_salmonella &>> $log_file || error ${LINENO} $(basename $0)
  echo "Trimming"
  trimming &>> $log_file || error ${LINENO} $(basename $0)
  echo "Assembly"
  assembly_idba &>> $log_file || error ${LINENO} $(basename $0)
  echo "Stats"
  assembly_stats_cov &>> $log_file || error ${LINENO} $(basename $0)
  echo "Virulence"
  plas_vir_vfdb &>> $log_file || error ${LINENO} $(basename $0)
  echo "Assembly MLST"
  assembly_mlst &>> $log_file || error ${LINENO} $(basename $0)
  echo "KmerFinder"
  kmer_finder &>> $log_file || error ${LINENO} $(basename $0)
  echo "Kraken2"
  kraken_tax &>> $log_file || error ${LINENO} $(basename $0)
  cd ..
fi

echo "OTHERS"

if [ -s "SCREENING/nosalm_id_ncbi.txt" ]; then
  file="SCREENING/nosalm_id_ncbi.txt"
else
  file="SCREENING/nosalm_id.txt"
fi

mkdir -p OTHERS

while read -r fastq
do
  find . -name "$fastq*fastq.gz" -type f -not -path "OTHERS/*" \
      -print0 | xargs -0 mv -t "OTHERS/" 2>/dev/null
done < <(cat $file | cut -f1)

cd OTHERS
trimming &>> $log_file || error ${LINENO} $(basename $0)
assembly_idba &>> $log_file || error ${LINENO} $(basename $0)
assembly_stats_cov &>> $log_file || error ${LINENO} $(basename $0)
plas_vir_vfdb &>> $log_file || error ${LINENO} $(basename $0)
assembly_mlst &>> $log_file || error ${LINENO} $(basename $0)
kmer_finder &>> $log_file || error ${LINENO} $(basename $0)
kraken_tax &>> $log_file || error ${LINENO} $(basename $0)
antibiotics &>> $log_file || error ${LINENO} $(basename $0)
cd ..

END

run_name=$(basename `pwd` | cut -d\_ -f1)
path_results=`find . -type d -name "RESULTS"`
mkdir -p RESULTS_${run_name}

for i in $path_results
do
  cp $i/* RESULTS_${run_name}
done

tmprm
end &>> $log_file
