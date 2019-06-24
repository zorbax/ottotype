#!/bin/bash

docker_cmd="docker run --rm -it -v $(pwd):/data -w /data"

for i in virulencefinder vfdb_full plasmidfinder
do
  $docker_cmd ariba ariba getref ${i} ${i}
  $docker_cmd ariba ariba prepareref -f /data/${i}.fa -m /data/${i}.tsv ${i}.prepareref
done

for r1 in *R1.fastq.gz
do
  r2="${r1/R1/R2}"
  name="${r1%%_R1*}"
  echo "$r1 $r2 $name"
  for k in virulencefinder vfdb_full plasmidfinder
  do
    $docker_cmd ariba ariba run /data/${k}.prepareref \
                   ${r1} ${r2} ${name}_${i} &> /dev/null
  done
done
