#!/bin/bash

display_usage(){
  echo -e "\nUsage:"
  echo -e "\t$(basename $0) SampleSheet.csv\n"
}

if [ $# -eq 1 ]; then
  display_usage
  exit 1
fi

bcl2fastq -R . --no-lane-splitting \
    -r $(nproc) -p $(nproc) -w $(nproc) \
    --fastq-compression-level 9 \
    --minimum-trimmed-read-length 70 \
    -l NONE -o fastq --sample-sheet $1
