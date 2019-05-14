#!/bin/bash

clean() {
  find -maxdepth 1 -name "*fastq.gz" -type f -or -type l | \
     rename 's/_L001//; s/_001//; s/_1.fastq/_S01_R1.fastq/ ;
     s/_2.fastq/_S01_R2.fastq/'

  echo -e "\n\n# The filenames were renamed with the ${FUNCNAME[0]} function" &>> $log_file
}
