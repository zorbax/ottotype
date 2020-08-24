#!/bin/bash

clean_filenames() {
    # Rename MiSeq (001) and NextSeq (L001, 001) files
    find . -maxdepth 1 -name "*fastq.gz" -type f -or -type l | \
        rename 's/_L001//; s/_001//; s/_1.fastq/_S01_R1.fastq/ ;
        s/_2.fastq/_S01_R2.fastq/; s/-//'

}

# shellcheck disable=SC2181
if [[ $? -eq 0 ]]; then
    echo -e "\n\n# The filenames were renamed with the ${FUNCNAME[0]} function"
else
    echo "Cleaning filenames failed"
fi