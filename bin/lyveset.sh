#!/bin/bash

#set -e

echo ""
echo "=============================="
echo " LyveSET - SSB-CNRDOGM v0.9.4 "
echo "=============================="
echo ""

display_usage() {
    echo -e "\nUsage:"
    echo -e "\t$(basename $0) -r reference.fna"
    echo -e "\t$(basename $0) -r reference.fna -d [Docker enabled]"
    echo -e "\t$(basename $0) -n [Only trim and combo-clean]\n"
}

tmpmk(){
    if [[ -O $PWD/tmp && -d $PWD/tmp ]]; then
        TMPDIR=$PWD/tmp
    else
        rm -rf $PWD/tmp 2> /dev/null
        mkdir -p $PWD/tmp
        TMPDIR=$(mktemp -d $PWD/tmp/XXXX)
    fi

    TMP=$TMPDIR
    TEMP=$TMPDIR
    export TMPDIR TMP TEMP
}

tmprm(){
    if [[ -O $TMPDIR && -d $TMPDIR ]]; then
        rm -rf "${TMPDIR:?}/"*
    fi
}

use_docker=false
no_ref=false

while getopts ":r:dnh" opt
do
    case $opt in
        r ) reference="$OPTARG"
            ;;
        d ) use_docker=true
            ;;
        n ) no_ref=true
            ;;
        h ) display_usage
            ;;
        ? ) echo "Invalid option -$OPTARG" >&2
            exit 1
            ;;
    esac
done

if [ $OPTIND -eq 1 ]; then
  display_usage
  exit 1
fi

ls *fastq.gz | grep -v SRR | cut -d\_ -f1,2 | sort | uniq > original_names.txt
#Rename SRR files
find -maxdepth 1 -type f -or -type l -regex  '.*\(_1.fastq.gz\|_2.fastq.gz\)$' | \
                        rename 's/_1.fastq/_R1.fastq/ ; s/_2.fastq/_R2.fastq/'
find -maxdepth 1 -name "*fastq.gz" -type f -or -type l | \
                        rename 's/_L001// ; s/_S[0-9]{1,}// ; s/_001//'

dir="set_$(basename `pwd`)"
docker_run="docker run --rm -it -v $(pwd):/data -u $(id -u):$(id -g) -w /data"

docker_clean() {
    docker ps --filter status=dead --filter status=exited -aq | \
                            xargs -r docker rm -v 2>/dev/null
    docker images --no-trunc | grep '<none>' | awk '{ print $3 }' | \
                            xargs -r docker rmi 2>/dev/null
}

if [ $use_docker == true ]; then
  tmpmk
  docker_clean
  if [ ! -d "$dir" ]; then
    $docker_run lyveset set_manage.pl --create set_$(basename `pwd`) 2>/dev/null
  fi

  cleaned=(`find . -maxdepth 1 -name "*cleaned.fastq.gz"`)
  interleave=(`find . -maxdepth 1 -name "*combo.fastq.gz"`)

  if [ ${#cleaned[@]} -gt 0 ]; then
    for i in $(ls *cleaned.fastq.gz | cut -d\_ -f1 | sort | uniq )
    do
      echo "Add reads $i"
      $docker_run lyveset set_manage.pl set_$(basename `pwd`) \
                  --add-reads ${i}_cleaned.fastq.gz &> /dev/null
      echo "Done"
    done
  elif [ ${#interleave[@]} -gt 0 ]; then
    for i in $(ls *combo.fastq.gz | cut -d\_ -f1 | sort | uniq )
    do
      echo "TrimClean $i"
      $docker_run lyveset run_assembly_trimClean.pl -i ${i}_combo.fastq.gz \
           -p 2 -quieter --nosingletons --min_avg_quality 24 --min_quality 30 \
           --bases_to_trim 20 --numcpus 2 --min_length 62 \
           -o ${i}_cleaned.fastq.gz &> /dev/null
      echo "Add reads"
      $docker_run lyveset set_manage.pl set_$(basename `pwd`) \
                     --add-reads ${i}_cleaned.fastq.gz &> /dev/null
      rm *combo.fastq.gz
    echo "Done"
    done
  else
    for i in $(ls *fastq.gz | cut -d\_ -f1 | sort | uniq )
    do
      echo "Trimming $i"
      $docker_run lyveset trimmomatic PE -phred33 -threads "$(nproc)" \
            ${i}_R1.fastq.gz ${i}_R2.fastq.gz \
            ${i}_R1.trim.fastq.gz ${i}.1U.trim.gz \
            ${i}_R2.trim.fastq.gz ${i}.2U.trim.gz \
            SLIDINGWINDOW:5:20 MINLEN:70 &> /dev/null
      echo "Combo"
      paste <(zcat ${i}_R1.trim.fastq.gz ) <(zcat ${i}_R2.trim.fastq.gz) | \
         paste - - - - | awk -v OFS="\n" -v FS="\t" '{ print $1, $3, $5, $7, \
         $2, $4, $6, $8 }' | pigz --best --processes "$(nproc)" > ${i}_combo.fastq.gz
      rm *trim*
      echo "TrimClean"
      $docker_run lyveset run_assembly_trimClean.pl -i ${i}_combo.fastq.gz \
          -p 2 -quieter --nosingletons --min_avg_quality 24 --min_quality 30 \
          --bases_to_trim 20 --numcpus 2 --min_length 62 \
          -o ${i}_cleaned.fastq.gz &>/dev/null
      rm *combo*
      echo "Add reads"
      $docker_run lyveset set_manage.pl set_$(basename `pwd`) \
                          --add-reads ${i}_cleaned.fastq.gz &> /dev/null
      echo "Done"
    done
  fi

  if [ $no_ref == false ]; then
    echo "Change reference"
    $docker_run lyveset set_manage.pl set_$(basename `pwd`) --change-reference "$reference"
    echo "Set"
    $docker_run lyveset launch_set.pl set_$(basename `pwd`) -ref "$reference" \
              --min_coverage 20 --min_alt_frac 0.95 --allowedFlanking 5 \
              --mask-phages 1 --numcpus "$(nproc)" &>/dev/null
    rm *cleaned.fastq.gz
  fi

  echo "DONE"
  tmprm
else
  tmpmk
  if [ ! -d "$dir" ]; then
    set_manage.pl --create set_$(basename `pwd`) 2>/dev/null
  fi

  cleaned=(`find . -maxdepth 1 -name "*cleaned.fastq.gz"`)
  interleave=(`find . -maxdepth 1 -name "*combo.fastq.gz"`)

  if [ ${#cleaned[@]} -gt 0 ]; then
    for i in $(ls *cleaned.fastq.gz | cut -d\_ -f1 | sort | uniq )
    do
      echo "Add reads $i"
      set_manage.pl set_$(basename `pwd`) \
                    --add-reads ${i}_cleaned.fastq.gz &> /dev/null
      echo "Done"
    done
  elif [ ${#interleave[@]} -gt 0 ]; then
    for i in $(ls *combo.fastq.gz | cut -d\_ -f1 | sort | uniq )
    do
      echo "TrimClean $i"
      run_assembly_trimClean.pl -i ${i}_combo.fastq.gz -p 2 -quieter \
         --nosingletons --min_avg_quality 24 --min_quality 30 \
         --bases_to_trim 20 --numcpus 2 --min_length 62 \
         -o ${i}_cleaned.fastq.gz &> /dev/null
      rm *combo.fastq.gz
      echo "Add reads"
      set_manage.pl set_$(basename `pwd`) \
                    --add-reads ${i}_cleaned.fastq.gz &> /dev/null
      echo "Done"
    done
  else
    for i in $(ls *fastq.gz | cut -d\_ -f1 | sort | uniq )
    do
      echo "Trimming $i"
      trimmomatic PE -phred33 -threads "$(nproc)" \
          ${i}_R1.fastq.gz ${i}_R2.fastq.gz \
          ${i}_R1.trim.fastq.gz ${i}.1U.trim.gz \
          ${i}_R2.trim.fastq.gz ${i}.2U.trim.gz \
          SLIDINGWINDOW:5:20 MINLEN:70 &> /dev/null
      echo "Combo"
      paste <(zcat ${i}_R1.trim.fastq.gz ) <(zcat ${i}_R2.trim.fastq.gz) | \
         paste - - - - | awk -v OFS="\n" -v FS="\t" '{ print $1, $3, $5, $7, \
         $2, $4, $6, $8 }' | pigz --best --processes "$(nproc)" > ${i}_combo.fastq.gz
      rm *trim*
      echo "TrimClean"
      run_assembly_trimClean.pl -i ${i}_combo.fastq.gz -p 2 -quieter --nosingletons \
           --min_avg_quality 24 --min_quality 30 --bases_to_trim 20 --numcpus 2 \
           --min_length 62 -o ${i}_cleaned.fastq.gz &>/dev/null
      rm *combo*
      echo "Add reads"
      #check if links exist don't add, if not remove all links and link again
      set_manage.pl set_$(basename `pwd`) --add-reads ${i}_cleaned.fastq.gz &> /dev/null
      echo "Done"
    done
  fi

  if [ $no_ref == false ]; then
    echo "Change reference"
    #check if links exist don't add, if not remove all links and link again
    set_manage.pl set_$(basename `pwd`) --change-reference "$reference"
    echo "Set"
    launch_set.pl set_$(basename `pwd`) -ref "$reference" --min_coverage 20 \
                --min_alt_frac 0.95 --allowedFlanking 5 --mask-phages 1 \
                --numcpus "$(nproc)" &>/dev/null
  fi
  echo "DONE"
  tmprm
fi
