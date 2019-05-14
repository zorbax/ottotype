#!/bin/bash

run_salmonella() {
  #echo -e "\n#Executing" ${FUNCNAME[0]} "\n"

  run_name=$(basename `pwd` | cut -d\_ -f1)

  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v &> /dev/null
  mkdir -p seqsero_$run_name $PWD/OUTPUT

  if [ -z $file ]; then
    file="salm_id.txt"
  fi
  for i in $(cat $file | cut -f1)
  do
    if [[ "find . -type l -name "$i*"" ]]; then
      R1=`ls -lt $i* | awk '{ print $NF }' | awk '($1 ~ /R1/) { print $1 }'`
      R2=`ls -lt $i* | awk '{ print $NF }' | awk '($1 ~ /R2/) { print $1 }'`
      docker run --rm -it -v $(pwd):/data -w /data seqsero SeqSero.py -m 2 -i $R1 $R2 &> /dev/null && \
      mv SeqSero_result* seqsero_$run_name
    else
      docker run --rm -it -v $(pwd):/data -w /data seqsero SeqSero.py -m 2 \
             -i $i\_R1.fastq.gz $i\_R2.fastq.gz &> /dev/null && mv SeqSero_result* seqsero_$run_name
    fi
  done

  find seqsero_$run_name -type f -name '*_result.txt' -exec cat {} \
                 > seqsero_$run_name/SeqSero\_$run_name\_serotype.txt \;
  cp seqsero_$run_name/SeqSero\_$run_name\_serotype.txt $PWD/OUTPUT

  echo "SeqSero: DONE"
  mkdir -p $PWD/RESULTS
  cat $PWD/OUTPUT/SeqSero\_$run_name\_serotype.txt | sed 's/potential monophasic variant of //; /^$/d' | \
       cut -d\_ -f1 | sed 's/\t/#/; /^$/d' | grep -i "input\|predicted" | grep -v formula | cut -d\# -f2 | \
       paste - - - | awk -F'\t' '{ print $1, $NF }' | sed 's/*//; s/N\/A/N\/A#/' | cut -d\# -f1 | \
       cut -d\( -f1 | sort | sed 's/ /\tSalmonella enterica subsp. enterica serovar /; s/*//' \
       > $PWD/RESULTS/seqsero_$run_name\_serotype.tsv

  cat  $PWD/OUTPUT/SeqSero\_$run_name\_serotype.txt | grep "Input\|antigen\|serotype" | grep -v -w are | \
       sed 's/ prediction//; s/O antigen/Antígeno/; s/Predicted antigenic profile/Perfil antigénico/;
       s/Predicted serotype/Serotipo/; s/(s)//; s/Input files:\t//; s/)//; s/\t/ /; /^$/d; s/*//;
       s/H[12] antigen(//; s/(O5-//; s/N\/A/N\/A\&/; s/.fastq/@/; s/ potential monophasic variant of //' | \
       cut -d\& -f1 | cut -d\@ -f1 | sed -E 's/_S[0-9]{1,}//' | grep -v "The predicted serotypes" | \
       paste - - - - - - | sort | sed 's/\t/\n/g' | cut -d\_ -f1 > $PWD/RESULTS/seqsero_$run_name\_serotype_sup.txt

  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v &> /dev/null

  mkdir -p srst2\_$run_name

  docker run --rm -it -v $(pwd):/data -w /data srst2 getmlst.py --species "Salmonella" &> /dev/null

  if [ -f "ARGannot.fasta" ]; then
    "ARGannot is already downloaded"
  else
    wget -q -nc https://raw.githubusercontent.com/CNRDOGM/bin/master/srst2/data/ARGannot.fasta
  fi

  docker run --rm -it -v $(pwd):/data -w /data srst2 srst2 --log --output /data/SRST2 \
         --input_pe *fastq.gz --forward R1 --reverse R2 --mlst_db Salmonella_enterica.fasta \
         --mlst_definitions senterica.txt --mlst_delimiter '_' --gene_db ARGannot.fasta \
         --threads $(nproc) &> /dev/null

  find . -maxdepth 1 -name "SRST2_*" -type f -exec mv {} srst2\_$run_name \;
  cp srst2\_$run_name/SRST2__compiledResults.txt $PWD/OUTPUT/srst2\_$run_name\_compiledResults.txt

  rm -f senterica.txt *tfa Salmonella* mlst* ARG* *bt2 *fai SRST2.log

  cat $PWD/OUTPUT/srst2\_$run_name\_compiledResults.txt | tail -n+2 | cut -d$'\t' -f 1-9 | \
      sed 's/[?*]//g; /^$/d' | perl -pe 's/_S[0-9]{1,}_//g' | sort \
      > $PWD/RESULTS/srst2_$run_name\_achtman.tsv

  cat $PWD/OUTPUT/srst2\_$run_name\_compiledResults.txt | tail -n+2 | \
      awk -v OFS='\t' -v f=2 -v t=13 '{
        for( i=1;i<=NF;i++ )
        if( i>=f&&i<=t )
          continue;
        else
          printf( "%s%s",$i,( i!=NF )?OFS:ORS )
      }' | sed 's/[?*]//g' | perl -pe 's/\t\-/#/g; s/\#{1,}//g;
      s/\_[0-9]{1,}//g' | perl -pe 's/_S[0-9]{1,}_//g; s/f{3,}/\tFAILED/' | \
      sed "s/''/-/" | sort > $PWD/RESULTS/srst2_$run_name\_argannot.tsv

  translate.py $PWD/RESULTS/srst2_$run_name\_argannot.tsv $PWD/RESULTS/antibiotics_$run_name.tsv
  translate.py $PWD/RESULTS/srst2_$run_name\_argannot.tsv $PWD/RESULTS/antibiotics_$run_name\_freq.tsv -freq

  while read line
  do
    id=`echo $line | cut -d ' ' -f1`
    st=`echo $line | awk '{ print $2 }'`
    cat $HOME/bin/Strain_Senterica.tsv | grep -w ^$st | \
      awk -v var="$id" -v OFS='\t' '{ print var, $0}' | awk 'BEGIN {FS="\t"} $10!="" {print}' \
      >> $PWD/RESULTS/srst2\_$run_name\_enterobase.tsv

    cat $HOME/bin/Strain_Senterica.tsv | grep -w ^$st | awk -v var="$id" -v OFS='\t' '{ print var, $0}' | \
                    awk 'BEGIN {FS="\t"} $10=="" {print}' >> null.tsv
  done < $PWD/RESULTS/srst2_$run_name\_achtman.tsv

  cat $PWD/RESULTS/srst2\_$run_name\_enterobase.tsv | awk '{ print $1, $2, $11, $10 }' | \
      sed 's/ /\nST:#/; s/ /\nST#Complex:#/; s/ /\nSerotipo:#/; s/#/ /g' \
      > $PWD/RESULTS/srst2\_$run_name\_enterobase_ST.tsv

  id1=`cat $PWD/RESULTS/srst2\_$run_name\_enterobase.tsv | awk '{ print $1}' | sort | uniq`
  id2=`cat null.tsv | awk '{ print $1}' | sort | uniq`
  diff <(echo "$id1") <(echo "$id2") | grep "\>" > $PWD/RESULTS/srst2\_$run_name\_enterobase_NF.tsv 2>/dev/null
  rm $PWD/null.tsv

  echo "SRST2: DONE"

  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v
  mkdir ariba_$run_name

  docker run --rm -it -v $(pwd):/data ariba ariba pubmlstget "Salmonella enterica" Salmonella &> /dev/null && \
  docker run --rm -it -v $(pwd):/data ariba ariba getref card card &> /dev/null && \
  docker run --rm -it -v $(pwd):/data ariba ariba prepareref -f card.fa -m card.tsv card.prepareref &> /dev/null

  for i in $(cat $file | cut -f1)
  do
    docker run --rm -it -v $(pwd):/data -w /data ariba ariba run /data/Salmonella/ref_db \
                              $i\_R1.fastq.gz $i\_R2.fastq.gz $i\_ariba &> /dev/null && \
    mv $i\_ariba ariba\_$run_name

    docker run --rm -it -v $(pwd):/data -w /data ariba ariba run /data/card.prepareref \
                             $i\_R1.fastq.gz $i\_R2.fastq.gz $i\_card &> /dev/null && \
    mv $i\_card ariba\_$run_name
  done

  rm -rf Salmonella card.*

  cd ariba_$run_name

  for i in *ariba
  do
    cd $i
    id=`echo $i | cut -d\_ -f1`
    if [ -e "mlst_report.tsv" ]; then
      cat mlst_report.tsv | tail -n+2 | awk -v var="$id" -v OFS='\t' '{ print var, $0}' | sort
      cd ..
    else
      cd ..
    fi
  done > ../OUTPUT/ariba_$run_name\_achtman.tsv

  cp ../OUTPUT/ariba_$run_name\_achtman.tsv ../RESULTS/ariba_$run_name\_achtman.tsv

  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v

  cd ..
  echo "ariba: DONE"
}
