#!/bin/bash

set -e

echo "===================================="
echo "Serotyping pipeline from SSB-CNRDOGM"
echo "               v0.9f                "
echo "===================================="

####   _____                 _   _
####  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___
####  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
####  |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
####  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
####


clean() {
  find -maxdepth 1 -name "*fastq.gz" -type f -or -type l | rename 's/_L001//; s/_001//';
}

small_samples() {

  mkdir $PWD/SMALL_SAMPLES

  for i in *fastq.gz
  do
    if [[ -L "$i" && -f "$i" ]]; then
      if [[ "$(( $( stat -Lc '%s' $i ) / 1024 / 1024 ))" -le "4" ]]; then
        small_links=`ls -l $i | awk '{ print $NF }'`
        ls -lh $small_links | awk '{ print $NF": " $5}' | \
        sed 's/_L001//; s/_001// ; s/.*\///; s/_R[12].fastq.gz//' \
        >> $PWD/SMALL_SAMPLES/small_samples_size.txt
      fi
    else
      find . -maxdepth 1 -size -4M -name "*fastq.gz" -exec ls -lh {} \; | \
           awk '{ print $NF": " $5}' | cut -d\/ -f2 | \
           sed 's/_R[12].fastq.gz//' >> $PWD/SMALL_SAMPLES/small_samples_size.txt
    fi
  done

  x=`ls -l $PWD/SMALL_SAMPLES/small_samples_size.txt | awk '{ print $5 }'`
  if [ $x == 0 ]; then
    rm -rf $PWD/SMALL_SAMPLES
    echo "Not found small samples in dataset. Great!" #This to log
  else
    small="$PWD/SMALL_SAMPLES/small_samples_size.txt"
    id_samples=`cat $small | awk '{ print $1 }' | cut -d\_ -f1,2 | cut -d\: -f1 | sort | uniq`
    echo $id_samples | cut -d\: -f1 > $PWD/SMALL_SAMPLES/small_samples_ids.txt
    n_samples=`cat $id_samples | wc -l`
    echo -e "
        The next $n_samples samples are too small and
        will be excluded from the further analysis:
        $(echo $id_samples)"

    for i in $id_samples
    do
      mv $i\_R*.fastq.gz $PWD/SMALL_SAMPLES
    done
  fi
}

screen_tax() {
  salmonella.py -d . -e _R1.fastq.gz | grep -v file | sed 's/_R[12].fastq.gz//' | \
      awk -F'\t' -v OFS='\t' '{if($6 == "100") { print $1, $6 > "salm_id.txt"} else { print $1, $6 > "nosalm_id.txt"}}'

  for i in *R1.fastq.gz
  do
    acc=$(minimap2 -t $(nproc) -x sr $HOME/bin/16S/NCBI.gz $i 2> /dev/null | awk -F "\t" '$12>0{print}' | \
          cut -f 6 | sort | uniq -c | sort -nr | head -n 1 | sed 's/^ *//' | cut -d ' ' -f 2)

    if [ -z "$acc" ]; then
      echo "No match found!"
    else
      echo "Reads: $( echo $i | cut -d\_ -f1,2)"
      echo "Top hit: $acc"
      hit=$(gzip -c -d -f $HOME/bin/16S/NCBI.gz | grep -m 1 -F $acc)
      echo "Description: $hit"
      binomial=$(echo $hit | cut -d ' ' -f 2,3)
      echo "Species: $binomial"
    fi
  done > screen_tax_raw_ncbi.txt

  cat screen_tax_raw_ncbi.txt | grep "Reads\|Species" | sed 's/ /#/' | cut -d\# -f2 | \
      grep -B1 Salmonella | sed 's/--//; /^$/d' | paste - - | sort -k1 > salm_id_ncbi.txt
  cat screen_tax_raw_ncbi.txt | grep "Reads\|Species" | sed 's/ /#/' | cut -d\# -f2 | \
      sed 's/--//; /^$/d' | paste - - | grep -v Salmonella | sed '/^$/d' | sort -k1 > nosalm_id_ncbi.txt

  if [ -s "salm_id.txt" ]; then
    #awk -F'\t' 'NR==FNR{salmid[$1]=$1;ncbi[$1]=$1;next};
    #      (1!=salmid[$1]){ print $1 }' <(sort -k1 salm_id.txt | cut -f1) <(sort -k2 salm_id_ncbi.txt | cut -f1) \
    #      > salm-like.txt 2> /dev/null
    grep -vwif <(sort -k1 salm_id.txt | cut -f1) <(cut -f1 salm_id_ncbi.txt) > salm-like.txt
  fi

  for i in *txt
  do
    if [ ! -s "$i" ]; then
      rm -f $i
    fi
  done
}

run_salmonella() {
  run_name=$(basename `pwd` | cut -d\_ -f1)

  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v
  mkdir -p SEQSERO_$run_name $PWD/OUTPUT

  if [ -z $file ]; then
    file="salm_id.txt"
  fi
  for i in $(cat $file | cut -f1)
  do
    if [[ "find . -type l -name "$i*"" ]]; then
      R1=`ls -lt $i* | awk '{ print $NF }' | awk '($1 ~ /R1/) { print $1 }'`
      R2=`ls -lt $i* | awk '{ print $NF }' | awk '($1 ~ /R2/) { print $1 }'`
      docker run --rm -it -v $(pwd):/data -w /data seqsero SeqSero.py -m 2 -i $R1 $R2 &> /dev/null && \
      mv SeqSero_result* SEQSERO_$run_name
    else
      docker run --rm -it -v $(pwd):/data -w /data seqsero SeqSero.py -m 2 \
             -i $i\_R1.fastq.gz $i\_R2.fastq.gz &> /dev/null && mv SeqSero_result* SEQSERO_$run_name
    fi
  done

  find SEQSERO_$run_name -type f -name '*_result.txt' -exec cat {} \
                 > SEQSERO_$run_name/SeqSero\_$run_name\_serotype.txt \;
  cp SEQSERO_$run_name/SeqSero\_$run_name\_serotype.txt $PWD/OUTPUT

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

  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v

  mkdir -p SRST2_$run_name

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

  find . -maxdepth 1 -name "SRST2_*" -type f -exec mv {} SRST2_$run_name \;
  cp SRST2_$run_name/SRST2__compiledResults.txt $PWD/OUTPUT/SRST2_$run_name\_compiledResults.txt

  rm -f senterica.txt *tfa Salmonella* mlst* ARG* *bt2 *fai SRST2.log

  cat $PWD/OUTPUT/SRST2_$run_name\_compiledResults.txt | tail -n+2 | cut -d$'\t' -f 1-9 | \
      sed 's/[?*]//g; /^$/d' | perl -pe 's/_S[0-9]{1,}_//g' | sort \
      > $PWD/RESULTS/srst2_$run_name\_achtman.tsv

  cat $PWD/OUTPUT/SRST2_$run_name\_compiledResults.txt | tail -n+2 | \
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
  mkdir ARIBA_$run_name

  docker run --rm -it -v $(pwd):/data ariba ariba pubmlstget "Salmonella enterica" Salmonella &> /dev/null && \
  docker run --rm -it -v $(pwd):/data ariba ariba getref card card &> /dev/null && \
  docker run --rm -it -v $(pwd):/data ariba ariba prepareref -f card.fa -m card.tsv card.prepareref &> /dev/null

  for i in $(cat $file | cut -f1)
  do
    docker run --rm -it -v $(pwd):/data -w /data ariba ariba run /data/Salmonella/ref_db \
                              $i\_R1.fastq.gz $i\_R2.fastq.gz $i\_ariba &> /dev/null && \
    mv $i\_ariba ARIBA\_$run_name

    docker run --rm -it -v $(pwd):/data -w /data ariba ariba run /data/card.prepareref \
                             $i\_R1.fastq.gz $i\_R2.fastq.gz $i\_card &> /dev/null && \
    mv $i\_card ARIBA\_$run_name
  done

  rm -rf Salmonella card.*

  cd ARIBA_$run_name

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

  docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v

  cd ..
  echo "ARIBA: DONE"
}

trimming() {

  mkdir 1U2U
  for i in $(ls *fastq.gz | cut -d\_ -f1,2 | sort | uniq )
  do
    trimmomatic PE -phred33 -threads $(nproc) $i\_R1.fastq.gz $i\_R2.fastq.gz \
       $i\_R1.trim.fastq.gz $i.1U.trim.fastq.gz \
       $i\_R2.trim.fastq.gz $i.2U.trim.fastq.gz \
       SLIDINGWINDOW:5:20 MINLEN:70 &> $i.trim.log
    mv *U.trim.fastq.gz 1U2U

    if [ $? -eq 0 ]; then
        rm $i.trim.log
    fi
  done
}

assembly() {

  for i in $(ls *trim.fastq.gz | grep -v U | cut -d\_ -f1,2 | sort | uniq)
  do
    mkdir -p $PWD/ASSEMBLY
    echo "$i"
    sga preprocess -q 25 -f 20 --pe-mode=1  $i\_R1.trim.fastq.gz $i\_R2.trim.fastq.gz > $i\_12.t.pp.fq 2> /dev/null
    median=`cat $i\_12.t.pp.fq | awk 'NR%4==2 { print length }' | head -10000 | Rscript -e 'summary (as.numeric (readLines ("stdin")))' | tail -n+2 | awk '{ print $3}'  | cut -d\. -f1`

    if [ $median -le 200 ]; then
      #mem=`cat /proc/meminfo | grep MemTotal | awk '{ print $2}'`
      #kb=$(echo "($mem*0.5)" | bc | cut -d\. -f1)
      sga index -t $(nproc) -a ropebwt $i\_12.t.pp.fq > $i.index.out 2> $i.index.err
      sga correct -t $(nproc) $i\_12.t.pp.fq -o $i\_12.t.pp.ec.fq > $i.correct.out 2> $i.correct.err
    else
      sga index -t $(nproc) -a sais $i\_12.t.pp.fq > $i.index.out 2> $i.index.err
      sga correct -t $(nproc) $i\_12.t.pp.fq -o $i\_12.t.pp.ec.fq > $i.correct.out 2> $i.correct.err
    fi

    sed -n '1~4s/^@/>/p;2~4p' $i\_12.t.pp.ec.fq > $i\_12.t.pp.ec.fa
    medianec=`cat $i\_12.t.pp.ec.fa | awk '$0 !~ /^>/ { print length }' | head -10000 | Rscript -e 'summary (as.numeric (readLines ("stdin")))' | tail -n+2 | awk '{ print $3}'  | cut -d\. -f1`

    if [[ $medianec -ge 128 ]]; then
      idba_ud500 -r $i\_12.t.pp.ec.fa -o $i --mink 35 --maxk 249 --num_threads $(nproc) --min_pairs 2 > /dev/null
    else
      idba_ud -r $i\_12.t.pp.ec.fa -o $i --mink 35 --maxk 124 --num_threads $(nproc) --min_pairs 2 > /dev/null
    fi

    if [ $? != 0 ]; then
       cp $i/contig.fa $i-assembly.fa
    else
       cp $i/scaffold.fa $i-assembly.fa
    fi

    rm -rf *pp.{fq,bwt,sai,rbwt} *out *err *rsai *ec.{fq,fa} $i
    mv $i-assembly.fa $PWD/ASSEMBLY
  done
}

assembly_stats_cov() {

  run_name=$(basename `pwd` | cut -d\_ -f1)
  echo -e "Assembly\tNumber_of_contigs\tCoverage\tAssembly_size\tLargest_contig\tN50\tN90" > assembly_$run_name.stats

  for i in $PWD/ASSEMBLY/*assembly.fa
  do
    name=`basename $i | cut -d\_ -f1 | cut -d\- -f1`
    reads=`basename $i | cut -d\- -f1`

    echo -e "Assembly:\t$name" | tee $name.stats
    bwa index -a bwtsw $i -p $name 2> /dev/null
    bwa mem -t $(nproc) $name $reads\_R1.trim.fastq.gz $reads\_R2.trim.fastq.gz 2> /dev/null |\
              samtools view -Sb -F4 - | samtools sort - $name.mapped.sorted
    samtools index $name.mapped.sorted.bam
    rm $name.{bwt,pac,ann,amb,sa}

    cov=$(samtools mpileup $name.mapped.sorted.bam 2> /dev/null | \
      awk '{ count++ ; sum += $4 } END { printf "%s\t%.2f\n", "Coverage:", sum/count }')
    cover=`echo $cov | awk '{ print $2 }'`

    cat $i | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }'| \
                                 sed '/^>/ d'| awk '{ print length($0) }' | sort -gr > $name\_contig_lengths.stat

    contigs_number=$(cat $name\_contig_lengths.stat | wc -l)
    assembly_size=$(paste -sd+ $name\_contig_lengths.stat | bc)
    awk 'BEGIN {sum=0} {sum= sum+$0; print sum}' $name\_contig_lengths.stat > $name\_contig_lengths_cum.stat

    awk -v var=$assembly_size 'BEGIN {FS=OFS=","} {for (i=1; i<=NF; i++) $i/=var;}1' \
                                    $name\_contig_lengths_cum.stat > $name\_cum_perc.stat

    paste $name\_contig_lengths.stat $name\_cum_perc.stat > $name\_matrix.stat

    N50=$(awk '$2 >= 0.50' $name\_matrix.stat | head -1 | awk '{ print $1 }')
    N90=$(awk '$2 >= 0.90' $name\_matrix.stat | head -1 | awk '{ print $1 }')
    large_contig=$(head -1 $name\_contig_lengths.stat)

    rm *stat
    echo -e "Number of contigs:\t$contigs_number" | tee -a $name.stats
    echo "$cov" | tee -a $name.stats
    echo -e "Assembly size:\t$assembly_size" | tee -a $name.stats
    echo -e "Largest contig:\t$large_contig" | tee -a $name.stats
    echo -e "N50:\t$N50" | tee -a $name.stats
    echo -e "N90:\t$N90" | tee -a $name.stats
    echo -e "$name\t$contigs_number\t$cover\t$assembly_size\t$large_contig\t$N50\t$N90" | tee -a assembly_$run_name.stats
    find . -maxdepth 1 -type f \( -name "*.bam" -o -name "*.bai" \) -exec rm {} \;
  done
}

assembly_mlst(){
  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p $run_name\_mlst OUTPUT
  cd $PWD/ASSEMBLY
  for i in *assembly.fa
  do
    mlst --csv $i --threads $(nproc) >> ../$run_name\_mlst/assembly_$run_name\_mlst.csv 2> /dev/null
  done
  cd ..
  cat $run_name\_mlst/assembly_$run_name\_mlst.csv | sed 's/-assembly.fa//; s/,/\t/g' | \
        sort  > $run_name\_mlst/assembly_$run_name\_mlst.tsv
  cp $run_name\_mlst/assembly_$run_name\_mlst.tsv OUTPUT/assembly_$run_name\_mlst.csv
}

kmer_finder(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir $run_name\_kmer

  DB="/mnt/disk1/bin/kmerfinder_DB/bacteria.organisms.ATGAC"

  for i in $PWD/ASSEMBLY/*assembly.fa
  do
    genome_name=`basename $i | cut -d\- -f1`
    findTemplate -i $i -t $DB -x ATGAC -w -o $run_name\_kmer/$genome_name.kmer
  done
}

kraken_tax(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p $run_name\_kraken2 OUTPUT RESULTS

  for i in $(ls *fastq.gz | grep -v trim | cut -d\_ -f1,2 | sort | uniq)
  do
    kraken2 --paired --gzip-compressed --threads $(nproc) \
      --db $YGGDRASIL --report $run_name\_kraken2/$i.kraken2-report.tsv \
      $i\_R1.fastq.gz $i\_R2.fastq.gz > /dev/null 2> $run_name\_kraken2/$i.kraken2.log

    # sponge > sudo apt-get install moreutils
    tax=`cat $run_name\_kraken2/$i.kraken2-report.tsv | awk -F'\t' '{if($1>5) print }' | \
         grep -P '\t[DPCOFGS]\t' | sed 's/D/K/' | awk -F'\t' '$4=tolower($4){ print $4"_", $6}' | \
         sed -E 's/[ ]{1,}/_/g' | tr  "\n" ";" | sed 's/;$/\n/'`
    echo -e "$i\t$tax" > $run_name\_kraken2/$i.tax.tsv
    cat $run_name\_kraken2/$i.kraken2.log | tail -2 | paste - - | sed "s/^  /$i\t/" | sponge $run_name\_kraken2/$i.kraken2.log
  done

  cat $run_name\_kraken2/*tax.tsv | awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1'| \
      sed -e 's/k__/#/; s/s__/#/; s/\(#\).*\(#\)//' | sed -E 's/_S[0-9]{1,}//' | \
      awk -F'\t' -v OFS='\t' '{gsub(";s_"," |",$2);gsub("_"," ",$2)}1' | \
      perl -pe 'if(/\#/){s/\ /\_/g}' > OUTPUT/$run_name\_kraken.tax.tsv
  cp OUTPUT/$run_name\_kraken.tax.tsv RESULTS/$run_name\_kraken.tax.tsv
}

####   __  __       _
####  |  \/  | __ _(_)_ __
####  | |\/| |/ _` | | '_ \
####  | |  | | (_| | | | | |
####  |_|  |_|\__,_|_|_| |_|
####

docker images --no-trunc | grep '<none>' | awk '{ print $3 }' | xargs -r docker rmi
memory=`awk '{ printf "%.2f", $2/1024/1024 ; exit}' /proc/meminfo | cut -d\. -f1`
run_name=$(basename `pwd` | cut -d\_ -f1)

echo "Clean"
clean
echo "Small samples"
small_samples
echo "Taxa screening"
screen_tax

echo "SALMONELLA"
if [ -s "salm_id.txt" ]; then
  mkdir -p SALMONELLA

  while read -r fastq
  do
    find . -type f -name "$fastq*fastq.gz" -exec mv -f {} SALMONELLA \; 2>/dev/null
  done < <(cat salm_id.txt | cut -f1)

  cd SALMONELLA

  file="../salm_id.txt"
  run_salmonella
  trimming
  assembly
  assembly_stats_cov
  cd ..
fi

echo "SALM-LIKE"

if [ -s "salm-like.txt" ]; then
  mkdir -p SALM-LIKE

  while read -r fastq
  do
    find . -name "$fastq*fastq.gz" -exec mv -f {} SALM-LIKE \; 2>/dev/null
  done < <(cat salm-like.txt | cut -f1)

  cd SALM-LIKE

  file="../salm-like.txt"
  echo "run_salmonella"
  run_salmonella
  echo "trimm"
  trimming
  echo "assembly"
  assembly
  echo "stats"
  assembly_stats_cov
  echo "Assembly mlst"
  assembly_mlst
  echo "kmerfinder"
  kmer_finder
  cd ..
fi

echo "NOSALM"

if [ -s "nosalm_id_ncbi.txt" ]; then

  mkdir -p OTHERS

  while read -r fastq
  do
    find . -name "$fastq*fastq.gz" -exec mv -f {} OTHERS \;
  done < <(cat nosalm_id_ncbi.txt | cut -f1)

  cd OTHERS

  file="../nosalm_id_ncbi.txt"
  echo "Trimming"
  trimming
  echo "Assembly"
  assembly
  echo "Assembly stats"
  assembly_stats_cov
  echo "Assembly mlst"
  assembly_mlst
  echo "Kmerfinder"
  kmer_finder
  echo "Kraken"
  kraken_tax
  cd ..
fi
