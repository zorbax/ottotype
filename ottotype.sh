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

trimming() {

  mkdir -p TRIMMING/1U2U
  for r1 in *R1.fastq.gz
  do
    r2=${r1/R1/R2}
    name=${r1%%_R1*}
    trimmomatic PE -phred33 -threads $(nproc) $r1 $r2 \
              TRIMMING/${name}_R1.trim.fastq.gz TRIMMING/1U2U/${name}.1U.trim.fastq.gz \
              TRIMMING/${name}_R2.trim.fastq.gz TRIMMING/1U2U/${name}.2U.trim.fastq.gz \
              SLIDINGWINDOW:4:20 MINLEN:125 &> ${name}.trim.log

    if [ $? -eq 0 ]; then
      rm *trim.log
    fi
  done
}

assembly() {

  for r1 in TRIMMING/*R1.trim.fastq.gz
  do
    mkdir -p ASSEMBLY

    r2="${r1/R1/R2}"
    name="${r1##*/}"; name="${name%%_R1*}"

    echo "${name}"
    sga preprocess -q 25 -f 20 --pe-mode=1 ${r1} ${r2} \
        > ${name}_12.t.pp.fq 2> /dev/null
    median=`cat ${name}_12.t.pp.fq | awk 'NR%4==2 { print length }' | head -20000 | \
            Rscript -e 'summary (as.numeric (readLines ("stdin")))' | tail -n+2 | \
            awk '{ print $3 }'  | cut -d\. -f1`

    if [ $median -le 200 ]; then
      sga index -t $(nproc) -a ropebwt ${name}_12.t.pp.fq \
          > ${name}.index.out 2> ${name}.index.err
      sga correct -t $(nproc) ${name}_12.t.pp.fq -o ${name}_12.t.pp.ec.fq \
          > ${name}.correct.out 2> ${name}.correct.err
    else
      sga index -t $(nproc) -a sais ${name}_12.t.pp.fq \
          > ${name}.index.out 2> ${name}.index.err
      sga correct -t $(nproc) ${name}_12.t.pp.fq -o ${name}_12.t.pp.ec.fq \
          > ${name}.correct.out 2> ${name}.correct.err
    fi

    sed -n '1~4s/^@/>/p;2~4p' ${name}_12.t.pp.ec.fq > ${name}_12.t.pp.ec.fa
    medianec=`cat ${name}_12.t.pp.ec.fa | awk '$0 !~ /^>/ { print length }' | \
              head -20000 | Rscript -e 'summary (as.numeric (readLines ("stdin")))' | \
              tail -n+2 | awk '{ print $3}'  | cut -d\. -f1`

    if [[ $medianec -ge 128 ]]; then
      idba_ud500 -r ${name}_12.t.pp.ec.fa -o ${name} --mink 35 --maxk 249 \
          --num_threads $(nproc) --min_pairs 2 > /dev/null
    else
      idba_ud -r ${name}_12.t.pp.ec.fa -o ${name} --mink 35 --maxk 124 \
          --num_threads $(nproc) --min_pairs 2 > /dev/null
    fi

    if [ $? != 0 ]; then
       cp ${name}/contig.fa ${name}-idba-assembly.fa
    else
       cp ${name}/scaffold.fa ${name}-idba-assembly.fa
    fi

    rm -rf *pp.{fq,bwt,sai,rbwt} *out *err *rsai *ec.{fq,fa} ${name}

    cat ${name}-idba-assembly.fa | sed ':a;N;/^>/M!s/\n//;ta;P;D' | \
         awk '/^>/ { getline seq } length(seq) >500 { print $0 "\n" seq }' \
         > ASSEMBLY/${name}-idba-assembly.fa
    rm ${name}-idba-assembly.fa
  done
}

assembly_spades() {

mkdir -p ASSEMBLY
memory=$(awk '{ printf "%.2f", $2/1024/1024 ; exit}' /proc/meminfo | cut -d\. -f1)

for r1 in TRIMMING/*R1.trim.fastq.gz
do
  r2="${r1/R1/R2}"
  name="${r1##*/}"; name="${name%%_R1*}"
  spades.py --pe1-1 $r1 --pe1-2 $r2 \
            --pe1-s TRIMMING/1U2U/${name}.1U.trim.fastq.gz \
            --pe1-s TRIMMING/1U2U/${name}.2U.trim.fastq.gz \
            -o ${name}_spades -t $(nproc) -m $memory &>/dev/null

  if [ $? != 0 ]; then
    cp ${name}/contig.fa ${name}-spades-assembly.fa
  else
    cp ${name}/scaffolds.fasta ${name}-spades-assembly.fa
  fi

find ${name}_spades -maxdepth 2 -type f -name 'scaffolds.fasta' -exec cp {} ${name}.tmp \;

cat ${name}.tmp | sed ':a;N;/^>/M!s/\n//;ta;P;D' | \
     awk '/^>/ { getline seq } length(seq) >500 { print $0 "\n" seq }' \
     > ASSEMBLY/${name}-spades-assembly.fa

rm -rf ${name}_spades ${name}.tmp
done
}

assembly_stats_cov() {

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p OUTPUT ASSEMBLY/Stats RESULTS
  echo -e "Assembly\tNumber_of_contigs\tCoverage\tAssembly_size\tLargest_contig\tN50\tN90" \
      > assembly_$run_name.stats

  for i in ASSEMBLY/*assembly.fa
  do
    assembly=ASSEMBLY/$(basename "$i" | cut -d\- -f3 --complement)
    name=`echo $assembly | cut -d\/ -f2`
    reads=TRIMMING/$(basename $i | cut -d\- -f1)

    echo -e "Assembly:\t${name}" | tee ${name}.stats
    bwa index -a bwtsw $i -p ${name} &> /dev/null
    bwa mem -t $(nproc) ${name} $reads\_R1.trim.fastq.gz $reads\_R2.trim.fastq.gz \
         2> /dev/null | samtools view -Sb -F4 - | \
         samtools sort -@ $(nproc) - -o ${name}.mapped.sorted.bam 2>/dev/null
    samtools index ${name}.mapped.sorted.bam
    rm ${name}.{bwt,pac,ann,amb,sa}

    cov=$(samtools mpileup ${name}.mapped.sorted.bam 2> /dev/null | \
      awk '{ count++ ; sum += $4 } END { printf "%s\t%.2f\n", "Coverage:", sum/count }')
    cover=`echo $cov | awk '{ print $2 }'`

    cat $i | awk '!/^>/ { printf "%s", $0 n = "\n" } /^>/
      { print n $0 n = "" } END { printf "%s", n }'| \
      sed '/^>/ d'| awk '{ print length($0) }' | sort -gr > ${name}_contig_lengths.stat

    contigs_number=$(cat ${name}_contig_lengths.stat | wc -l)
    assembly_size=$(paste -sd+ ${name}_contig_lengths.stat | bc)
    awk 'BEGIN {sum=0} {sum= sum+$0; print sum}' ${name}_contig_lengths.stat \
        > ${name}_contig_lengths_cum.stat

    awk -v var=$assembly_size 'BEGIN {FS=OFS=","} {for (i=1; i<=NF; i++) $i/=var;}1' \
                                ${name}_contig_lengths_cum.stat > ${name}_cum_perc.stat

    paste ${name}_contig_lengths.stat ${name}_cum_perc.stat > ${name}_matrix.stat

    N50=$(awk '$2 >= 0.50' ${name}_matrix.stat | head -1 | awk '{ print $1 }')
    N90=$(awk '$2 >= 0.90' ${name}_matrix.stat | head -1 | awk '{ print $1 }')
    large_contig=$(head -1 ${name}_contig_lengths.stat)

    rm *stat
    echo -e "Number of contigs:\t$contigs_number\n$cov" >> ${name}.stats
    echo -e "Assembly size:\t$assembly_size" >> ${name}.stats
    echo -e "Largest contig:\t$large_contig" >> ${name}.stats
    echo -e "N50:\t$N50\nN90:\t$N90" >> ${name}.stats
    mv ${name}.stats ASSEMBLY/Stats
    echo -e "${name}\t$contigs_number\t$cover\t$assembly_size\t$large_contig\t$N50\t$N90" \
        >> assembly_$run_name.stats
    find . -maxdepth 1 -type f \( -name "*.bam" -o -name "*.bai" \) -exec rm {} \;
  done
  mv assembly_$run_name.stats OUTPUT/ && cp OUTPUT/assembly_$run_name.stats RESULTS/
}

assembly_mlst(){
  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p ASSEMBLY_MLST_${run_name} OUTPUT RESULTS
  cd ASSEMBLY

  for i in *assembly.fa
  do
    mlst --csv $i --threads $(nproc) \
        >> ../ASSEMBLY_MLST_${run_name}/assembly_mlst_${run_name}.csv 2> /dev/null
  done

  cd ..

  cat ASSEMBLY_MLST_${run_name}/assembly_mlst_${run_name}.csv | \
      sed 's/-assembly.fa//; s/,/\t/g' | sort | \
      tee OUTPUT/assembly_mlst_${run_name}.tsv \
      RESULTS/assembly_mlst_${run_name}.tsv > /dev/null

  cat OUTPUT/assembly_mlst_${run_name}.tsv | sed 's/^/\n>/; s/\t/\tST /2; s/\t/\n/g' \
      > RESULTS/assembly_mlst_${run_name}.list
}

kmer_finder(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p KMERFINDER_${run_name}/{raw,genome}_${run_name} RESULTS/

  KmerFinder_DB=/mnt/disk1/bin/KmerFinder_DB
  kraw="KMERFINDER_$run_name/raw_$run_name"
  kgenome="KMERFINDER_$run_name/genome_$run_name"

  for r1 in *R1.fastq.gz
  do
    name="${r1%%_R1*}"
    docker run --rm -it -v $KmerFinder_DB:/database -v $(pwd):/workdir \
          -u $(id -u):$(id -g) -w /data kmerfinder -i /workdir/${r1} \
          -o /workdir/$kraw/${name} -db /database/bacteria.ATG \
          -tax /database/bacteria.name -x &> /dev/null
  done

  for i in ASSEMBLY/*assembly.fa
  do
    genome_name="${i##*/}"; name="${name%%-*}"
    echo $genome_name
    #docker run --rm -it -v $KmerFinder_DB:/database -v $(pwd):/workdir \
    #      -u $(id -u):$(id -g) -w /data kmerfinder -i /workdir/$i \
    #      -o /workdir/$kgenome/$genome_name -db /database/bacteria.ATG \
    #      -tax /database/bacteria.name -x &> /dev/null
  done

  for i in *R1.fastq.gz
  do
    name="${i%%_R1*}"
    echo "# RAW_${name}"
    cat $kraw/${name}/results.txt | tail -n+2 | awk -F'\t' -v OFS='\t' '{ print $NF, $3, $13}'
    echo
    echo "# GENOME_${name}"
    cat $kgenome/${name}/results.txt | tail -n+2 | awk -F'\t' -v OFS='\t' '{ print $NF, $3, $13}'
    echo
  done > RESULTS/kmerfinder_${run_name}.txt
}

kraken_tax(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p KRAKEN2_${run_name} OUTPUT RESULTS

  YGGDRASIL="/mnt/disk2/bin/Kraken2/yggdrasil"

  for r1 in *R1.fastq.gz
  do
    r2="${r1/R1/R2}"
    name="${r1%%_R1*}"

    echo "$i"  # --use-mpa-style
    docker run --rm -it -v $YGGDRASIL:/database -v $(pwd):/data \
         -u $(id -u):$(id -g) -w /data kraken2 kraken2 --paired \
         --gzip-compressed --threads $(nproc) --db /database \
         --report KRAKEN2_${run_name}/${name}.kraken2-report.tsv \
         ${r1} ${r2} > /dev/null 2> KRAKEN2_${run_name}/${name}.kraken2.log

    cat KRAKEN2_${run_name}/${name}.kraken2.log | tail -2 | paste - - | \
        sed "s/^  /${name}\t/" | sponge KRAKEN2_${run_name}/${name}.kraken2.log
    # rs and sponge > sudo apt-get install moreutils rs
    tax=$(cat KRAKEN2_${run_name}/${name}.kraken2-report.tsv | \
          awk -F'\t' '{if($1>5) print }' | grep -P '\t[S]\t'| \
          awk -F'\t' '{ print $6, "#"$1}' | sed -E 's/^ {1,}//; s/[ ]{1,}/_/g' | \
          rs -TeC | sed 's/_#_/ /g')

    echo -e "${name}\t$tax" >> OUTPUT/kraken2_${run_name}.tax.tsv
  done

  cp OUTPUT/kraken2_${run_name}.tax.tsv RESULTS/kraken2_${run_name}.tax.tsv
}

antibiotics(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p OUTPUT RESULTS ANTIBIOTICS_$run_name

  repo="https://raw.githubusercontent.com/CNRDOGM"
  wget -q -nc $repo/srst2/master/data/ARGannot_r3.fasta
  docker_cmd="docker run --rm -it -v $(pwd):/data -w /data"

  $docker_cmd srst2 srst2 --log --output /data/SRST2 --input_pe *fastq.gz \
     --forward R1 --reverse R2 --gene_db ARGannot_r3.fasta --threads $(nproc) &> /dev/null

  find . -maxdepth 1 -name "SRST2_*" -type f -not -path "ANTIBIOTICS_$run_name/*" \
      -exec mv {} ANTIBIOTICS_$run_name/ \;

  cat ANTIBIOTICS_$run_name/SRST2__genes__ARGannot_r3__results.txt | \
      tail -n +2 | sed 's/[?*]//g' | sed -E 's/_S[0-9]{1,}_//' | \
      perl -pe "s/\t\-/#/g; s/\#{1,}//g; s/\_[0-9]{1,}//g; s/\'\'/\-/g;
      s/f{3,}//" | sort > RESULTS/antibiotics_${run_name}_argannot.tsv

  translate.py RESULTS/antibiotics_${run_name}_argannot.tsv RESULTS/antibiotics_r3_${run_name}.tsv
  translate.py RESULTS/antibiotics_${run_name}_argannot.tsv RESULTS/antibiotics_r3_${run_name}_freq.tsv -freq

  $docker_cmd ariba ariba getref card card &> /dev/null && \
  $docker_cmd ariba ariba prepareref -f card.fa -m card.tsv card.prepareref &> /dev/null

  for r1 in *R1.fastq.gz
  do
    r2="${r1/R1/R2}"
    name="${r1%%_R1*}"
    $docker_cmd ariba ariba run /data/card.prepareref \
                ${r1} ${r2} ${name}_card &> /dev/null && \
    mv ${name}_card ANTIBIOTICS_$run_name
  done

  rm -rf card.* *ARGannot*{bt2,fasta,fai} SRST2.log
  $docker_cmd ariba ariba summary --preset all out.summary \
          $(find ARIBA_PIPELINE/ -type f -name report.tsv -printf "%p ")
}

plasmids(){

  plasmid_db="/mnt/disk1/bin/plasmidid_db/plasmid.complete.nr100.fna"
  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir PLASMIDS_${run_name}

  for i in *R1.fastq.gz
  do
    hit=$(minimap2 -t $(nproc) -x sr $plasmid_db $i -t $(nproc) 2> /dev/null | \
            awk -F "\t" '$12 > 0 { print }' | cut -f 6 | sort | uniq -c | \
            sort -nr | head -1 | sed 's/^ *//')

    plasmid_acc=$(echo $hit | cut -d ' ' -f2)
    sample_name=$( echo $i | cut -d\_ -f1,2)

    if [ -z "$plasmid_acc" ]; then
      echo -e "$sample_name\tNF"
    else
      description=$(cat $plasmid_db | grep -m 1 -F $plasmid_acc | \
                    cut -d ' ' -f1 --complement | tr -d '>')
      reads=$(echo $hit | cut -d ' ' -f1)
      echo -e "$sample_name\t$reads\t$plasmid_acc\t$description"
    fi
  done > PLASMIDS_${run_name}/plasmid_candidates_${run_name}.tsv

  memory=`awk '{ printf "%.2f", $2/1024 ; exit}' /proc/meminfo | cut -d\. -f1`
  file=PLASMIDS_${run_name}/plasmid_id_NF_${run_name}.tsv

  ln -f /mnt/disk1/bin/plasmidid_db/plasmid.complete.nr100.fna .

  for r1 in *R1.fastq.gz
  do
    r2="${r1/R1/R2}"
    name="${r1%%_R1*}"
    cp ASSEMBLY/$name-idba-assembly.fa ${name}.fna
    grep -q "$name" PLASMIDS_${run_name}/plasmid_candidates_${run_name}.tsv

    if [ $? -eq 0 ]; then
      docker run --rm -it -v $(pwd):/data -w /data \
            -u $(id -u):$(id -g) plasmidid plasmidID.sh \
            -1 ${r1} -2 ${r2} -T $(nproc) \
            -d /data/plasmid.complete.nr100.fna -M $memory \
            -c ${name}.fna --no-trim -s ${name} -g PLASMIDS_${run_name}
      rm ${name}.fna

      out="PLASMIDS_${run_name}/${name}/images/${name}_summary.png"
      if [ ! -f "$out" ]; then
        rm -rf PLASMIDS_${run_name}/${name}
        echo -e "${name}\tNF" >> $file
      fi
    else
      :
    fi
  done

  if [ ! -s "$file" ]; then
    rm $file
  fi
  rm plasmid.complete.nr100.fna

}

ecoli_type(){

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p {SRST2Ec_$run_name,OUTPUT,RESULTS}

  repo="https://raw.githubusercontent.com/CNRDOGM"
  docker_cmd="docker run --rm -it -v $(pwd):/data -u $(id -u):$(id -g) -w /data"

  wget -q -nc $repo/srst2/master/data/EcOH.fasta

  $docker_cmd srst2 srst2 --log --output /data/SRST2_EcOH \
      --input_pe *fastq.gz --forward R1 --reverse R2 \
      --gene_db EcOH.fasta --threads $(nproc) &> /dev/null

   wget -q -nc $repo/srst2/master/data/{LEE_mlst.fasta,LEE_profiles.txt}
   $docker_cmd srst2 srst2 --log --output /data/SRST2_LEE \
       --input_pe *fastq.gz --forward R1 --reverse R2 \
       --mlst_db LEE_mlst.fasta --mlst_definitions LEE_profiles.txt \
       --threads $(nproc) &> /dev/null

  wget -q -nc "$repo/srst2/master/data/ARGannot_r3.fasta" -O ARGannot.fasta
  $docker_cmd srst2 srst2 --log --output /data/SRST2_ARG \
       --input_pe *fastq.gz --forward R1 --reverse R2 \
       --gene_db ARGannot.fasta --threads $(nproc) &> /dev/null

  $docker_cmd srst2 getmlst.py --species "Escherichia coli#1" &> /dev/null && \
  $docker_cmd srst2 srst2 --log --output /data/SRST2_Ec1 \
       --input_pe *fastq.gz --forward R1 --reverse R2 \
       --mlst_db Escherichia_coli#1.fasta --mlst_definitions ecoli.txt \
       --mlst_delimiter '_' --threads $(nproc) &> /dev/null

  $docker_cmd srst2 getmlst.py --species "Escherichia coli#2" &> /dev/null && \
  $docker_cmd srst2 srst2 --log --output /data/SRST2_Ec2 \
       --input_pe *fastq.gz --forward R1 --reverse R2 \
       --mlst_db Escherichia_coli#2.fasta --mlst_definitions ecoli_2.txt \
       --mlst_delimiter '_' --threads $(nproc) &> /dev/null

  find . -maxdepth 1 -iname "*results.txt" -type f -exec mv {} SRST2Ec_${run_name} \;
  rm -f EcOH* LEE* ARG* *tfa SRST2_*.{bam,pileup} Escherichia* ARG* ecoli* *log
  #mv SRST*log LOGS

  cat SRST2Ec_${run_name}/SRST2_EcOH__fullgenes__EcOH__results.txt | tail -n+2 | \
     sort -nk 1 | sed -E 's/_S[0-9]{1,}_//; s/_//' | awk -v OFS='\t' '{ print $1, $3, $NF }' | \
     sed -E 's/\t[A-Z0-9]{1,}.[0-9]{1};/\t/; s/\t/\tGen: /; s/;/: /' | \
     awk -F'\t' -v OFS='\t' '{x=$1; $1=""; a[x]=a[x]$0 } END { for (x in a) print x,a[x] }' | \
     sort -nk 1 | sed -E 's/\t{1,}Gen/\nGen/g; s/\t/; /g; s/polyermase/polymerase/' \
     > RESULTS/srst2Ec_${run_name}_EcOH.tsv

  cat SRST2Ec_${run_name}/SRST2_ARG__genes__ARGannot__results.txt | sed 's/[?*]//g' | \
     perl -pe 's/\t\-/#/g; s/\#{1,}//g;s/\_[0-9]{1,}//g'| perl -pe 's/_S[0-9]{1,}_//g;
     s/_//; s/f{3,}/\tFAILED/' | sed "s/''/-/" | tail -n+2 > RESULTS/srst2Ec_${run_name}_argannot.tsv

  translate.py RESULTS/srst2Ec_${run_name}_argannot.tsv RESULTS/antibioticsEc_${run_name}.tsv
  translate.py RESULTS/srst2Ec_${run_name}_argannot.tsv RESULTS/antibioticsEc_${run_name}_freq.tsv -freq

  cat SRST2Ec_${run_name}/SRST2_Ec1__mlst__Escherichia_coli#1__results.txt | \
      cut -d$'\t' -f 1-9 | sed 's/[?*]//g; /^$/d' | perl -pe 's/_S[0-9]{1,}_//g;
      s/_//' | sort | sed 's/ST/ST1/' > SRST2Ec_${run_name}/srst2Ec_${run_name}_achtman#1.tsv

  cat SRST2Ec_${run_name}/SRST2_Ec2__mlst__Escherichia_coli#2__results.txt | \
      cut -d$'\t' -f 1-9 | sed 's/[?*]//g; /^$/d' | perl -pe 's/_S[0-9]{1,}_//g;
      s/_//' | sort | sed 's/ST/ST2/' > SRST2Ec_${run_name}/srst2Ec_${run_name}_achtman#2.tsv

  paste SRST2Ec_${run_name}/srst2Ec_${run_name}_achtman#1.tsv \
      SRST2Ec_${run_name}/srst2Ec_${run_name}_achtman#2.tsv | \
      awk -v OFS='\t' '$1 == $10 { $10="" ; print}' \
      > RESULTS/srst12Ec_${run_name}_achtman12.tsv
  mv SRST2Ec_${run_name} OUTPUT
}

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
  assembly
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
  assembly
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
assembly
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
