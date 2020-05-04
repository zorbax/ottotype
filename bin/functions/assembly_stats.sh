#!/bin/bash

bwa mem -t $(nproc) ${i} \
         ${i} | samtools view -Sb -F4 - | \
         samtools sort -@ $(nproc) - -o ${i%%_*}.mapped.sorted.bam 2>/dev/null

assembly_stats_cov() {

  run_name=$(basename `pwd` | cut -d\_ -f1)
  mkdir -p OUTPUT ASSEMBLY/Stats RESULTS
  echo -e "Assembly\tNumber_of_contigs\tCoverage\tAssembly_size\tLargest_contig\tN50\tN90" \
      > assembly_$run_name.stats

  for i in ASSEMBLY/*assembly.fa
  do
    assembly=ASSEMBLY/$(basename "$i" | cut -d\- -f3 --complement)
    name=$(echo $assembly | cut -d\/ -f2)
    reads=TRIMMING/$(basename $i | cut -d\- -f1)

    echo -e "Assembly:\t${name}" | tee ${name}.stats
    bwa index -a bwtsw $i -p ${name}
    bwa mem -t $(nproc) ${name} \
         ${reads}_R1.trim.fastq.gz ${reads}_R2.trim.fastq.gz | \
         samtools view -Sb -F4 - | \
         samtools sort -@ $(nproc) - -o ${name}.mapped.sorted.bam 2>/dev/null
    samtools index ${name}.mapped.sorted.bam
    rm ${name}.{bwt,pac,ann,amb,sa}

    cov=$(samtools mpileup ${name}.mapped.sorted.bam 2> /dev/null | \
          awk '{ count++ ; sum += $4 } END { printf "%s\t%.2f\n", "Coverage:", sum/count }')
    cover=$(echo $cov | awk '{ print $2 }')
    cat $i | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ \
        { print n $0; n = "" } END { printf "%s", n }'| \
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
