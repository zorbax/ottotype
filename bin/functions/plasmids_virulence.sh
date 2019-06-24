#!/bin/bash

: <<'END'
#DTU
mkdir plasmidfinder && mv plasmidfinder.* plasmidfinder
ariba prepareref -f /home/dell/bin/ariba/virulencefinder.fa -m /home/dell/bin/ariba/virulencefinder.tsv output_directory
ariba prepareref -f /home/dell/bin/ariba/plasmidfinder.fa -m /home/dell/bin/ariba/plasmidfinder.tsv output_directory
ariba prepareref -f /home/dell/bin/ariba/vfdb.fa -m /home/dell/bin/ariba/vfdb.tsv output_directory

ariba prepareref -f virulencefinder.fa -m virulencefinder.tsv virulencefinder.prepareref
ariba prepareref -f plasmidfinder.fa -m plasmidfinder.tsv plasmidfinder.prepareref
ariba prepareref -f vfdb.fa -m vfdb.tsv vfdb.prepareref

for i in $(ls *.gz | grep -v fastq | cut -d\. -f1 | cut -d\_ -f1 | sort | uniq)
do
  name=`echo $i| cut -d\_ -f1`
  ariba run $HOME/bin/ariba/card/out.card.prepareref $i\_R1.fqtrim.gz $i\_R2.fqtrim.gz $name\_card
done

for i in $(ls *.gz | grep -v fqtrim | cut -d\. -f1 | cut -d\_ -f1 | sort | uniq)
do
  name=`echo $i| cut -d\_ -f1`
  ariba run $HOME/bin/ariba/card/out.card.prepareref $i\_R1.fastq.gz $i\_R2.fastq.gz $name\_card
done

for i in $(ls *.gz | grep -v fastq | cut -d\. -f1 | cut -d\_ -f1 | sort | uniq)
do
  name=`echo $i| cut -d\_ -f1`
  ariba run $HOME/bin/ariba/vfdb/vfdb.prepareref $i\_R1.fqtrim.gz $i\_R2.fqtrim.gz $name\_virulence_vfdb
done

for i in $(ls *.gz | grep -v fqtrim | cut -d\. -f1 | cut -d\_ -f1 | sort | uniq)
do
  name=`echo $i| cut -d\_ -f1`

  ariba run $HOME/bin/ariba/vfdb/vfdb.prepareref $i\_R1.fastq.gz $i\_R2.fastq.gz $name\_virulence_vfdb
done
END
