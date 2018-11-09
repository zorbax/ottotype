#!/bin/bash

dbNCBI="$HOME/bin/16S/NCBI.gz"
dbRDP="$HOME/bin/16S/RDP.gz"
dbSILVA="$HOME/bin/16S/SILVA.gz"

if [ -f "$dbNCBI" ]; then
  esearch -db nucleotide -query '33175[BioProject] OR 33317[BioProject]' \
    | efetch -db nuccore -format fasta | gzip -9 > $HOME/bin/16S/NCBI.gz
elif [ -f "$dbRDP"]; then
  wget --no-check-certificate https://rdp.cme.msu.edu/download/current_Bacteria_unaligned.fa.gz
    gunzip -c current_Bacteria_unaligned.fa.gz \
    | bioawk -cfastx '/\(T\)/{print ">" $name " " $comment "\n" toupper($seq)}' | gzip -9 \
    > $HOME/bin/16S/RDP.gz
elif [ -f "$dbSILVA" ]; then
  wget https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz
  gunzip -v SILVA_132_SSURef_Nr99_tax_silva.fasta.gz \
    | bioawk -cfastx \
    '$comment ~ /^Bacteria;|^Archaea;/ \
    && $comment !~ /(;unidentified|Mitochondria;|;Chloroplast|;uncultured| sp\.)/ \
    { sub(/^.*;/,"",$comment);
      gsub("U","T",$seq);
      print ">" $name " " $comment "\n" $seq }' | seqtk seq -l 60 -U > SILVA.tmp1
   cd-hit-est -i SILVA.tmp1 -o SILVA.tmp2 -c 1.0 -T 0 -M 2000 -d 250
   cp SILVA.tmp2 $HOME/bin/16S/SILVA && rm -f SILVA.tmp1 SILVA.tmp2 SILVA.tmp2.clstr
   cd $HOME/bin/16S/SILVA && gzip -9 SILVA && cd $HOME
fi
