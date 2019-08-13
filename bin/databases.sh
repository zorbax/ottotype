#!/bin/bash

display_usage(){
  echo -e "\nUsage:"
  echo -e "\t$(basename $0) -p "
}

if [ $# -le 1 ]; then
  display_usage
  exit 1
fi

if [ -z "$var" ]; then
  echo "\$var is empty"
else
  echo "\$var is NOT empty"
fi

while getopts ":p:" opt
do
  case $opt in
    p) path="$OPTARG"
    ;;
    ?) echo "Invalid option -$OPTARG" >&2
       exit 1
    ;;
  esac
done



if [ ! -d "$HOME/bin" ]; then
  mkdir -p $HOME/bin

  cat <<EOF >> $HOME/.bashrc

if [ -d "\$HOME/bin" ] ; then
  export PATH="\$HOME/bin:\$PATH"
fi

EOF
else
  perl -e 'exit(!(grep(m{^$ENV{HOME}/bin$},split(":", $ENV{PATH}))) > 0)'
  if [ $? != 0 ]; then
    cat <<EOF >> $HOME/.bashrc

if [ -d "\$HOME/bin" ] ; then
  export PATH="\$HOME/bin:\$PATH"
fi

EOF
  fi
fi

. ~/.bashrc

dbNCBI="$HOME/bin/16S/NCBI.gz"
dbRDP="$HOME/bin/16S/RDP.gz"
dbSILVA="$HOME/bin/16S/SILVA.gz"
dbKraken="$HOME/bin/Kraken2/yggdrasil"
dbKmerfinder="$HOME/bin/KmerFinder_DB/bacteria.organisms.ATGAC"
plasmid_db="$HOME/bin/plasmidid_db/plasmid.complete.nr100.fna"


mkdir -p $HOME/bin/{16S,KmerFinder_DB}

if [ ! -f "$dbNCBI" ]; then
  esearch -db nucleotide -query '33175[BioProject] OR 33317[BioProject]' \
  | efetch -db nuccore -format fasta | gzip -9 > $HOME/bin/16S/NCBI.gz
elif [ ! -f "$dbRDP" ]; then
  wget --no-check-certificate https://rdp.cme.msu.edu/download/current_Bacteria_unaligned.fa.gz
  gunzip -c current_Bacteria_unaligned.fa.gz \
  | bioawk -cfastx '/\(T\)/{print ">" $name " " $comment "\n" toupper($seq)}' | gzip -9 \
  > $HOME/bin/16S/RDP.gz
elif [ ! -f "$dbSILVA" ]; then
  silva="https://www.arb-silva.de/fileadmin/silva_databases/release_32/Exports/"
  wget $silva/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz
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
elif [ ! -f "$dbKmerfinder" ]; then
  dbKmerfinder2="/mnt/disk1/bin/Kmerfinder_DB/bacteria.organisms.ATGAC"
  if [ ! -f "$dbKraken2" ]; then
    echo "The Kmerfinder Database will use ~13GB of disk space during creation"
    echo "and will take a long time to download it."
    read -p "Do you want to proceed? [y/n] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
      KmerFinder_DB=$HOME/bin/KmerFinder_DB
      current_dir=$(pwd)
      server="ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder"
      wget $server/version/latest/{config,bacteria*} -P $KmerFinder_DB
      tar -xzf $KmerFinder_DB/bacteria.tar.gz -C $KmerFinder_DB
      rm $KmerFinder_DB/*.tar.gz && cd $KmerFinder_DB
      if [ ! -f "bacteria.name" ]; then
        find . -type f -name "bacteria.*" -exec mv {} . \;
        rm -rf srv/ && cd $current_dir
      fi
    fi
  fi
elif [ ! -d "$dbKraken" ]; then
  dbKraken2="/mnt/disk1/bin/Kraken2/yggdrasil"
  if [ ! -d "$dbKraken2" ]; then
    echo "The Kraken Database will use ~100GB of disk space during creation"
    echo "and will take a long time to build it."
    read -p "Do you want to proceed? [y/n] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
      current_dir=$(pwd)
      mkdir $HOME/bin/Kraken2/yggdrasil && cd "$_"
      kraken2-build --standard --threads $(nproc) --db yggdrasil
      cd current_dir
    else
      exit 1
    fi
  fi
elif [! -d ""]; then
  cd $HOME/bin/plasmidid_db
  mkdir db_plasmids

  server="ftp://ftp.ncbi.nlm.nih.gov/refseq/release"
  wget -nv -P db_plasmids/ $server/plasmid/plasmid*genomic.fna.gz

  zcat db_plasmids/plasmid.*.genomic.fna.gz | \
       perl -pe 'if(/\>/){s/\n/\t/}; s/\n//; s/\>/\n\>/' | \
       grep "complete" | grep "Salmonella" | grep -v "CDS" | \
       sed 's/\t/\n/' > plasmid.complete.$(date +%F).fna

  cd-hit-est -i plasmid.complete.2019-04-23.fna \
             -o plasmid.complete.nr100.fna \
             -c 1 -T $(nproc) -M 400000
  mv plasmid.complete.nr100.fna $HOME/bin/plasmidid_db/ && cd ..
  rm -rf db_plasmids
fi
