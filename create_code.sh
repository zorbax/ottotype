#!/bin/bash

grep < ARGannot_r3.fasta \> | sed 's/ /;/' | \
  tr ';' '\t' | awk -v OFS='\t' '{ print $5, $4}' | \
  sed 's/Agly/AGly/g' | tee >( sort -k1,2 > code.raw.tsv)\
  >( cut -f1 | sort | uniq > categories.txt)\
  > /dev/null

samples=(
         "AGly:Aminoglucósidos" "Bla:Betalactámicos"
         "Colistin:Colistina"   "Fcd:Ácido_Fusídico"
         "Fcyn:Fosfomicina"     "Flq:Fluoroquinolonas"
         "Gly:Glicopéptidos"    "MLS:Macrólidos-Lincosamidas-Estreptograminas"
         "Ntmdz:Nitroimidazol"  "Oxzln:Oxazolidinonas"
         "Phe:Fenicoles"        "Rif:Rifampicina"
         "Sul:Sulfonamidas"     "Tet:Tetraciclinas"
         "Tmt:Trimetroprima"
        )

while read line
do
  col=$(echo $line | cut -f1)
  gene=$(echo $line | cut -f2)
  for k in "${samples[@]}"
  do
    key="${k%%:*}"
    val="${k##*:}"

    if [ ${col} == ${key} ] ; then
      echo -e "${val}\t${gene}"
    fi
  done
done < <(cat code.raw.tsv)
