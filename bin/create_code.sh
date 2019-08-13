#!/bin/bash

grep < $1 \> | sed 's/ /;/' | \
  tr ';' '\t' | awk -v OFS='\t' '{ print $5, $4}' | \
  sed 's/Agly/AGly/g' | tee >( sort -k1,2 \
  > antibiotics_code.v3.tsv) >( cut -f1 | sort | \
  uniq > categories.txt) > /dev/null

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

for i in $(cat categories.txt)
do
  for k in "${samples[@]}"
  do
     key="${k%%:*}"
     val="${k##*:}"

     if [ ${i} == ${key} ] ; then
       sed -i "s/${i}/${val}/" antibiotics_code.v3.tsv
     fi
  done
done

rm categories.txt
