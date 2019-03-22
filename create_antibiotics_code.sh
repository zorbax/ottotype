#!/bin/bash  

grep < ARGannot_r3.fasta \> | sed 's/ /;/' | \
    tr ';' '\t' | awk -v OFS='\t' '{ print $5, $4}' | \
    tee >( sort -k1,2 | sed 's/Agly/AGly/' > code.raw.tsv) \
    >( cut -f1 | sort | uniq > categories.txt) \
    > /dev/null

samples=(
         AGly:Aminoglucósidos
         Agly:Aminoglucósidos
         Bla:Betalactámicos
         Colistin:Colistina
         Fcd:
         Fcyn:Fosfomicina
         Flq:Fluoroquinolonas
         Gly:Glicopéptidos
         MLS:Macrólidos-Lincosamidas-Estreptograminas
         Ntmdz:Nitroimidazol
         Oxzln:Oxazolidinonas
         Phe:Fenicoles
         Rif:Rifampicina
         Sul:Sulfonamidas
         Tet:Tetraciclinas
         Tmt:Trimetroprima
        )

for i in "${samples[@]}"
do
    key="${i%%:*}"
    val="${i##*:}"
while read categorie
do
    if key == categorie
        sed -i 's/categorie/value/' code.raw.tsv
    else
        echo"not foud this gen in samples categories" (print value of $samples[@]) 
    fi
done < categories.txt

fi
=======
#!/bin/bash - 

grep < ARGannot_r3.fasta \> | sed 's/ /;/' | \
    tr ';' '\t' | awk -v OFS='\t' '{ print $5, $4}' | \
    tee >( sort -k1,2 > antibiotics_code.v3.raw.tsv) \
    >( cut -f1 | sort | uniq > antibiotics_categories.txt) \
    > /dev/null

>>>>>>> be4d3252e79051f44768936f8a02fc512587876f
