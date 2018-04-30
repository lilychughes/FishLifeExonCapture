#!/bin/sh

for f in *.atram.fasta;
do
/storage/home/fishlab/FishLifeData/lily-scripts/exonerate-2.2.0/bin/exonerate --model coding2coding --geneticcode 2 -q $f -t Reference_barcode.fasta --ryo ">%qi%qd\n%qas\n" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 1 > ${f%.*}.exonerate.fasta;
done
