#!/bin/sh

for f in *.atram.fasta;
do
exonerate --model coding2coding --geneticcode 2 -q $f -t Reference_barcode.fasta --ryo ">%qi%qd\n%qas\n" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 1 > ${f%.*}.exonerate.fasta;
done
