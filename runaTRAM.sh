#!/bin/sh

# assumes atram is in your path!
# using velvet as the assembler for aTRAM

for f in *.initial.combined.fa;
do
	if [  -e ${f%.*.*.*.*.*.*.*.*}.${f%.*}.atram.log  ];
	then 
		echo aTRAM assembly already complete for $f;
	else
		atram.py -b ${f%.*.*.*.*.*.*.*} -q $f -a velvet -o exon -i 10;
	fi;
done		 	



for f in *.initial.combined.fa;
do
	if [  -e exon.${f%.*.*.*.*.*.*.*}_${f%.*}.filtered_contigs.fasta  ]
	then 
		python2.7 getLongest.py -f exon.${f%.*.*.*.*.*.*.*}_${f%.*}.filtered_contigs.fasta -o ${f%.*}.atram.fasta;
	fi;
done		


	
# may add a clean-up step to remove extra aTRAM files
