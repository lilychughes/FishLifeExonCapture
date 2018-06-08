#!/bin/sh

# assumes atram is in your path!
# using velvet as the assembler for aTRAM

for f in *.initial.combined.fa;
do
	if [  -e ${f%.*.*.*.*.*.*.*}.${f%.*}.atram.log  ];
	then 
		echo $f aTRAM assembly already completed;
	else
		atram.py -b ${f%.*.*.*.*.*.*.*} -q $f -a velvet -o exon;
	fi;
done		 	

	
# may add a clean-up step to remove extra aTRAM files
