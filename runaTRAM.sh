#!/bin/sh

# assumes atram is in your path!
# using velvet as the assembler for aTRAM

for f in *.initial.refs.fasta;
do
	if [  -e ${f%.*.*.*.*.*}.*.atram.log  ];
	then 
		echo $f aTRAM assembly already completed;
	else
		atram.py -b ${f%.*.*.*.*.*} -Q $f -a velvet -o exon;
	fi;
done		 	

	
# may add a clean-up step to remove extra aTRAM files
#Pleuronectidae_Glyptocephalus_cynoglossus_KU1474.trimmed.fastq.initial.refs.fasta