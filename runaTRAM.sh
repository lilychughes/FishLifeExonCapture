#!/bin/sh

# assumes atram is in your path!
# using velvet as the assembler for aTRAM

for f in *.initial.longest.fasta;
do
	if [  -e ${f%.*.*.*.*.*.*.*}.${f%.*}.atram.log  ];
		then echo aTRAM assembly already complete for $f;
	else
		atram.py -b ${f%.*.*.*.*.*.*.*} -q $f -a velvet -o exon -i 10;
	fi;
done		 	
	
Bothidae_Bothus_lunatus_USNMT154.Bothidae_Bothus_lunatus_USNMT154.fastq.COI.fq.initial.longest.atram.log	