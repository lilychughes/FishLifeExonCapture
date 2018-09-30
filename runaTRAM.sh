#!/bin/sh

# assumes atram is in your path!
# using velvet as the assembler for aTRAM; velvet must also be in the path

for directory in *;
do
if [  -d $directory  ];
then
	if [  ! -e $directory.aTRAM.txt  ];
	then
		echo $directory aTRAM assembly started > ../$directory.aTRAM.txt
		cd $directory;
		atram_preprocessor.py -b $directory --mixed-ends $directory.trimmed.fastq;
		for f in *.initial.combined.fa;
			do
				if [  ! -e ${f%.*.*.*.*.*.*.*}.${f%.*}.atram.log  ];
					then 
					atram.py -b ${f%.*.*.*.*.*.*} -q $f -a velvet -o exon -i 10;
				fi;
			done;
		cd ../;
		echo aTRAM assembly completed $directory >> $directory.aTRAM.txt;
	fi;	
fi;
done	
	 	

	
# may add a clean-up step to remove extra aTRAM files
