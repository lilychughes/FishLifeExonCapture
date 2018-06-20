#!/bin/sh

# assumes atram is in your path!
# using velvet as the assembler for aTRAM; velvet must also be in the path

for directory in *;
do
if [  -d $directory  ];
then
	if [ ! -e $directory.started.txt  ];
	then
		echo $directory.started.txt found!;
	else	
		cd $directory;
		echo $directory aTRAM assembly started > ../$directory.started.txt
		for f in *.initial.combined.fa;
			do
			if [  ! -e ${f%.*.*.*.*.*.*.*}.${f%.*}.atram.log  ];
				then 
				atram.py -b ${f%.*.*.*.*.*.*.*} -q $f -a velvet -o exon;
			fi;
		cd ../;
		echo aTRAM assembly completed $directory > $directory.completed.txt;
		fi;	
	fi;
done	
	 	

	
# may add a clean-up step to remove extra aTRAM files
