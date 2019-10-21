#!/bin/sh

# assumes atram is in your path!
# using velvet as the assembler for aTRAM; trinity must also be in the path

for directory in *;
do
if [  -d $directory  ];
then
	if [  ! -e $directory.step4.aTRAM.txt  ];
	then
		echo $directory aTRAM assembly started > $directory.step4.aTRAM.txt
		cd $directory;
		gunzip *rmdup.fastq.gz;
		atram_preprocessor.py -b ${directory%/} --mixed-ends *.fastq;
		ls *blast* > preprocess_files.txt;
			if [  ! -s preprocess_files.txt  ];
			then
				echo atram_preprocessor.py is not running correctly. Check aTRAM_preprocessor log files in $directory. >> ../$directory.step4.aTRAM.txt
			else
				for f in *.initial.combined.fa;
				do
					if [  ! -e ${f%.*.*.*.*.*.*.*}.${f%.*}.atram.log  ];
					then 
						atram.py -b ${directory%/} -q $f -a trinity -o trinity -i 5 --cpus 6;
					fi;
				done;
			fi;	
		gzip *fastq;
		ls *filtered_contigs.fasta > atram_list.txt
			if [  ! -s atram_list.txt  ];
			then
				echo aTRAM did not assemble any contigs. Check aTRAM log files in $directory and check all dependencies are loaded (blast, Trinity, samtools) >> ../$directory.step4.aTRAM.txt;
			else 			
				echo aTRAM assembly completed $directory >> ../$directory.step4.aTRAM.txt
			fi;
		cd ../;
	fi;	
fi;
done	