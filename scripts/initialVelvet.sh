#!/bin/sh

# velvet assembly of mapped reads by locus
# this longest contig will be the input to aTRAM
# assemble mapped reads into a starter contig using velvet, and get the longest contig
# deletes empty fq files where no reads were recovered


for directory in *;
do
	if  [  -d $directory  ];
	then
		if  [  ! -e $directory.step3.initialVelvet.txt  ];
		then
		echo Starting initial assemblies > $directory.step3.initialVelvet.txt;
		cd $directory/;
			for f in *.fq;
			do
				if [  -s $f  ];
				then
				velveth $f.initial 29 -short -fastq $f;
				velvetg $f.initial;
				python ../../FishLifeExonCapture/scripts/getLongest.py -f $f.initial/contigs.fa -o $f.initial.combined.fa;
				else
				rm $f;
				fi;
			done;
			ls -d *initial > velvet_directories.txt;
			ls *.initial.combined.fa > velvet_output.txt;
				if [  ! -s velvet_directories.txt  ];
				then
				echo Initial assemblies failing. Check that Step 2 completed successfully & velevet installation >> ../$directory.step3.initialVelvet.txt;
				elif [  ! -s velvet_output.txt  ];
				then
				echo Initial assemblies failing. Check your Biopython installation. >> ../$directory.step3.initialVelvet.txt;
				else
				echo Completed initial assemblies >> ../$directory.step3.initialVelvet.txt;	
				fi;
		cd ../;
		fi;		
	fi;
done

