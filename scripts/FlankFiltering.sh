#!/bin/sh


#Filtering loci with flanks 
#You have to complete step 5 before running this, because it relies on the output of step 5
#This step is optional


for directory in *;
do
	
	if  [  -d $directory  ];
	then
		if [  ! -e $directory.step5b.flankfiltering.txt  ];
		then
			echo Filtering exons started $directory > $directory.step5b.flankfiltering.txt;
		cd $directory;
			
			# Get contigs with reading frames & flanking regions in the correct orientation
			for f in *cdhit;
			do
				python ../../FishLifeExonCapture/scripts/filterFlanks.py -f $f.exonerate.fasta -l $f -o $f.flanks -t $directory;
			done	

			
			# Filter again with CD-HIT
			# Uses a slightly more relaxed threshold, because the flanking regions are more variable than the exon
			# The exon -c parameter is set at 0.99
			# For the flanks, it is relaxed to 0.97
			for f in *flanks;
			do
				cd-hit-est -i $f -o ${f%.*}.filtered_flanks.fa -c 0.97;
			done	
			
			# Check which files have less than one, or more than one contig
			# Files with no contigs are removed
			# Files with more than one contig are labeled as 'failed'
			for f in *filtered_flanks.fa;
			do
				count=$(grep -c ">" $f)
				if [  $count -gt 1  ];
				then
				mv $f $f.failed;
				elif [  $count -eq 0  ];
				then
				rm $f;
				fi;
			done

			# Count the number of loci that passed filters
			ls *filtered_flanks.fa > passedflanks.txt;
			ls *filtered_flanks.fa.failed > failedflanks.txt;
			passed=$(grep -c "trinity" passedflanks.txt);
			failed=$(grep -c "trinity" failedflanks.txt);
			
			
		cd ../;
			echo $passed loci passed all filters >> $directory.step5b.flankfiltering.txt;
			echo $failed loci failed filters >> $directory.step5b.flankfiltering.txt;			
			echo Filtering exons completed $directory >> $directory.step5b.flankfiltering.txt;	
		fi;
	fi;	
done	