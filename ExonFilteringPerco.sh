#!/bin/sh


#Filtering with CD-HIT & Exonerate


for directory in *;
do
	
	if  [  -d $directory  ];
	then
		if [  ! -e $directory.step5.exonfiltering.txt  ];
		then
			echo Filtering exons started $directory > $directory.step5.exonfiltering.txt;
		cd $directory;
			
			for f in *filtered_contigs.fasta;
			do
			# change the similarity threshold here at -c if you want
				cd-hit-est -i $f -o $f.cdhit -c 1;
			done;
			
			# Filter nuclear exons with exonerate
			while read -r exon;
			do
				filterfile=$( echo trinity*$exon.*.cdhit );
				exonerate --model coding2genome -t $filterfile -q ../../FishLifeExonCapture/ReadingFramesPercomorph/$exon.fasta --ryo ">%ti\t%qab-%qae\n%tas" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 2 > $filterfile.exonerate.fasta;
				sed -i 's/-- completed exonerate analysis//g' $filterfile.exonerate.fasta;
			done < ../../FishLifeExonCapture/ExonList.txt
			

			# Filter mitochondrial exons with exonerate
			while read -r exon;
			do
				filterfile=$( echo trinity*$exon.*.cdhit );
				exonerate --model coding2genome -t $filterfile -q ../../FishLifeExonCapture/ReadingFramesPercomorph/$exon.fasta --ryo ">%ti\t%qab-%qae\n%tas" --geneticcode 2 --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 2 > $filterfile.exonerateMito.fasta;
			sed -i 's/-- completed exonerate analysis//g' $filterfile.exonerateMito.fasta;
			done < ../../FishLifeExonCapture/MitochondrialExonList.txt
			
			# Get contigs in the correct orientation
        	for f in trinity*.exonerate*fasta;
	    	do
	        	python2.7 ../../FishLifeExonCapture/filterExons.py -f $f -o $f.exonerate_filtered.final_contigs.fa -t $directory;
		    done
			
			
		cd ../;
			echo Filtering exons completed $directory > $directory.step5.exonfiltering.txt;	
		fi;
	fi;	
done	