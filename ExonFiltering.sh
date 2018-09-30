#!/bin/sh


#Filtering with CD-HIT & Exonerate


for directory in *;
do
	
	if  [  -d $directory  ];
	then
		if [  ! -e $directory.exonfiltering.txt  ];
		then
		echo Filtering exons started $directory > $directory.exonfiltering.txt;
		cd $directory;
			
			for f in *filtered_contigs.fasta;
			do
			# change the similarity threshold here at -c if you want
			cd-hit-est -i $f -o $f.cdhit -c 0.98;
			done;
			
			while read -r exon;
			do
			filterfile = $( echo *$exon*.cdhit );
			exonerate --model coding2genome -t $filterfile -q ../../ReadingFrames/$exon.fasta --ryo ">%ti\t%qab-%qae\n%tas" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 2 > $filterfile.exonerate.fasta;
			sed -i 's/-- completed exonerate analysis//g' $filterfile.exonerate.fasta;
			done < ../../FishLifeExonCapture/ExonList.txt
			
			for f in *exonerate.fasta;
			do
			taxon=${f%_*_*_*_*.*.*.*.*.*.*.*}
			python ../../filterExons.py -f $f.exonerate.fasta -o $f.exonerate_filtered.fa -t $taxon -l 100 -c 1.5;
			done;
			
		cd ../;
		echo Filtering exons completed $directory > $directory.exonfiltering.txt;	
		fi;
	fi;	
done	