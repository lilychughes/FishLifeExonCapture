#!/bin/sh
cp /storage/home/fishlab/FishLifeData/lily-scripts/FishLifeExonCapture/ReadingFrames/*.fasta .

for j in *.db;
do
	for f in E*.fasta;
	do	
		if [  -e exon.${j%.*.*.*}_${j%.*.*.*}.trimmed.fastq.${f%.*}.fq.initial.combined.filtered_contigs.fasta  ];
		then
			for i in exon*${f%.*}*filtered_contigs.fasta;
			do
				exonerate --model coding2coding -q $i -t $f --ryo ">%qi%qd\n%qas\n" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 1 > ${i%.*}.exonerate.fasta;
			done;	
		fi;
	done;
done;

