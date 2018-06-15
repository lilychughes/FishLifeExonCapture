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
				if [  -e ${i%.*}.exonerate.fasta  ];
				then
					echo reading frames already found $i;
				else		
					exonerate --model coding2coding -q $i -t $f --ryo ">%qi%qd\n%qas\n" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 1 > ${i%.*}.exonerate.fasta;
					sed -i 's/-- completed exonerate analysis//g' ${i%.*}.exonerate.fasta;
				fi;
			done;	
		fi;
	done;
done;

