#!/bin/sh


cp /storage/home/fishlab/FishLifeData/lily-scripts/FishLifeExonCapture/ReadingFrames/*.fasta .

for f in E*.fasta;
do
	for i in exon*${f%.*}*filtered*fasta;
	do
		if [  -e ${i%.*}.exonerate.fasta  ];
		then 
		echo Reading Frames found for $i;
		else
		exonerate --model coding2coding -q $i -t $f --ryo ">%qi%qd\n%qas\n" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 1 > ${i%.*}.exonerate.fasta;
		fi;
	done;
done;

rm E*.fasta