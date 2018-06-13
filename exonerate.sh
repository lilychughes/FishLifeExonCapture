#!/bin/sh


for f in /storage/home/fishlab/FishLifeData/lily-scripts/FishLifeExonCapture/ReadingFrames/*.fasta;
do
	for i in *${f%.*}*filtered*fasta;
	do
		if [  -e ${i%.*}.exonerate.fasta  ];
		then 
		echo Reading Frames found for $i;
		else
		exonerate --model coding2coding --geneticcode 2 -q $i -t $f --ryo ">%qi%qd\n%qas\n" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 1 > ${i%.*}.exonerate.fasta;
	done;
done;
