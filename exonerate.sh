#!/bin/sh
cp /storage/home/fishlab/FishLifeData/lily-scripts/FishLifeExonCapture/ReadingFrames/*.fasta .

for f in E*.fasta;
do
	if [  -e exon*${f%.*}*filtered*fasta  ];
	then
		for i in exon*${f%.*}*filtered*fasta;
		do
		exonerate --model coding2coding -q $i -t $f --ryo ">%qi%qd\n%qas\n" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 1 > ${i%.*}.exonerate.fasta;
		done;
	fi;
done;

