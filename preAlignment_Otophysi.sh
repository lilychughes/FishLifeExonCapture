#!/bin/sh

# make a directory to store alignments
mkdir Alignments/


# cat COI (G0001)
for directory in *;
do
if [  -d $directory  ];
then
cat $directory/trinity*G0001*final_contigs.fa >> Alignments/G0001.unaligned.fasta;
fi;
done

# cat individual exons into files for alignment
while read -r exon;
do
	for directory in *;
	do
		if [  -d $directory  ];
		then
		cat $directory/trinity*.$exon.*final_contigs.fa >> Alignments/$exon.unaligned.fasta;
		fi;
	done;
done < ../FishLifeExonCapture/OtophysiExons.txt



