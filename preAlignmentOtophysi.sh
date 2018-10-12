#!/bin/sh

# make a directory to store alignments
mkdir Alignments/


# cat individual exons into files for alignment
while read -r exon;
do
	for directory in *;
	do
		if [  -d $directory  ];
		then
		cat $directory/exon*$exon*filtered.fa >> Alignments/$exon.unaligned.fasta;
		cat $directory/exon*Mito*filtered.fa >> Alignments/G0001.unaligned.fasta
		fi;	
	done;
done < ../FishLifeExonCapture/OtophysiExons.txt

