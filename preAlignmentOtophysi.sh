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
		cat $directory/*$exon*fa >> Alignments/$exon.unaligned.fasta;
		fi;
	cat $directory/*Mito*fa >> Alignments/G0001.unaligned.fasta
	done;
done < ../FishLifeExonCapture/OtophysiExons.txt

