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
		contigs=$(grep -o ">" $directory/trinity*.$exon.*final_contigs.fa | wc -l); 
		if [  $contigs -eq 1  ]
		then
		cat $directory/trinity*.$exon.*final_contigs.fa >> Alignments/$exon.unaligned.fasta;
		fi;
		fi;
	done;
done < ../FishLifeExonCapture/ExonList.txt

# cat individual mitochondrial exons into files for alignment
while read -r exon;
do
	for directory in *;
	do
		if [  -d $directory  ];
		then
		contigs=$( grep -c ">" $directory/trinity*.$exon.*final_contigs.fa )
		if [  $contigs -eq 1  ]
		then
		cat $directory/trinity*.$exon.*final_contigs.fa >> Alignments/$exon.unaligned.fasta;
		fi;
		fi;
	done;
done < ../FishLifeExonCapture/MitochondrialExonList.txt

