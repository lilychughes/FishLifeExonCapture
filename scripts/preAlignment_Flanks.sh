#!/bin/sh

# make a directory to store alignments
mkdir Alignments_Flanks/


# cat individual exons into files for alignment
while read -r exon;
do
	for directory in *;
	do
		if [  -d $directory  ];
		then
		contigs=$(grep -o ">" $directory/trinity*.$exon.*filtered_flanks.fa | wc -l); 
		if [  $contigs -eq 1  ]
		then
		cat $directory/trinity*.$exon.*filtered_flanks.fa >> Alignments_Flanks/$exon.unaligned.fasta;
		fi;
		fi;
	done;
done < ../FishLifeExonCapture/scripts/ExonList.txt

# cat individual mitochondrial exons into files for alignment
while read -r exon;
do
	for directory in *;
	do
		if [  -d $directory  ];
		then
		contigs=$( grep -c ">" $directory/trinity*.$exon.*filtered_flanks.fa )
		if [  $contigs -eq 1  ]
		then
		cat $directory/trinity*$exon*.filtered_flanks.fa >> Alignments_Flanks/$exon.unaligned.fasta;
		fi;
		fi;
	done;
done < ../FishLifeExonCapture/scripts/MitochondrialExonList.txt

