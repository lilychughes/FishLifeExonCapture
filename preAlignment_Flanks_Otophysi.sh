#!/bin/sh

# make a directory to store alignments
mkdir Alignments/


# cat COI (G0001)
for directory in *;
do
if [  -d $directory  ];
then
cat $directory/exon*G0001*filtered_flanks.fa >> Alignments_Flanks/G0001.unaligned.fasta;
fi;
done

# cat individual exons into files for alignment
while read -r exon;
do
	for directory in *;
	do
		if [  -d $directory  ];
		then
		cat $directory/trinity*.$exon.*filtered_flanks.fa >> Alignments_Flanks/$exon.unaligned.fasta;
		fi;
	done;
done < ../FishLifeExonCapture/OtophysiExons.txt



