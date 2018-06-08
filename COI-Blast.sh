#!/bin/sh

# blasting COI
# uses COI BOLD sequences for Actinopterygii (downloaded April 27th, 2018)

for f in *exonerate.fasta;
do
	if [  -e ${f%.*}.txt  ];
		then echo Blast search already completed!;
	else
		blastn -db BOLDFishCOI.fas -query $f -outfmt 6 -out ${f%.*}.txt -max_target_seqs 10;
	fi;
done		