#!/bin/sh

# adapter file needs to be in the same directory for now. Users should change 'adapters.fa' to the appropriate file.

for f in *.fastq;
do
	if [  -e $f.trimlog  ];
		then echo Reads already trimmed.;
	else;
		trimmomatic SE -threads 4 -phred33 -trimlog $f.trimlog $f ${f%.*}.trimmed.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:31;
	fi;
done		