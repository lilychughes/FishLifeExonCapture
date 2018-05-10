#!/bin/sh

# adapter file needs to be in the same directory for now. Users should change 'adapters.fa' to the appropriate file.
# currently points to where trimmomatic is installed on Makaria; this would need to be changed for other systems.

for f in *.fastq;
do
	if [  -e $f.trimlog  ];
		then echo Reads already trimmed.;
	else
		java -jar /storage/apps/trimmomatic/0.36/trimmomatic-0.36.jar SE -threads 4 -phred33 -trimlog $f.trimlog $f ${f%.*}.trimmed.fq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:31;

	fi;
done		