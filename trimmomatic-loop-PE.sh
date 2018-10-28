#!/bin/sh

for directory in *;
do
if [  -d $directory  ];
then
	if [  ! -e $directory.trimming.txt  ];
	then
	echo Trimming started > $directory.trimming.txt;
	cd $directory/;
	for fastq in *R1.fastq.gz;
	do
#### If you are running this on a different system, change the path to the trimmomatic jar file and adapters fasta file below
	java -jar /storage/apps/trimmomatic/0.36/trimmomatic-0.36.jar PE -threads 4 -phred33 -trimlog $directory.trimlog $fastq ${fastq%_*.*.*}_R2.fastq.gz ${fastq%_*.*.*}_R1.trimmed.fastq.gz ${fastq%_*.*.*}_rem1.fastq.gz ${fastq%_*.*.*}_R2.trimmed.fastq.gz ${fastq%_*.*.*}_rem2.fastq.gz ILLUMINACLIP:../adapters.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:31;
	echo Trimming completed > ../$directory.trimming.txt;
	done;
	cd ../;	
	fi;
fi;
done	
