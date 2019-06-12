#!/bin/sh

# adapter file needs to be in the same directory for now. Users should change 'adapters.fa' to the appropriate file.
# currently points to where trimmomatic is installed on Makaria; this would need to be changed for other systems.

for directory in *;
do
if [  -d $directory  ];
then
	if [  ! -e $directory.step1.trimming.txt  ];
	then
	echo Trimming started > $directory.step1.trimming.txt;
	cd $directory/;
	### If you are running this on a different system, change the path to the trimmomatic jar file and adapters fasta file below
	java -jar /c1/apps/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads 4 -phred33 -trimlog $directory.trimlog $directory.fastq $directory.trimmed.fastq ILLUMINACLIP:../adapters.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:31;
	echo Trimming completed > ../$directory.trimming.txt;
	tar -cvzf $directory.fastq.tar.gz $directory.fastq;
	if [  -s $directory.fastq.tar.gz  ];
	then
	rm $directory.fastq;
	fi;
	echo Untrimmed reads compressed $directory.fastq.tar.gz  > ../$directory.step1.trimming.txt;
	cd ../;	
	fi;
fi;
done	

