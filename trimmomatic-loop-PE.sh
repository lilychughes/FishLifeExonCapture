#!/bin/sh


###########################
# Put the full path of the Trimmomatic .jar file on your system between the quotes below
###########################
jarfile="/path/to/jarfile.jar"


###########################
# Put the full path of the adapters file on your system between the quotes below
###########################
adapters="/path/to/adapters.fa"




for directory in *;
do
if [  -d $directory  ];
then
	if [  ! -e $directory.step1.trimming.txt  ];
	then
	echo Trimming started >> $directory.step1.trimming.txt;
	cd $directory/;
		for fastq in *R1.fastq.gz;
		do
		java -jar $jarfile PE -threads 4 -phred33 -trimlog $directory.trimlog $fastq ${fastq%_*.*.*}_R2.fastq.gz ${fastq%_*.*.*}_R1.trimmed.fastq.gz ${fastq%_*.*.*}_rem1.fastq.gz ${fastq%_*.*.*}_R2.trimmed.fastq.gz ${fastq%_*.*.*}_rem2.fastq.gz ILLUMINACLIP:$adapters:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:31 > ../$directory.step1.trimming.txt 2>&1;
			if [  -s ${fastq%_*.*.*}_R1.trimmed.fastq.gz  ];
			then
			echo Trimming completed >> ../$directory.step1.trimming.txt;
			else
			echo Trimming failed, please check your Trimmomatic installation >> ../$directory.step1.trimming.txt;
			fi;
		done;
	cd ../;	
	fi;
fi;
done	
