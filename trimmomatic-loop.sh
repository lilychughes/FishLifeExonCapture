#!/bin/sh

# adapter file needs to be in the same directory for now. Users should change 'adapters.fa' to the appropriate file.
# currently points to where trimmomatic is installed on Makaria; this would need to be changed for other systems.

for directory in *;
do
if [  -d $directory  ];
then
	if [  ! -e $directory.trimming.txt  ];
	then
	echo Trimming started > $directory.trimming.txt;
	cd $directory/;
	### If you are running this on a different system, change the path to the trimmomatic jar file and adapters fasta file below
	java -jar /c1/apps/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads 4 -phred33 -trimlog $f.trimlog $f ${f%.*}.trimmed.fq ILLUMINACLIP:../adapters.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:31;
	cd ../;
	echo Trimming completed > $directory.trimming.txt;
	fi;
fi;
done	

