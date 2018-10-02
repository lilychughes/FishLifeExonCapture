#!/bin/sh



# all the things aTRAM needs to run in c1 environment
module load velvet
module load sqlite
module load aTRAM/2.0

# special python3 environment set up for aTRAM
source $aTRAM



# Altered slightly for Colonial One, because the sqlite databases won't work on the /lustre/ file system
# This just moves these files *temporarily* to the /scratch space on the compute node
# If you need to send the files somewhere other than /scratch, just replace /scratch below


# may add a clean-up step to remove extra aTRAM files

##########################
### MUST BE RUN IN SAMPLES DIRECTORY TO WORK
BASEDIR=$( pwd )
###

### requested by C1 staff
lfs setstripe $BASEDIR -c 1


for directory in *;
do
if [  -d $directory  ];
then
	if [  ! -e $directory.aTRAM.txt  ];
	then
		echo $directory aTRAM assembly started > $directory.aTRAM.txt
		cd $directory;
		cp $directory.trimmed.fastq /scratch/
		cp *fa /scratch/
		cd /scratch/
		# if you want to change atram_preprocessor.py options, change here:
		atram_preprocessor.py -b $directory --mixed-ends $directory.trimmed.fastq;
		for f in *.initial.combined.fa;
		do
			if [  ! -e ${f%.*.*.*.*.*.*.*.*}.${f%.*}.atram.log  ];
			then 
			# if you want to change atram.py options, change here:
			atram.py -b ${f%.*.*.*.*.*.*.*} -q $f -a velvet -o exon -i 10;
			mv *fasta $BASEDIR/$directory/;
			fi;
		done;
		mv *log $BASEDIR/$directory/;
		cd 	$BASEDIR/$directory/;
		rm /scratch/*
		cd ../;
		echo aTRAM assembly completed $directory >> $directory.aTRAM.txt;
	fi;	
fi;
done	
