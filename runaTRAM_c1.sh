#!/bin/sh



# all the things aTRAM needs to run in c1 environment
module load trinity
module load sqlite
module load aTRAM/2.0
module load blast+

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
	if [  ! -e $directory.step4.aTRAM.txt  ];
	then
		echo $directory aTRAM assembly started > $directory.step4.aTRAM.txt
		cd $directory;
		cp $directory.rmdup.fastq /scratch/;
		cp *fa /scratch/;
		cd /scratch/;
		mkdir temp;
		# if you want to change atram_preprocessor.py options, change here:
		atram_preprocessor.py -b $directory --mixed-ends *fastq -t /scratch/temp/ --fastq;
		for f in *.initial.combined.fa;
		do
			if [  ! -e ${f%.*.*.*.*.*.*.*.*}.${f%.*}.step4.aTRAM.log  ];
			then 
			# if you want to change atram.py options, change here:
			atram.py -b $directory -q $f -a trinity -o trinity -i 5;
			mv *fasta $BASEDIR/$directory/;
			fi;
		done;
		mv *log $BASEDIR/$directory/;
		cd 	$BASEDIR/$directory/;
		rm /scratch/*;
		rm -r /scratch/temp/
		cd ../;
		echo aTRAM assembly completed $directory >> $directory.step4.aTRAM.txt;
	fi;	
fi;
done	
