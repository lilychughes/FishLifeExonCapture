#!/bin/sh



# all the things aTRAM needs to run in the Pegasus environment
module load trinity
module load samtools
module load blast+

# aTRAM needs to be locally installed on Pegasus, and added to your path
# I have a miniconda3 installation, and it seems to work fine for installing dependencies of aTRAMup


# Altered slightly for Pegasus, because the sqlite databases won't work on the /lustre/ file system
# This just moves these files *temporarily* to the /local space on the compute node
# If you need to send the files somewhere other than /local, just replace /local below


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
		gunzip *rmdup.fastq.gz;
		cp $directory.rmdup.fastq /local/;
		cp *fa /local/;
		cd /local/;
		atram_preprocessor.py -b ${directory%/} --mixed-ends *.fastq;
		ls *blast* > preprocess_files.txt;
			if [  ! -s preprocess_files.txt  ];
			then
				echo atram_preprocessor.py is not running correctly. Check aTRAM_preprocessor log files in $directory. >> $BASEDIR/$directory.step4.aTRAM.txt
			else
				for f in *.initial.combined.fa;
				do
						atram.py -b ${directory%/} -q $f -a trinity -o trinity -i 5 --cpus 6;
						mv *fasta $BASEDIR/$directory/;
						mv *log $BASEDIR/$directory/;
				done;
			fi;	
		rm /local/*;
		cd $BASEDIR/$directory;	
		gzip *fastq;
		ls *filtered_contigs.fasta > atram_list.txt
			if [  ! -s atram_list.txt  ];
			then
				echo aTRAM did not assemble any contigs. Check aTRAM log files in $directory and check all dependencies are loaded >> $BASEDIR/$directory.step4.aTRAM.txt;
			else 			
				echo aTRAM assembly completed $directory >> $BASEDIR/$directory.step4.aTRAM.txt;
			fi;
		cd $BASEDIR;
	fi;	
fi;
done	