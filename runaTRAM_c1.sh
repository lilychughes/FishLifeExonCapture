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
		gunzip *rmdup.fastq.gz;
		cp $directory.rmdup.fastq /scratch/;
		cp *fa /scratch/;
		cd /scratch/;
		mkdir temp;
		atram_preprocessor.py -b ${directory%/} --mixed-ends *.fastq -t /scratch/temp/;
		ls *blast* > preprocess_files.txt;
			if [  ! -s preprocess_files.txt  ];
			then
				echo atram_preprocessor.py is not running correctly. Check aTRAM_preprocessor log files in $directory. >> $BASEDIR/$directory.step4.aTRAM.txt
			else
				for f in *.initial.combined.fa;
				do
					if [  ! -e ${f%.*.*.*.*.*.*.*}.${f%.*}.atram.log  ];
					then 
						atram.py -b ${directory%/} -q $f -a trinity -o trinity -i 5 --cpus 6;
						mv *fasta $BASEDIR/$directory/;
						mv *log $BASEDIR/$directory/;
						rm /scratch/*;
						rm -r /scratch/temp/;
					fi;
				done;
			fi;	
		cd $BASEDIR/$directory;	
		gzip *fastq;
		ls *filtered_contigs.fasta > atram_list.txt
			if [  ! -s atram_list.txt  ];
			then
				echo aTRAM did not assemble any contigs. Check aTRAM log files in $directory and check all dependencies are loaded (blast, Trinity, samtools) >> $BASEDIR/$directory.step4.aTRAM.txt;
			else 			
				echo aTRAM assembly completed $directory >> $BASEDIR/$directory.step4.aTRAM.txt;
			fi;
		cd $BASEDIR;
	fi;	
fi;
done	