#!/bin/sh

# all the things aTRAM needs to run
module load velvet
module load sqlite
module load parallel
module load aTRAM/2.0

# special python3 environment set up for aTRAM
source $aTRAM



# Altered slightly for Colonial One, because the sqlite databases won't work on the /lustre/ file system
# This just moves these files *temporarily* to the /scratch space on the compute node
# If you need to send the files somewhere other than /scratch, just replace /scratch below


# may add a clean-up step to remove extra aTRAM files


mkdir completed_database_files

for f in *.db;
do
	if [  -e ${f%.*.*}*atram.log  ];
	then
		echo aTRAM running for ${f%.*.*};
	else
		mv $f /scratch;
		mv ${f%.*.*}*blast* /scratch;
		for i in ${f%.*.*}*.fa;
		do
			atram.py -b /scratch/${f%.*.*} -q $i -a velvet -o exon;
			mv /scratch/${f%.*.*.*}.trimmed.${i%.*}.atram.log .;
		done
		mv /scratch/$f completed_database_files/;
		mv /scratch/${f%.*.*}*blast* completed_database_files/;
	fi;
done		 
