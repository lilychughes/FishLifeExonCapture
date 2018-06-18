#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lilychughes@gwu.edu
#SBATCH --job-name=atram
#SBATCH -n 1
#SBATCH -p debug
#SBATCH --output=atram.out
#SBATCH --error=atram.err
#SBATCH -t 01:00:00



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


for f in *.db;
do
	if [  -e ${f%.*.*}*atram.log  ];
	then
		echo aTRAM running for ${f%.*.*};
	else
		for i in /lustre/groups/ortilab/FishLife/180327_K00242_0377_BHTGVMBBXX-RB-BF03-L2/${f%.*.*}*.fa;
		do
			atram.py -b ${f%.*.*} -q $i -a velvet -o /lustre/groups/ortilab/FishLife/180327_K00242_0377_BHTGVMBBXX-RB-BF03-L2/exon;
		done
	fi;
done		 
