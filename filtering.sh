#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lilychughes@gwu.edu
#SBATCH --job-name=filter
#SBATCH -n 1
#SBATCH -p defq
#SBATCH --output=filter.out
#SBATCH --error=filter.err
#SBATCH -t 14-00:00:00

module load exonerate
module load mafft
module load cd-hit


for directory in *;
do
	if [  -d $directory  ];
	then
		if  [  -e $directory.filtering.started.txt  ];
		then
		echo Filtering of $directory already started
		else
		echo Filtering of $directory started > $directory.filtering.started.txt;
		cd $directory
			rename .trimmed_ _ *
			for f in *filtered_contigs.fasta;
			do
			cd-hit-est -i $f -o $f.cdhit -c 0.99;
			exonerate --model coding2coding -q $f.cdhit -t ../../ReadingFrames/$directory.fasta --ryo ">%qi%qd\n%qas\n" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 2 > $f.exonerate.fasta;
			sed -i 's/-- completed exonerate analysis//g' $f.exonerate.fasta;
			python ../../StopCodonFilter.py -f $f.exonerate.fasta -o ${f%.*}.exonerate_filtered.fasta;
			taxon=${f%_*_*_*_*.*.*.*.*.*.*.*}
			python ../../filterExons.py -f ${f%.*}.exonerate_filtered.fasta -o $f.fa -t $taxon  -e $directory -l 100 -c 1.75;
			done;
		cat *.fa > $directory.filtered.fasta;
		sed -i 's/|E.*$//g' $directory.filtered.fasta;
		perl ../../tx.pl -i $directory.filtered.fasta -p F -o $directory.tx;
		cd ../
		echo Filtering of $directory completed > $directory.filtering.completed.txt;
		fi;
	fi;	 	
done			

