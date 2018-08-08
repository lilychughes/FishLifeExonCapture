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


# loop of cd-hit (clustering at 99% identity), exonerate, and a filter for length,strand orientation, and coverage
# puts all passing sequences into locus.filtered.fasta
# locus.filtered.fasta can then be aligned by preferred method.

for directory in *;
do
	if [  -d $directory  ];
	then
		if  [  ! -e $directory.filtering.started.txt  ];
		then
		echo Filtering of $directory started > $directory.filtering.started.txt;
		cd $directory
			rename .trimmed_ _ *
			for f in *filtered_contigs.fasta;
			do
			cd-hit-est -i $f -o $f.cdhit -c 0.99;
			exonerate --model coding2genome -t $f.cdhit -q ../../ReadingFrames/$directory.fasta --ryo ">%ti\t%qab-%qae\n%tas" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 2 > $f.exonerate.fasta;
			sed -i 's/-- completed exonerate analysis//g' $f.exonerate.fasta;
			taxon=${f%_*_*_*_*.*.*.*.*.*.*.*}
			python ../../filterExons.py -f ${f%.*}.exonerate.fasta -o $f.fa -t $taxon -l 100 -c 1.75;
			done;
		cat *.fa > $directory.filtered.fasta;
		cd ../
		echo Filtering of $directory completed > $directory.filtering.completed.txt;
		fi;
	fi;	 	
done