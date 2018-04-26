#!/bin/sh

# velvet assembly of mapped reads by locus
# COI used for testing

# run velvet on the reads that mapped to COI and pull out the longest contig that was assembled
# this longest contig will be the input to aTRAM
# getLongest.py needs to be in the path


# assemble mapped reads into a starter contig using velvet, and get the longest contig
for f in *COI.fq;
do
	if [  -e $f.initial.longest.fasta  ];
		then echo Initial assembly of $f reads already completed.;
	else
		velveth $f.initial 29 -short -fastq $f;
		velvetg $f.initial;
		python2.7 getLongest.py -f $f.initial/contigs.fa -o $f.initial.longest.fasta;		
	fi;
done


# remove other velvet output (not used in future steps) if fasta file exists
for f in *.initial;
do
	if [  -e $f.longest.fasta  ];
		then rm -r $f;
	fi;
done