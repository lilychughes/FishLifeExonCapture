#!/bin/sh

# velvet assembly of mapped reads by locus
# COI used for testing

# run velvet on the reads that mapped to COI and pull out the longest contig that was assembled
# this longest contig will be the input to aTRAM



# assemble mapped reads into a starter contig using velvet, and get the longest contig
for f in *.fq;
do
	if [  -e $f.initial  ];
		then echo Initial assembly of $f reads already completed.;
	else
		velveth $f.initial 29 -short -fastq $f;
		velvetg $f.initial;
	fi;
done

