#!/bin/sh

# velvet assembly of mapped reads by locus
# COI used for testing

# run velvet on the reads that mapped to COI and pull out the longest contig that was assembled
# this longest contig will be the input to aTRAM

for f in *COI.fq;
do
velveth $f.initial 29 -short -fastq $f;
velvetg $f.initial;
python getLongest.py -f $f.initial/contigs.fa -o $f.initial.longest.fasta;
done

