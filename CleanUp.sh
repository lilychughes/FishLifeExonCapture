#!/bin/sh

# removes intermediate files
# leaves all raw and trimmed fastq files, output of aTRAM, and filtered exonerate output

for directory in *;
do
if [  -d $directory  ];
then
rm $directory/*log;
rm $directory/*.mapped.bam;
rm -r $directory/*.initial;
rm $directory/*clstr;
rm $directory/*cdhit;
rm $directory/*all_contigs*;
rm $directory/*exonerate.fasta;
rm $directory/*.fq;
rm $directory/*initial.combined.fa
fi;
done

