#!/bin/sh

for directory in *;
do
if [  -d $directory  ];
then
taxon=$(echo ${directory%/})
echo $taxon;
if [  ! -e $directory.dedup.txt  ];
then
echo $directory deduplication started > $directory.dedup.txt;
cd $directory/;
gunzip *trimmed.fastq.gz;
ls *.trimmed.fastq > inputlist.txt;
../../FastUniq/source/fastuniq -i inputlist.txt -o $taxon_R1.trimmed_dedup.fastq -p $taxon_R2.trimmed_dedup.fastq;
gzip *trimmed.fastq;
cd ../;
echo $directory deduplication completed > $directory.dedup.txt;
fi;
fi;
done