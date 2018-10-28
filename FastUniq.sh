#!/bin/sh

for directory in *;
do
if [  -d $directory  ];
then
cd $directory/;
gunzip *trimmed.fastq;
ls *.trimmed.fastq > inputlist.txt;
../../FastUniq/source/fastuniq -i inputlist.txt -o ${directory%/}_R1.trimmed_dedup.fastq -p ${directory%/}_R2.trimmed_dedup.fastq;
gzip *trimmed.fastq;
cd ../;
echo $directory deduplication completed;
fi;
done