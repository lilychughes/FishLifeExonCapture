#!/bin/sh

for directory in *;
do
if [  -d $directory  ];
then
if [  ! -e $directory.dedup.txt  ];
then
echo $directory deduplication started > $directory.dedup.txt;
gunzip $directory/*trimmed.fastq.gz;
ls $directory/*.trimmed*
ls $directory/*.trimmed.fastq > $directory/inputlist.txt;
../FastUniq/source/fastuniq -i $directory/inputlist.txt -o $directory/${directory%/}_R1.trimmed_dedup.fastq -p $directory/${directory%/}_R2.trimmed_dedup.fastq;
gzip $directory/*trimmed.fastq;
echo $directory deduplication completed > $directory.dedup.txt;
fi;
fi;
done

