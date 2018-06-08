#!/bin/sh

# this only needs to be run once in the entire pipeline
# uses trimmed reads
# this assumes mixed end data. change if you're using single-end or two paired-end files. see aTRAM docs.

rename .trimmed.fq .trimmed.fastq *.fq

for f in *.trimmed.fastq;
do
	if [  -e ${f%.*}.atram_preprocessor.log  ];
	    then echo Reads already preprocessed! See ${f%.*}.atram_preprocessor.log;
	else    
        atram_preprocessor.py -b ${f%.*} --mixed-ends $f;
    fi;    
done
