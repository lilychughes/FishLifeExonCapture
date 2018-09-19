#!/bin/sh

# this only needs to be run once in the entire pipeline
# uses trimmed reads
# this assumes mixed end data. change if you're using single-end or two paired-end files. see aTRAM docs.


for directory in *;
do
if [  -d $directory  ];
then
	if [  -e $directory/$directory.trimmed.atram_preprocessor.log  ];
	    then echo Reads already preprocessed! See ${f%.*}.atram_preprocessor.log;
	else    
        cd $directory;
        atram_preprocessor.py -b $directory --mixed-ends $directory.trimmed.fastq;
        cd ../;
    fi;
fi;        
done
