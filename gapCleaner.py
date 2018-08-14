#!/usr/bin/env python

import re
from sys import argv
import argparse
import os
from Bio import AlignIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="Requires python 2.7 and Biopython. Removes columns with single taxon insertions.")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file')
args, unknown = parser.parse_known_args()

# read in the alignment in fasta format
alignment = AlignIO.read(args.fasta,"fasta")

# get the number of columns in the alignment
length = alignment.get_alignment_length()

# find the columns that are not composed of mostly gaps (where three or fewer taxa have sequence data or 'NNNN' sequences)
goodColumns = []

for x in range(0,length):
	column = alignment[:,x]
	if column.count("-") < (len(alignment)-3):
		slice = alignment[:,x:x+1]
		goodColumns.append(slice)

# make an empty alignment to populate
		
newAlignment = alignment[:,0:0]

for column in goodColumns:
	newAlignment = newAlignment+column

# change to uppercase

upperAlignment = []

for sequence in newAlignment:
	seqUP = sequence.seq.upper()
	recUP = SeqRecord(seqUP, id=sequence.id, description='')
	upperAlignment.append(recUP)

# write the cleaned alignment to a new file

out = open(args.output, "w")
for item in upperAlignment:
	out.write(item.format("fasta"))
out.close()