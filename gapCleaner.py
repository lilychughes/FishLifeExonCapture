#!/usr/bin/env python

import re
from sys import argv
import argparse
import os
from Bio import AlignIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

parser = argparse.ArgumentParser(description="Requires python 2.7 and Biopython. Removes columns with single taxon insertions by default. If -c is specified, this will also remove sequences with more than this percentage of gaps")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file')
parser.add_argument('-c', '--coverage', dest = 'coverage', type = float, default = 0.50, required = False, help = 'Maximum proportion of gaps allowed. Default = 0.50')
args, unknown = parser.parse_known_args()

# read in the alignment in fasta format
alignment = AlignIO.read(args.fasta,"fasta")

# get the number of columns in the alignment
length = alignment.get_alignment_length()

# get the number of taxa in the alignment
numTaxa = len(alignment)

# find the columns that are not composed of mostly gaps (where three or fewer taxa have sequence data or 'NNNN' sequences)
goodColumns = []

for x in range(0,length):
	column = alignment[:,x]
	if column.count("-") < (numTaxa-3):
		slice = alignment[:,x:x+1]
		goodColumns.append(slice)


# make an empty alignment to populate
		
goodColumnsAlignment = alignment[:,0:0]
for column in goodColumns:
	goodColumnsAlignment = goodColumnsAlignment+column

newlength = goodColumnsAlignment.get_alignment_length()

# drop taxa with more than the allowed percentage of NNNs
goodTaxa = []

for record in goodColumnsAlignment:
	seq = record.seq
	gaps = seq.count("-")
	Ns = seq.count("N")
	if gaps < (newlength*args.coverage) and Ns < 4:
		goodTaxa.append(record)

TaxAlignment =  MultipleSeqAlignment(goodTaxa)

# write the cleaned alignment to a new file

AlignIO.write(TaxAlignment, args.output, "fasta")