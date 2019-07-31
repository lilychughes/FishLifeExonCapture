#!/usr/bin/env python

import re
from sys import argv
import argparse
import os
from Bio import Seq
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy

parser = argparse.ArgumentParser(description='Requires python 2.7 and Biopython. Flags nucleotide sequences based on distance, based on a user-specified threshold.')
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
parser.add_argument('-d', '--distance', dest = 'distance', type = float , default = 0.6 , required = False, help = 'Distance threshold. Sequences with average distances above this threshold will be flagged. Default = 0.5')
args, unknown = parser.parse_known_args()


# read in the fasta-formatted alignment
alignment = AlignIO.read(args.fasta, "fasta")

# Remove Ns (distance calculator can't deal with them well)

NewSeqs = []

for record in alignment:
	original = str(record.seq.upper()).replace("N","-").replace("R","-").replace("Y","-").replace("W","-").replace("S","-").replace("K","-").replace("M","-")
	new = Seq(original)
	NewRecord = SeqRecord(new, id = record.id, description = '')
	NewSeqs.append(NewRecord)

NewAlignment = MultipleSeqAlignment(NewSeqs)	

# set the method of distance calculation 
calculator = DistanceCalculator('blastn')
# get pairwise distances for all sequences in the alignment
dm = calculator.get_distance(NewAlignment)

# get a list of sequence names
names = []
for sequence in alignment:
	names.append(sequence.id)

# get the mean sequence distance for each sequence
meanDist = []
for x in range(0,len(alignment)):
	mean = (sum(dm[x]))/float(len(alignment))
	meanDist.append(mean)

# make a numpy formatted array
meanArray = numpy.array(meanDist)

# make a dictionary of sequences and their means
dictionary = dict(zip(names,meanDist))

# print sequence names that have more than two standard deviations difference in their mean distance
for item in dictionary.keys():
	if float(dictionary[item]) > args.distance:
		print(item)
		