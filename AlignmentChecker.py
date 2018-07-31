#!/usr/bin/env python

import re
from sys import argv
import argparse
import os
from Bio import Seq
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
import numpy

parser = argparse.ArgumentParser(description='Requires python 2.7 and Biopython. Flags sequences based on distance')
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
args, unknown = parser.parse_known_args()


# read in the fasta-formatted alignment
alignment = AlignIO.read(args.fasta, "fasta")

# set the method of distance calculation 
calculator = DistanceCalculator('identity')
# get pairwise distances for all sequences in the alignment
dm = calculator.get_distance(alignment)

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
	if float(dictionary[item]) > ((numpy.mean(meanArray))+(numpy.mean(meanArray)+((numpy.std(meanArray))*2))):
		print(item)
		