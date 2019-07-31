#!/usr/bin/env python

import re
from sys import argv
import argparse
import os
from Bio import AlignIO
from Bio import Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser(description="Requires python 2.7 and Biopython. This script cleans alignments for single-taxon insertions, short sequences, gappy edges, but does not consider reading frames. Use the AlignmentCleanerCodon.py script if you want to maintain reading frames starting on the first codon position during trimming. Removal of gappy sequences is controlled with the -c parameter.")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to clean.')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file.')
parser.add_argument('-c', '--coverage', dest = 'coverage', type = float, default = 0.50, required = True, help = 'Maximum proportion of gaps allowed in a sequence. Default = 0.50')
parser.add_argument('-t', '--trim', dest = 'trim', type = float, default = 0.6, required = True, help = 'The maximum proportion of gaps to allow at the edges of the alignment. Edges with more than this proportion of gaps will be trimmed. Setting this value to 1 disables trimming. Default = 0.6')
args, unknown = parser.parse_known_args()

# read in the alignment in fasta format
rawAlignment = AlignIO.read(args.fasta,"fasta")

# find the columns that are not composed of mostly gaps (where two or fewer taxa have sequence data or 'NNNN' sequences)
goodColumns = []

for x in range(0,rawAlignment.get_alignment_length()):
    column = rawAlignment[:,x]
    if column.count("-") < (len(rawAlignment)-2):
        slice = rawAlignment[:,x:x+1]
        goodColumns.append(slice)
# make an empty alignment to populate

goodColumnsAlignment = rawAlignment[:,0:0]
for column in goodColumns:
    goodColumnsAlignment = goodColumnsAlignment+column

# Remove gappy edges where the average is more than 60% gaps. This script does not take

def TrimEdges(alignment):
    """Trims gappy edges of alignments with more than the default 60% gaps given an alignment object"""
    percentages = []
    goodColumnsIndices = []
    for x in range(0,alignment.get_alignment_length()):
        column = alignment[:,x]
        gapPerc = float(column.count("-"))/(len(alignment))
        percentages.append(gapPerc)
    for i in range(0,len(percentages)-1):
        if percentages[i] < args.trim:
           goodColumnsIndices.append(i)
    cleanedAlignment = alignment[:,goodColumnsIndices[0]:goodColumnsIndices[len(goodColumnsIndices)-1]]
    return(cleanedAlignment)

trimmedAlignment = TrimEdges(goodColumnsAlignment)

# write the cleaned alignment to a new file

AlignIO.write(trimmedAlignment, args.output, "fasta")