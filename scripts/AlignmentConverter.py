#!/bin/usr/env python

from sys import argv
import argparse
from Bio import AlignIO
from Bio.Alphabet import Gapped
from Bio import Alphabet

parser = argparse.ArgumentParser(description="Requires python 2.7 or 3 and Biopython 1.48 or higher. Use 'phylip-relaxed' if you have taxon names longer than 10 characters or they will be truncated. ")
parser.add_argument('-f', '--file' , dest = 'file' , type = str , default= None , required= True, help = 'Alignment file')
parser.add_argument('-i', '--infmt' , dest = 'infmt', type = str, default= None, required= True, help = 'Input format of alignment file')
parser.add_argument('-o', '--outfmt' , dest = 'outfmt', type = str, default= None, required= True, help = 'Output format of alignment file')
parser.add_argument('-n', '--outfile', dest = 'outfile', type = str, default=None, required= True, help = 'Name of file to write re-formatted alignment to')
args, unknown = parser.parse_known_args()

input = open(args.file)

output = open(args.outfile, "w")

alignment = AlignIO.parse(input, args.infmt, alphabet=Alphabet.generic_dna)
AlignIO.write(alignment, args.outfile, args.outfmt)

input.close()
output.close()
