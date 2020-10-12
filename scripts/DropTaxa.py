#!/usr/bin/env python

import re
from sys import argv
import argparse

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

parser = argparse.ArgumentParser(description="Written under python 2.7 and requires Biopython. Will also work with python 3!")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
parser.add_argument('-t', '--taxa', dest = 'taxa', type = str, default = None, required = True, help = 'List of taxa to remove from alignment, as a text file')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file')
parser.add_argument('-a', '--alignment', dest = 'alignment', type = str, default = "False", required = False, help = 'Add -a True if this is an alignment file in fasta format and you want to remove any gaps caused by dropping sequences.')
args, unknown = parser.parse_known_args()

reference = open(args.fasta)
seq_dict = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))
reference.close()
allTaxa = seq_dict.keys()
toRemove = open(args.taxa)

TaxatoRemove = []

for line in toRemove:
    taxon = line.strip("\n")
    TaxatoRemove.append(taxon)

keepers =[]

for item in allTaxa:
    if item not in TaxatoRemove:
        keepers.append(seq_dict[item])

# clean out gaps created by removing sequences if sequences are part of a multiple sequence alignment
if 'True' in args.alignment:
    rawAlignment = MultipleSeqAlignment(keepers)
    goodColumns = []
    for x in range(0,rawAlignment.get_alignment_length()):
        column = rawAlignment[:,x]
        if column.count("-") < (len(rawAlignment)-2):
            slice = rawAlignment[:,x:x+1]
            goodColumns.append(slice)
    goodColumnsAlignment = rawAlignment[:,0:0]
    for column in goodColumns:
        goodColumnsAlignment = goodColumnsAlignment+column
    AlignIO.write(goodColumnsAlignment, args.output, "fasta" )
# otherwise, just write the sequences
else:
    SeqIO.write(keepers, args.output, "fasta")

