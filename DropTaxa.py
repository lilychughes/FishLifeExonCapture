#!/usr/bin/env python

import re
from sys import argv
import argparse

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


parser = argparse.ArgumentParser(description="Requires python 2.7 and Biopython. ")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
parser.add_argument('-t', '--taxa', dest = 'taxa', type = str, default = None, required = True, help = 'List of taxa to remove from alignment, as a text file')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file')
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


outfile = open(args.output, "w")

for entry in keepers:
    outfile.write(entry.format("fasta"))

outfile.close()
