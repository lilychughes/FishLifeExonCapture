#!/usr/bin/env python

import re
from sys import argv
import argparse
import os

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser(description="Requires python 2.7 and Biopython. Filters the longest, highest coverage exon from the exonerate output, and renames the sequence simply as >Taxon")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file')
parser.add_argument('-t', '--taxon', dest = 'taxon', type = str, default = None, required = True, help = 'Name of taxon')
parser.add_argument('-l', '--length', dest = 'length', type = int, default = 100, required = True, help = 'Minimum sequence length to write to file. Default = 100 bp')
parser.add_argument('-m', '--mito', dest = 'mito', type = str, default = 'False', required = False , help = 'If set to True, uses vertebrate mitochondiral genetic code. Default = False')

args, unknown = parser.parse_known_args()


# read in sequences
input = open(args.fasta)
records = list(SeqIO.parse(input, "fasta"))
input.close()
records.sort(key=lambda r: -len(r))

# filter just the sequences in the correct orientation (known from the reference sequence used with exonerate)
oriented = []

for record in records:
	coordinates = record.description.split("\t")[1]
	start = int(coordinates.split("-")[0])
	end = int(coordinates.split("-")[1])
	if end > start:
		oriented.append(record)	

# filter for length
longest = []

if len(oriented) > 0:
	for record in oriented:
		if len(record.seq) >= args.length:
			longest.append(record)

longest.sort(key=lambda r: -len(r))

filteredSeq = SeqRecord(longest[0].seq, id=args.taxon, description='')
SeqIO.write(filteredSeq, args.output, "fasta")
