#!/usr/bin/env python

import re
from sys import argv
import argparse
import os

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="Requires python 2.7 and Biopython. Filters the exonerate output for sequences on the correct strand, and renames the sequence simply as >Taxon")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file')
parser.add_argument('-t', '--taxon', dest = 'taxon', type = str, default = None, required = True, help = 'Name of taxon')

args, unknown = parser.parse_known_args()


# read in sequences
input = open(args.fasta)
records = list(SeqIO.parse(input, "fasta"))
input.close()

# filter just the sequences in the correct orientation (known from the reference sequence used with exonerate)
oriented = []

for record in records:
	coordinates = record.description.split("\t")[1]
	start = int(coordinates.split("-")[0])
	end = int(coordinates.split("-")[1])
	if end > start:
		oriented.append(record)	

filtered = []

for record in oriented:
	filteredSeq = SeqRecord(record.seq, id=args.taxon, description='')
	filtered.append(filteredSeq)

SeqIO.write(filtered, args.output, "fasta")
