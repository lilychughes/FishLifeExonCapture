#!/usr/bin/env python

import re
from sys import argv
import argparse
import os

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="Requires python 2.7 and Biopython. Filters the longest, highest coverage exon from the exonerate output, and renames the sequence simply as >Taxon")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file')
parser.add_argument('-t', '--taxon', dest = 'taxon', type = str, default = None, required = True, help = 'Name of taxon')
parser.add_argument('-l', '--length', dest = 'length', type = int, default = 100, required = True, help = 'Minimum sequence length to write to file')
parser.add_argument('-c', '--coverage', dest = 'coverage', type = float, default = 5.0, required = True, help = 'Minimum coverage required to write a sequence to a file')

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
	start = coordinates.split("-")[0]
	end = coordinates.split("-")[1]
	if end > start:
		oriented.append(record)	

# filter for stop codons 

noStops = []

for record in oriented:
	if len(record.seq) % 3 == 0 and "*" not in record.seq.translate():
		noStops.append(record)


# filter for length
longest = []

if len(noStops) > 0:
	for record in noStops:
		if len(record.seq) >= args.length:
			longest.append(record)

# filter for coverage (assumes Velvet-style header)

covered = []

if len(longest) > 0:
	for record in longest:
		seqCoverage = float(record.id.split("_")[6])
		if seqCoverage >= args.coverage:
			covered.append(record)


# write the output 
# if no sequences passed filters, print message to screen
if len(covered) == 0:
	print("No sequences passed filters in "+args.fasta) 			
# print the longest remaining sequence
elif len(covered) >= 1:
	filteredSeq = SeqRecord(covered[0].seq, id=args.taxon, description='')
	SeqIO.write(filteredSeq, args.output, "fasta")
