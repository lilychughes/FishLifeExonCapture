#!/usr/bin/env python

import re
from sys import argv
import argparse
import os

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="Requires python 2.7 and Biopython. Filters the longest, highest coverage exon from the exonerate output, and renames the sequence simply as >Taxon|Exon")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file')
parser.add_argument('-t', '--taxon', dest = 'taxon', type = str, default = None, required = True, help = 'Name of taxon')
parser.add_argument('-e', '--exon', dest = 'exon', type = str, default = None, required = True, help = 'Name of locus')
parser.add_argument('-l', '--length', dest = 'length', type = int, default = 100, required = True, help = 'Minimum sequence length to write to file')
parser.add_argument('-c', '--coverage', dest = 'coverage', type = float, default = 5.0, required = True, help = 'Minimum coverage required to write a sequence to a file')

args, unknown = parser.parse_known_args()

seqID = args.taxon+"|"+args.exon

input = open(args.fasta)
records = list(SeqIO.parse(input, "fasta"))
input.close()
records.sort(key=lambda r: -len(r))

if len(records) == 1 and len(records[0].seq) >= args.length:
	seq0 = records[0]
	desc0 = seq0.description.split(" ")[1]
	coverage0 = (float(desc0.split("_")[6]))
	if coverage0 > args.coverage:
		filteredSeq = SeqRecord(records[0].seq, id = seqID, description = records[0].description)
		SeqIO.write(filteredSeq, args.output, "fasta")
		print("Filtered exon for "+args.taxon+" "+args.exon)
	else:
		print("No sequences passed filters for "+args.taxon+" locus "+args.exon)
elif len(records) == 1 and len(records[0].seq) < args.length:
	print("No sequences passed filters for "+args.taxon+" locus "+args.exon)
elif len(records) > 1:		
	seq0 = records[0]
	seq1 = records[1]
	desc0 = seq0.description.split(" ")[1]
	desc1 = seq1.description.split(" ")[1]
	coverage0 = (float(desc0.split("_")[6]))
	coverage1 = (float(desc1.split("_")[6]))
	if len(seq0.seq) < args.length and len(seq1.seq) < args.length:
		print("No sequences passed filters for "+args.taxon+" locus "+args.exon)
	elif coverage0 < args.coverage and coverage1 < args.coverage:
		print("No sequences passed filters for "+args.taxon+" locus "+args.exon)	
	elif coverage0 >= coverage1 and len(seq0.seq) >= len(seq1.seq):
		filteredSeq = SeqRecord(seq0.seq, id = seqID, description = seq0.description)
		SeqIO.write(filteredSeq, args.output, "fasta")
		print("Filtered exon for "+args.taxon+" "+args.exon)
	elif coverage1 > coverage0 and len(seq1.seq) >= len(seq0.seq):
		filteredSeq = SeqRecord(seq1.seq, id = seqID, description = seq1.description)
		SeqIO.write(filteredSeq, args.output, "fasta")
		print("Filtered exon for "+args.taxon+" "+args.exon)
	elif coverage0 < coverage1 and len(seq0.seq) > len(seq1.seq):
		filteredSeq = SeqRecord(seq0.seq, id = seqID, description = seq0.description)
		SeqIO.write(filteredSeq, args.output, "fasta")
		print("Filtered exon for "+args.taxon+" "+args.exon)
