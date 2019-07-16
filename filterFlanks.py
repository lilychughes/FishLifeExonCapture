#!/usr/bin/env python

import re
from sys import argv
import argparse
import os
import sys

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="Requires python 2.7 and Biopython. Filters the exonerate output for sequences on the correct strand, and renames the sequence simply as >Taxon")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file')
parser.add_argument('-t', '--taxon', dest = 'taxon', type = str, default = None, required = True, help = 'Name of taxon')
parser.add_argument('-l', '--flanks', dest = 'flanks', type = str, default = None, required = False, help = 'Complete assembled sequences')
args, unknown = parser.parse_known_args()


# read in sequences
records = list(SeqIO.parse(args.fasta, "fasta"))

# filter just the sequences in the correct orientation (known from the reference sequence used with exonerate)
oriented = []
	
for record in records:
	coordinates = record.description.split("\t")[1]
	start = int(coordinates.split("-")[0])
	end = int(coordinates.split("-")[1])
	if end > start:
		oriented.append(record)	

WSdict = Bio.SeqIO.to_dict(SeqIO.parse(args.flanks,"fasta"))

genomic = []
	
for sequence in oriented:
	genomicSeq = WSdict[sequence.id]
	if sequence.seq in genomicSeq.seq:
		newRecord = SeqRecord(genomicSeq.seq, id=args.taxon, description='')
		genomic.append(newRecord)
	else:
		genomicSeqRev = genomicSeq.seq.reverse_complement()	
		newRecord = SeqRecord(genomicSeqRev, id=args.taxon, description='')
		genomic.append(newRecord)
			
SeqIO.write(genomic, args.output, "fasta")		



