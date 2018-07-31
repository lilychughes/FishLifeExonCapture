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

parser = argparse.ArgumentParser(description='Requires python 2.7 and Biopython. Filters stop codons, and sequences not divisible by 3')
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file')
args, unknown = parser.parse_known_args()


records = list(SeqIO.parse(args.fasta, "fasta"))

filtered = []

for record in records:
	if len(record.seq) % 3 == 0 and "*" not in record.seq.translate():
		filtered.append(record)
		
SeqIO.write(filtered, args.output, "fasta")		
			