#!/usr/bin/env python

import re
from sys import argv
import argparse
import os

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


parser = argparse.ArgumentParser(description="Requires python 2.7 and Biopython. ")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file')
args, unknown = parser.parse_known_args()



records = list(SeqIO.parse(args.fasta, "fasta"))
records.sort(key=lambda r: -len(r))
SeqIO.write(records[0], args.output, "fasta")
