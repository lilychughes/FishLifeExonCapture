#!/usr/bin/env python

import re
from sys import argv
import argparse
import os

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


parser = argparse.ArgumentParser(description="Requires python 2.7 and Biopython. Combines sequences with a string of 100 Ns in between.")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta to process')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file')
args, unknown = parser.parse_known_args()



records = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))


dummySeq = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'

doneSeqs = []
combinedSeq = []

for record in records.keys():
	if len(records) - len(doneSeqs) > 1:
		combinedSeq.append((str(records[record].seq))+dummySeq)
		doneSeqs.append(record)
	else:
		combinedSeq.append(str(records[record].seq))
		doneSeqs.append(record)


newSeq = Seq(''.join(combinedSeq))

newSeqRecord = SeqRecord(newSeq, id=args.fasta+".initial.combined.contig", description = '')

SeqIO.write(newSeqRecord, args.output, "fasta")

