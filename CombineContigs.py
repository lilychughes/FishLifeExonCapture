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


dummySeq = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
dummySeqRecord = SeqRecord(dummySeq)

doneSeqs = []

combinedSeq = []
for record in records:
	if len(doneSeqs) - len(records) > 1:
		combinedSeq.append(record.seq)
		combinedSeq.append(dummySeqRecord.seq)
		doneSeqs.append(seq)
	else:
		combinedSeq.append(record.seq)	

newSeq = ''.join(combinedSeq)

newSeqRecord = SeqRecord(newSeq, id=args.fasta.combined.contig.fasta, description = '')

SeqIO.write(newSeqRecord, args.fasta+"combined.contig.fasta", "fasta")

