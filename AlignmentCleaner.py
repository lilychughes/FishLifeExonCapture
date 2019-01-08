#!/usr/bin/env python

import re
from sys import argv
import argparse
import os
from Bio import AlignIO
from Bio import Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser(description="Requires python 2.7 and Biopython. This script cleans alignments for single-taxon insertions, gappy sequences, more than one stop codon, and converts a single stop codon to 'NNN'. Remove of gappy sequences is controlled with the -c parameter.")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to clean.')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file.')
parser.add_argument('-c', '--coverage', dest = 'coverage', type = float, default = 0.50, required = False, help = 'Maximum proportion of gaps allowed in a sequence. Default = 0.50')
args, unknown = parser.parse_known_args()

# read in the alignment in fasta format
alignment = AlignIO.read(args.fasta,"fasta")

# get the number of columns in the alignment
length = alignment.get_alignment_length()

# get the number of taxa in the alignment
numTaxa = len(alignment)

# find the columns that are not composed of mostly gaps (where three or fewer taxa have sequence data or 'NNNN' sequences)
goodColumns = []

for x in range(0,length):
	column = alignment[:,x]
	if column.count("-") < (numTaxa-3):
		slice = alignment[:,x:x+1]
		goodColumns.append(slice)


# make an empty alignment to populate
		
goodColumnsAlignment = alignment[:,0:0]
for column in goodColumns:
	goodColumnsAlignment = goodColumnsAlignment+column

newlength = goodColumnsAlignment.get_alignment_length()

# drop taxa with more than the allowed percentage of NNNs, and more than 1 stop codon
goodTaxa = []

def RemoveStops(sequence):
    """ Replaces stop codons (standard genetic code) with 'NNN' in a nucleotide seq object and returns a new seq object"""
    TAG = 'TAG'
    TGA = 'TGA'
    TAA = 'TAA'
    clean_seq = Seq('')
    codons = [sequence[x:x+3] for x in range(0, len(sequence), 3)]
    for codon in codons:
        if TAG in codon or TGA in codon or TAA in codon:
            clean_seq = clean_seq+Seq("NNN")
        else:
            clean_seq = clean_seq+codon
    return clean_seq

for record in goodColumnsAlignment:
	seq = record.seq
	gaps = seq.count("-")
	Ns = seq.count("N")
	transeq = str(record.seq).replace("-","N")
	transeqrec = Seq(transeq)
	translated = transeqrec.translate()
	stops = translated.count("*")
	if gaps < (newlength*args.coverage) and Ns < 4:
		if stops < 1:
			goodTaxa.append(record)
		elif stops == 1:
			replace_stops = RemoveStops(seq)
			new_record = SeqRecord(replace_stops, id=record.id, description='')
			goodTaxa.append(new_record)
		elif stops > 1:
			pass	

TaxAlignment =  MultipleSeqAlignment(goodTaxa)

# write the cleaned alignment to a new file

AlignIO.write(TaxAlignment, args.output, "fasta")