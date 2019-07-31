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

parser = argparse.ArgumentParser(description="Requires python 2.7 and Biopython. This script cleans alignments for single-taxon insertions, short sequences, gappy edges, sequences with more than one stop codon, and converts a single stop codon to 'NNN'. Removal of gappy sequences is controlled with the -c parameter.")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to clean.')
parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file.')
parser.add_argument('-c', '--coverage', dest = 'coverage', type = float, default = 0.50, required = True, help = 'Maximum proportion of gaps allowed in a sequence. Default = 0.50')
parser.add_argument('-m', '--mitochondrial', dest = 'mitochondrial', type = str, default = 'False', required = False, help = 'If set to True, replaces single stop codons with NNN using Vertebrate Mitochondrial genetic code instead of the standard genetic code. Default = False')
parser.add_argument('-t', '--trim', dest = 'trim', type = float, default = 0.6, required = True, help = 'The maximum proportion of gaps to allow at the edges of the alignment. Edges with more than this proportion of gaps will be trimmed. Setting this value to 1 disables trimming. Default = 0.6')
args, unknown = parser.parse_known_args()

# read in the alignment in fasta format
rawAlignment = AlignIO.read(args.fasta,"fasta")

# find the columns that are not composed of mostly gaps (where two or fewer taxa have sequence data or 'NNNN' sequences)
goodColumns = []

for x in range(0,rawAlignment.get_alignment_length()):
    column = rawAlignment[:,x]
    if column.count("-") < (len(rawAlignment)-2):
        slice = rawAlignment[:,x:x+1]
        goodColumns.append(slice)
# make an empty alignment to populate

goodColumnsAlignment = rawAlignment[:,0:0]
for column in goodColumns:
    goodColumnsAlignment = goodColumnsAlignment+column

# Remove gappy edges using a user-specified gap threshold. This trimmer is codon-aware. Use the AlignmentCleanerSimple.py script if you are not trying to maintain reading frames.

def TrimEdges(alignment):
    """Trims gappy edges of alignments with more than the default 60% gaps given an alignment object"""
    codons = [alignment[:,x:x+3] for x in range(0, alignment.get_alignment_length(), 3)]
    percentages = []
    goodCodonsIndices = []
    for codon in codons:
        gapPerc = float(codon[:,0].count("-")+codon[:,1].count("-")+codon[:,2].count("-"))/(len(codon)*3)
        percentages.append(gapPerc)
    for i in range(0,len(percentages)-1):
        if percentages[i] < args.trim:
           goodCodonsIndices.append(i)
    goodCodons = codons[goodCodonsIndices[0]:goodCodonsIndices[len(goodCodonsIndices)-1]]
    cleanedAlignment =  alignment[:,0:0]
    for codon in goodCodons:
        cleanedAlignment = cleanedAlignment+codon
    return(cleanedAlignment)

trimmedAlignment = TrimEdges(goodColumnsAlignment)


# drop taxa with more than 1 stop codon, and replace single stop codons with 'NNN'
goodTaxa = []

def RemoveStops(sequence):
    """ Replaces stop codons (standard or vertebrate mitochondrial genetic code) with 'NNN' in a nucleotide seq object and returns a new seq object"""
    clean_seq = Seq('')
    codons = [sequence[x:x+3] for x in range(0, len(sequence), 3)]
    if 'True' in args.mitochondrial:
        for codon in codons:
            if '---' in codon :
                clean_seq = clean_seq+codon
            elif codon.translate(table="Vertebrate Mitochondrial") == '*':
                clean_seq = clean_seq+Seq("NNN")
            else:
                clean_seq = clean_seq+codon
    else:
        for codon in codons:
            if '---' in codon:
                clean_seq = clean_seq+codon
            elif codon.translate() == '*':
                clean_seq = clean_seq+Seq("NNN")
            else:
                clean_seq = clean_seq+codon
    return clean_seq

def StopCounter(sequence):
    """Counts the number of stop codons in a Seq object for standard or vertebrate mitochondrial genetic code"""
    count = 0
    codons = [sequence[x:x+3] for x in range(0, len(sequence), 3)]
    if 'True' in args.mitochondrial:
        for codon in codons:
            if '---' in codon :
                pass
            elif codon.translate(table="Vertebrate Mitochondrial") == '*':
                count = count+1
            else:
                pass
    else:
        for codon in codons:
            if '---' in codon :
                pass
            elif codon.translate() == '*':
                count = count+1
            else:
                pass
    return(count)


for record in trimmedAlignment:
    seq = record.seq
    gaps = seq.count("-")
    stops = StopCounter(seq)
    if gaps < (trimmedAlignment.get_alignment_length()*args.coverage):
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