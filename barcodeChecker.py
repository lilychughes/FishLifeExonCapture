#!/usr/bin/env python

import re
from sys import argv
import argparse
import os

parser = argparse.ArgumentParser(description="Requires python 2.7, but testing with other versions")
parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment')
parser.add_argument('-b', '--blast', dest = 'blast', type = str, default = None, required = True, help = 'Blast tabluar output')
args, unknown = parser.parse_known_args()

species = args.fasta.split("_")[2:3]

BlastInput = open(args.blast)

BestHit = BlastInput.readline()

BlastInput.close()

table = BestHit.split("\t")

barcodeSP = table[1]

percID = table[2]

if species in barcodeSP:
	print "Barcode match! "+species+" with "+str(percID)+"% identity match to BOLD COI database."
else:
	print "Mismatch! Species labeled "+species+" but closest COI match is "+barcodeSP+" with "+str(percID)+"% identity."
		