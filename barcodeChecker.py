#!/usr/bin/env python

import re
from sys import argv
import argparse
import os

parser = argparse.ArgumentParser(description="Requires python 2.7, but testing with other versions")
parser.add_argument('-b', '--blast', dest = 'blast', type = str, default = None, required = True, help = 'Blast tabluar output')
args, unknown = parser.parse_known_args()

species = args.blast.split(".")[0]

BlastInput = open(args.blast)

hits = []
idents = []

for line in BlastInput:
	barcodeSeq = line.split("\t")[1]
	percID = line.split("\t")[2]
	hits.append(barcodeSeq)
	idents.append(percID)

BlastInput.close()


dict = dict(zip(hits,idents))

# species ID matches and percent identity is >97%
bestMatch = []

# for matches that have high identity (>97%), but list a different species
goodMatches = []

# for matches >90%, 
mediumMatches = []

# for matches <90%
badMatches = []

for item in dict.keys():
	barcodeSP = item.split("|")[1]
	if barcodeSP in species and float(dict[item]) >= 97.0:
		bestMatch.append(item)
	elif float(dict[item]) >= 97.0:
	    goodMatches.append(item)	
	elif float(dict[item]) >= 90.0 and float(dict[item]) < 97.0:
		mediumMatches.append(item)
	elif float(dict[item]) < 90.0:
		badMatches.append(item)		

if len(bestMatch) > 0:
	key = bestMatch[0]
	perc = dict[key]
	print "Barcode match! "+species+" matches "+key+" with "+perc+"% identity."
elif len(goodMatches) > 0:
	key = goodMatches[0]
	perc = dict[key]
	print "No exact match. Best hit for "+species+" matches "+key+" with "+perc+"% identity."
elif len(mediumMatches) > 0:
	key = mediumMatches[0]
	perc = dict[key]
	print "No exact match. Best hit for "+species+" matches "+key+" with "+perc+"% identity."
elif len(badMatches) > 0:
	key = badMatches[0]
	perc = dict[key]
	print "No exact match. Best hit for "+species+" matches "+key+" with "+perc+"% identity."

		
		