#!/usr/bin/env python

import argparse
import re
from sys import argv

###### Documentation

parser = argparse.ArgumentParser(description="Requires python 2.7 and the python package ete3. Testing with python 3 is ongoing. Checks monophyly of group specified by a list of names. Please note that all trees are midpoint-rooted in this script. This script was written by Lily Hughes, lilychughes@gmail.com. ")
parser.add_argument('-t', '--tree' , dest = 'tree' , type = str , default= None , required= True, help = 'Newick-formatted file')
parser.add_argument('-c', '--clade' , dest = 'clade' , type = str , default= None , required= True, help = 'Text file with taxon names you want to test for monophyly, with one taxon name per line.')
args, unknown = parser.parse_known_args()

from ete3 import Tree

#load newick tree
tree = Tree(args.tree)

#get midpoint rooting
og = tree.get_midpoint_outgroup()
tree.set_outgroup(og)

#read a list of ingroup taxa whose monophyly you want to test

clade = open(args.clade)
target_clade = []

for line in clade:
	line = line.strip("\n")
	target_clade.append(line)

clade.close()


#prune taxa that are not in each gene from target clade
pruned_target_clade = []

leaves = []

for leaf in tree:
	leaves.append(leaf.name)
	
for leaf in leaves:
	if leaf in target_clade:
		pruned_target_clade.append(leaf)

		
#check the monophyly
mono = tree.check_monophyly(values=pruned_target_clade, target_attr="name")

if True in mono:
	print(args.clade+" is monophyletic in "+args.tree)
else:
	print(args.clade+" is NOT monophyletic in "+args.tree)
	print("Intruders into "+args.clade+" are:")
	for leaf in mono[2]:
		print(leaf)	
				