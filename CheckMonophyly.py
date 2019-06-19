#!/usr/bin/env python

import argparse
import re
from sys import argv

###### Documentation

parser = argparse.ArgumentParser(description="Requires python 2.7 and ete3. Checks monophyly of group specified by name. This script was written by Lily Hughes, lilychughes@gmail.com. ")
parser.add_argument('-t', '--tree' , dest = 'tree' , type = str , default= None , required= True, help = 'Newick-formatted file')
parser.add_argument('-c', '--clade' , dest = 'clade' , type = str , default= None , required= True, help = 'Name assigned to clade that you want to test for monophyly')
args, unknown = parser.parse_known_args()

from ete3 import Tree

#load newick tree
tree = Tree(args.tree)

#get midpoint rooting
og = tree.get_midpoint_outgroup()
tree.set_outgroup(og)

#make a list of ingroup taxa whose monophyly you want to test
target_clade = []
for leaf in tree:
	if args.clade in leaf.name:
		target_clade.append(leaf.name)
		
#check the monophyly
mono = tree.check_monophyly(values=target_clade, target_attr="name")

if True in mono:
	print(args.clade+" is monophyletic in "+args.tree)
else:
	print(args.clade+" is NOT monophyletic in "+args.tree)
	print("Intruders into "+args.clade+" are:")
	for leaf in mono[2]:
		print(leaf)	