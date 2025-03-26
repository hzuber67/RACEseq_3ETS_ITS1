#!/usr/bin/env python

# Copyright (c) 2016-2024 Institute of Plant Molecular Biology, CNRS
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# - Helene Zuber <helene.zuber@ibmp-cnrs.unistra.fr>

###############################################################################################################
# Script for RACEseq analysis 
# GOAL: - Last filtering of 3 prime tail to remove false tail due to mismatches in the 3' end extremities
#		- Sequences were reverse-complemented
###############################################################################################################

###import needed python and biopython modules
import re, regex
import os, sys
from Bio.Seq import Seq
from fractions import Fraction
from math import *


####################################################
### function to reverse complement
####################################################

def rev(a,b):
	obj =Seq(a)
	global rev_Encoded
	global rev_Modification
	rev_Encoded = obj.reverse_complement()
	obj =Seq(b)
	rev_Modification = obj.reverse_complement()

####################################################
#function to check extension length>5
####################################################

def check(x):
	SearchString1 = '^(.+)\t1:N:0:\d+\t(.+)\t([ATCGN]{15})\t(.?\d+)\t([ATCGN]+)\t([ATCG]{0,30})\t(\d{1,2})$'
	global Countnomod
	global CountmU
	global CountmA
	global CountmC
	global CountmG
	global CountU
	global CountA
	global CountC
	global CountG
	global Count2
	global Tag
	if re.compile(SearchString1).search(x):
		Result=re.search(SearchString1, x)
		Name = Result.group(1)
		Gene= Result.group(2)
		Random = Result.group(3)
		End= Result.group(4)
		Encoded = Result.group(5)
		Modification = Result.group(6)
		Modification_size = int(Result.group(7))
		j= ceil(Fraction(len(Modification)*2,3))#70% 
		MonoT = '(A){%i}' % (len(Modification))
		MonoA = '(T){%i}' % (len(Modification))
		MonoC = '(G){%i}' % (len(Modification))
		MonoG = '(C){%i}' % (len(Modification))
		
		SearchStringT = '(.*A.*){%i}' % (j)
		SearchStringA = '(.*T.*){%i}' % (j)
		SearchStringC = '(.*G.*){%i}' % (j)
		SearchStringG = '(.*C.*){%i}' % (j)
		rev(Encoded, Modification) #reverse complement sequences
		if len(Modification) == 0:
			Countnomod +=1
			Tag="No-mod"
			File.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Random, End, rev_Encoded, rev_Modification, Modification_size, Tag))
		else:
			if regex.findall(SearchStringT, Modification):
				if regex.findall(MonoT, Modification):
					Tag="onlyU"
					File.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Random, End, rev_Encoded, rev_Modification, Modification_size, Tag))
					CountmU +=1
				else:
					Tag="U-rich"
					File.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Random, End, rev_Encoded, rev_Modification, Modification_size, Tag))
					CountU +=1
			elif regex.findall(SearchStringA, Modification):
				if regex.findall(MonoA, Modification):
					Tag="onlyA"
					File.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Random, End, rev_Encoded, rev_Modification, Modification_size, Tag))
					CountmA +=1
				else:
					Tag="A-rich"
					File.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Random, End, rev_Encoded, rev_Modification, Modification_size, Tag))
					CountA +=1
			elif regex.findall(SearchStringC, Modification):
				if regex.findall(MonoC, Modification):
					Tag="onlyC"
					File.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Random, End, rev_Encoded, rev_Modification, Modification_size, Tag))
					CountmC +=1
				else:
					Tag="C-rich"
					File.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Random, End, rev_Encoded, rev_Modification, Modification_size, Tag))
					CountC +=1
			elif regex.findall(SearchStringG, Modification):
				if regex.findall(MonoG, Modification):
					Tag="onlyG"
					File.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Random, End, rev_Encoded, rev_Modification, Modification_size, Tag))
					CountmG +=1
				else:
					Tag="G-rich"
					File.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Random, End, rev_Encoded, rev_Modification, Modification_size, Tag))
					CountG +=1
			else:
				Tag="Artefact"
				File2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Random, End, rev_Encoded, rev_Modification, Modification_size, Tag))
				Count2 +=1
	

####################################################
#Main function
####################################################

def main() :
	global Countnomod
	global CountmU
	global CountmA
	global CountmC
	global CountmG
	global CountU
	global CountA
	global CountC
	global CountG
	global Count2
	global File
	global File2
	Count1= 0
	Count2= 0
	CountmU= 0
	CountmA= 0
	CountmC= 0
	CountmG= 0
	CountU= 0
	CountA= 0
	CountC= 0
	CountG= 0
	Countnomod=0
	readfn = sys.argv[1]
	outfn = sys.argv[2]
	outfn2 = sys.argv[3]
	Infile = open(readfn,'rt')
	File = open(outfn,'wt')
	File2 = open(outfn2,'wt')
	File.write("Read.ID\tName\tRandom\t3'end\tencoded sequence\tmodification\tmodification.size\tclassification\n") 
	File2.write("Read.ID\tName\tRandom\t3'end\tencoded sequence\tmodification\tmodification.size\tclassification\n")
	for Line in Infile:
		Count1+=1
		check(Line)
	print ("Read number after filtering:")
	print ("Total\tno tail\tonlyU\tonlyA\tonlyC\tonlyG\tU-rich\tA-rich\tC-rich\tG-rich\tlost")
	print ("%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i" % (Count1, Countnomod, CountmU, CountmA, CountmC, CountmG, CountU, CountA, CountC, CountG, Count2 ))


if __name__ == "__main__" :
    main()
