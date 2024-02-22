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

#################################################################################
# Script for RACEseq analysis 
# GOAL: - extract reads 2 that correspond to target RNA
#       - extract read 2 with delimiter 'GTCAG' (delimiter is reverse complement as in read 2) 
#       - remove delimiter and randomized sequence from read 2
#		- in this version two mismatch are accepted in the delimiter sequence 'GT..G'
# 		- random sequence  information is kept in the sequence ID
#################################################################################

import re
import sys, os
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip

#####################################################################
### functions definition
#####################################################################

Delimiter = '(^[ACGT]{15,16}GT..G)'

def sort(a):
	global Counter
	global ID
	Counter= 0
	for title, seq, qual in FastqGeneralIterator(a):
		key = title.split(None,1)[0] #extract fastq sequence from read2
		value = title.split(None,1)[1]
		ID[key] = value
		Counter += 1

def delim(b):
	global Count
	global Count2
	global ID
	global output
	Count = 0
	Count2 = 0
	for title, seq, qual in FastqGeneralIterator(b):
			key= title.split(None,1)[0]
			if key in ID:#extract fastq sequence from read2
				Count += 1
				if re.compile(Delimiter).search(seq): #look for seq with delimiter and remove delimiter
					Result=re.search(Delimiter, seq)
					Delim=str(Result.group(1))
					x = len(Delim)
					random= seq[0:15]
					seq= seq[x:]
					qual= qual[x:]
					output.write("@%s\t%s\t%s\n%s\n+\n%s\n" % (key, ID[key], random, seq, qual))
					Count2 += 1			
#####################################################################
### main function
#####################################################################

def main() :
	print ("Read1 \tRead2 \tRead2 with delimiter")
	global Counter
	global Count
	global Count1
	global ID
	global output
	ID={}
	readfn = sys.argv[1]
	readfn2 = sys.argv[2]
	files = gzip.open(readfn,'rt')
	sort(files)
	files2 = gzip.open(readfn2,'rt')
	outfile = sys.argv[3]
	output = gzip.open(outfile, 'wt')
	delim(files2)
	print ("%i\t%i\t%i" % (Counter, Count, Count2))
	output.close()

		
if __name__ == "__main__" :
    main()



