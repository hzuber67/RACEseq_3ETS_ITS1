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
# GOAL: -trim read 2 from sequence containing adaptor
#		-keep only read2 longer than 19nt after trimming
#################################################################################

import re
import gzip
import sys, os
from Bio.SeqIO.QualityIO import FastqGeneralIterator


#####################################################################
### functions definition
#####################################################################

PCRprimer = '(^.+)(GATCGTCGGACTGTAGAA).+$'

def delim(b):
	global output
	Count = 0
	Count1 = 0
	Count2 = 0
	Count3 = 0
	for title, seq, qual in FastqGeneralIterator(b):
			if re.compile(PCRprimer).search(seq):
				Result=re.search(PCRprimer, seq)
				Delim2=str(Result.group(1))
				x = len(Delim2)
				seq= seq[:x]
				qual= qual[:x]
				Count += 1
				if len(seq) > 19: # 
					output.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
					Count1 += 1
			else:
				Count2 += 1
				if len(seq) > 19:
					output.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
					Count3 += 1
	print ('%i\t%i\t%i\t%i' % (Count, Count1, Count2, Count3))		
			
			
#####################################################################
### main function
#####################################################################

def main():
	print ("trimmed\ttrimmed>19\tuntrimmed\tuntrimmed>19")
	readfn = sys.argv[1]
	outfn = sys.argv[2]
	files = gzip.open(readfn,'rt')
	global output
	output = gzip.open(outfn, 'wt')
	delim(files)
	output.close()


if __name__ == "__main__" :
    main()




