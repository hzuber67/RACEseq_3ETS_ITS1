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
# GOAL: extract reads that match to rRNA
#################################################################################


###import needed python and biopython modules
import regex
import os, sys
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
from collections import defaultdict

# the data directory
readfn = sys.argv[1] # deduplicated read1 fastq files
outfn = sys.argv[2] #to adjust directory where out files will be saved


# Sequence to be used for research
DICT =  {"%s" % (sys.argv[3]): "%s" % (sys.argv[4])}

	
###Main program
def main() :
	print ('All\tTarget')
	COUNT={"%s" % (sys.argv[3]): 0}	#create dictionary for count
	Infile = gzip.open(readfn,'rt')
	outfile = gzip.open(outfn,'wt')
	Countall = 0
	global Count
	Count = 0
	for title, seq, qual in FastqGeneralIterator(Infile):
		Countall += 1
		global key
		for key in DICT:
			if regex.findall('(%s){s<=1}'%(DICT[key]), seq): # 1 mm acccepted
				outfile.write("@%s\t%s\n%s\n+\n%s\n" % (title, key, seq, qual))
				COUNT[key] += 1
	outfile.close()
	for key in COUNT:
		print ('%i\t%i' % (Countall, COUNT[key]))

		
if __name__ == "__main__" :
    main()


