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
# GOAL: - determine the 3 prime extremities of the target RNAs and identify untemplated tail
#				1) reads 2 are mapped to the reference sequence of the corresponding target RNA
#       		2) to identify reads with untemplated tails and map their 3 prime end position, the sequences 
#				   of the unmatched reads 2 are  successively trimmed from their 3 prime end
###############################################################################################################

###required modules: 
import os, sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import regex, re
from Bio import SeqIO
from Bio.Seq import Seq
import gzip

########################################################################################################################
#function for mapping

def mapping():
	#give file names
	global Counttotal 
	global Countmapped 
	global Countunmapped
	Counttotal = 0
	Countmapped = 0
	Countunmapped = 0
	for title, seq, qual in FastqGeneralIterator(readfn_file):
		Done = 'FALSE'
		Counttotal += 1
		Gene = title.split('\t')[-2]
		Target_seq = DICT[Gene]
		Target_seq_rev = str(Seq(Target_seq).reverse_complement())
		Length = len(Target_seq_rev)
		obj1 =seq[0:20]
		stringfilter = '(%s){e<=0}'% (obj1)
		# First look for perfect mapping using 30 nt of the reference sequence
		if regex.compile(stringfilter).search(Target_seq_rev):
			research = regex.compile(stringfilter, regex.BESTMATCH)
			Countmapped += 1
			for match in research.finditer(Target_seq_rev):
				output_file.write("%s\t%s\t%s\t\t0\n" % (title, Length-match.start()-11, seq[0:30]))# cordinates are corrected in order to have the position 1 for the first ETS nt
				Done = 'TRUE'
		# Second look for progressively trimming
		else:
			for i in range(0, 31): #from 0 to 30
				x= i 
				objend =seq[x:x+1]
				obj =seq[x+1:]
				if len(seq[x:]) < 20:
					Done = 'FALSE'
					break
				else:	
					string= '(%s){e<=0}(%s){i<=1,d<=1,s<=2,e<=2}'% (objend,obj) # mismatches are not accepted at the 3' end
					if regex.compile(string).search(Target_seq_rev):
						research = regex.compile(string, regex.BESTMATCH)
						Countmapped += 1
						for match in research.finditer(Target_seq_rev):
							output_file.write("%s\t%s\t%s\t%s\t%i\n" % (title, Length-match.start()-11, seq[x:x+30], seq[:x], x))
							Done = 'TRUE'
					if Done == 'TRUE':
						break
			if Done == 'FALSE':
				other_file.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
				Countunmapped +=1
				Done = 'TRUE'
	print ('%s\t%i\t%i\t%i\t%i\n' % (readfn.split('/')[-1], x, Counttotal, Countmapped, Countunmapped ))





########################################################################################################################
#main program
########################################################################################################################
### header for files

def main() :
	global readfn
	global readfn_file
	global output_file
	global DICT
	global other_file
	readfn = sys.argv[1] #file with trimmed read2
	output_fn = sys.argv[2]
	other = sys.argv[3]
	filedir = sys.argv[4]

	# open file with reference mRNA sequences
	#file with the precursor sequence
	DICT = {}
	Target_sequences = SeqIO.parse(open(filedir),'fasta')
	for fasta in Target_sequences:
		DICT[fasta.id]=str(fasta.seq)
	readfn_file=gzip.open(readfn,'rt')
	output_file=open(output_fn,'wt')
	print ('File\tPosition\ttotal\tmapped\tunmapped') # write header for counting files
	other_file = open(other,'wt')
	mapping()

if __name__ == "__main__" :
    main()

