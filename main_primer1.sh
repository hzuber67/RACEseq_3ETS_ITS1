#!/bin/sh
#
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


#This script allows for RACE-seq analyses. This is the version for rRNA analyses. Input files fastq files generated from MiSeq 
#Require tools:
#python v3
#biopython
#regex 


##########################################
###Configuration part to be adjusted:
##########################################
###Data directories:
Root="/Users/hzuber/Documents/Run_MiSeq/RACEseq_scripts/RACEseq_3ETS/Python3/python_scripts/" # Directory of python scripts
Datadir_R1="/Users/hzuber/Documents/Toulouse_rRNA/Papier_Maxime/RACE_analysis/ETS/NGS267/R1/" #directory with read1 fastq files
Datadir_R2="/Users/hzuber/Documents/Toulouse_rRNA/Papier_Maxime/RACE_analysis/ETS/NGS267/R2/" #directory with read2 fastq files
Outdir="/Users/hzuber/Documents/Toulouse_rRNA/Papier_Maxime/RACE_analysis/ETS/Primer1_inside/Results_NGS267/" #output directory

########################################################
###Configuration part that doesn't need to be adjusted:
########################################################
### Data directories: 
Outdir_d=$Outdir"dedupl/" #output directory for deduplicated reads
Outdir_dR1=$Outdir_d"R1/" #output directory for deduplicated read 1 
Outdir_dR2=$Outdir_d"R2/" #output directory for deduplicated read 2 
Outdir_targetR1=$Outdir_d"target_R1/" #output directory for target read 1 
Outdir_targetR2=$Outdir_d"target_R2_delim/" #output directory for target read 2 with delimiter
Outdir_trimmedR2=$Outdir_d"trimmed_R2/" #output directory for trimmed read 2 
Outdir_ext=$Outdir"ext_analysis/" #output directory for extension analysis
Outdir_other=$Outdir_ext"other/" #output directory for non mapped read 2 
Outdir_results=$Outdir_ext"results/" #output directory for extension analysis results
Outdir_final=$Outdir"process_results/" #final output directory (after last filtering and reverse complementing)
Outdir_lost=$Outdir_ext"filtered_out_results/" #reads lost during last filtering 

###Program and scripts:
removeduplicate=$Root"removeduplicate_Q10.py"
targetresearch=$Root"target_research.py"
RNAdelimiter=$Root"RNA_delimiter.py"
trimming=$Root"trimming.py"
mappingloops=$Root"mapping_loops_light.py"
extension_analysis=$Root"extension_analysis.py"

##information & sequences needed for pipeline
Target="3ETS"
Sequence="TCCCGCCTCCTCCCCGTTCACC" # 21 nt primer PCR2-Fwd-inside_3’ETS+ill + 1 nt
Reference_seq="/Users/hzuber/Documents/Toulouse_rRNA/Papier_Maxime/RACE_analysis/Ref_seq_3ETS.txt" #reference sequence for 5,8S

##########################################
# Create Output directories
##########################################
if [ ! -d $Outdir ]
then
	mkdir $Outdir
	echo $Outdir": doesn't exist. Created now"
else
	echo $Outdir": exists"
fi

if [ ! -d $Outdir_d ]
then
	mkdir $Outdir_d
	echo $Outdir_d": doesn't exist. Created now"
else
	echo $Outdir_d": exists"
fi

if [ ! -d $Outdir_dR1 ]
then
	mkdir $Outdir_dR1
	echo $Outdir_dR1": doesn't exist. Created now"
else
	echo $Outdir_dR1": exists"
fi

if [ ! -d $Outdir_dR2 ]
then
	mkdir $Outdir_dR2
	echo $Outdir_dR2": doesn't exist. Created now"
else
	echo $Outdir_dR2": exists"
fi
if [ ! -d $Outdir_targetR1 ]
then
	mkdir $Outdir_targetR1
	echo $Outdir_targetR1": doesn't exist. Created now"
else
	echo $Outdir_targetR1": exists"
fi
if [ ! -d $Outdir_targetR2 ]
then
	mkdir $Outdir_targetR2
	echo $Outdir_targetR2": doesn't exist. Created now"
else
	echo $Outdir_targetR2": exists"
fi
if [ ! -d $Outdir_trimmedR2 ]
then
	mkdir $Outdir_trimmedR2
	echo $Outdir_trimmedR2": doesn't exist. Created now"
else
	echo $Outdir_trimmedR2": exists"
fi
if [ ! -d $Outdir_ext ]
then
	mkdir $Outdir_ext
	echo $Outdir_ext": doesn't exist. Created now"
else
	echo $Outdir_ext": exists"
fi
if [ ! -d $Outdir_other ]
then
	mkdir $Outdir_other
	echo $Outdir_other": doesn't exist. Created now"
else
	echo $Outdir_other": exists"
fi
if [ ! -d $Outdir_results ]
then
	mkdir $Outdir_results
	echo $Outdir_results": doesn't exist. Created now"
else
	echo $Outdir_results": exists"
fi
if [ ! -d $Outdir_final ]
then
	mkdir $Outdir_final
	echo $Outdir_final": doesn't exist. Created now"
else
	echo $Outdir_final": exists"
fi
if [ ! -d $Outdir_lost ]
then
	mkdir $Outdir_lost
	echo $Outdir_lost": doesn't exist. Created now"
else
	echo $Outdir_lost": exists"
fi

####################################################################################
##Pipeline
####################################################################################
for i in ${Datadir_R1}/*.fastq.gz;do
	prefix=$(basename $i _L001_R1_001.fastq.gz) # prefix to be used for all output names
	for j in ${Datadir_R2}/*.fastq.gz;do
		prefix_R2=$(basename $j _L001_R2_001.fastq.gz) #to be adjusted according to the name of your fastq file
		if [ $prefix_R2 == $prefix ]; then #analyse R1 file and its corresponding read 2 file
			echo "*******************FILES*******************"
			echo $(basename $i)
			echo $(basename $j)
			
			#step 1: "removeduplicate.py" deduplication of sequences with identical nucleotides 
			#in read 1 (insert) and the 1st to 15th cycle in read 2 
			#(randomized bases in 3’ adapter)
			Outname_dR1=${Outdir_dR1}/${prefix}.fastq.gz
			Outname_dR2=${Outdir_dR2}/${prefix}.fastq.gz
			echo "_______Step1_deduplication_______"
			time $removeduplicate $i $j $Outname_dR1 $Outname_dR2
			
			#step 2: "target_research.py" extraction of read1 that match to target  
			Infilename=$Outname_dR1
			Outname=${Outdir_targetR1}/${prefix}.fastq.gz
			echo "_______Step2_target research_______"
			time $targetresearch $Infilename $Outname $Target $Sequence
			
			#step 3: "RNA_delimiter.py" extracting read 2 with delimiter 
			#and removing delimiter + randomized sequence  
			Infilename_R1=$Outname
			Infilename_R2=$Outname_dR2
			Outname=${Outdir_targetR2}/${prefix}.fastq.gz
			echo "_______Step3_delimiter research_______"
			time $RNAdelimiter $Infilename_R1 $Infilename_R2 $Outname 
			
			#step 4: "trimming.py" removing 5' PCR primer sequence in read 2  
			Infilename=$Outname
			Outname=${Outdir_trimmedR2}/${prefix}.fastq.gz
			echo "_______Step4_trimming of the PCR primer sequence_______"
			time $trimming $Infilename $Outname 
			
			#step 5: "mapping_loops.py" mapping of reads 2 to the rRNA reference sequence
			Infilename=$Outname
			Outname=${Outdir_results}/${prefix}.txt
			Outname_other=${Outdir_other}/${prefix}.fastq
			echo "_______Step5_mapping of reads 2 to the rRNA reference sequence_______"
			time $mappingloops $Infilename $Outname $Outname_other $Reference_seq
			
			#step 6: "extension_analysing.py" filtering 3' tails : 3’ modifications longer than 6, 10 and 15 nt 
			#were considered only if they contained at least one stretch of AAA or TTT, 
			#two stretches of AAA or TTT, or three stretches of AAA or TTT, respectively
			Infilename=$Outdir_results/${prefix}.txt
			Outname_final=${Outdir_final}/${prefix}.txt
			Outname_lost=${Outdir_lost}/${prefix}.txt
			echo "_______Step6_filtering 3' tails and reverse complementing sequences_______"
			time $extension_analysis $Infilename $Outname_final $Outname_lost 
		fi
	done
done
