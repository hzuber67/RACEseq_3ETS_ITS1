# RACEseq_rRNA : Script for RACEseq analysis (pre-rRNA)
*Version 022024*

Copyright (c) 2016-2024 Institute of Plant Molecular Biology, Strasbourg, CNRS

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Hélène Zuber <helene.zuber@ibmp-cnrs.unistra.fr>


## Dependencies:
*  Python 3
*  Python librairies: 
	- Biopython (v1.79) http://biopython.org/wiki/Download
	- regex (v 2022.9.13) https://pypi.python.org/pypi/regex
		

## How to run the pipeline :
	
- This repository contain python scripts and bash pipelines  for running the pipeline.

- One shell script is provided for each  pre-RNA targets and used primer. Lines 37 to 40 of the shell script need to be adjusted according to your working directories:
	
```bash
Root="/Users/hzuber/Documents/Run_MiSeq/RACEseq_scripts/RACEseq_3ETS/Python3/python_scripts/" # Directory of python scripts
Datadir_R1="/Users/hzuber/Documents/Toulouse_rRNA/Papier_Maxime/RACE_analysis/ETS/NGS267/R1/" #directory with read1 fastq files
Datadir_R2="/Users/hzuber/Documents/Toulouse_rRNA/Papier_Maxime/RACE_analysis/ETS/NGS267/R2/" #directory with read2 fastq files
Outdir="/Users/hzuber/Documents/Toulouse_rRNA/Papier_Maxime/RACE_analysis/ETS/Primer1_inside/Results_NGS267/" #output directory
```

- Lines 166 and 168 of the shell script need to be adjusted according to the name of your files:
```bash
prefix=$(basename $i _L001_R1_001.fastq.gz) # prefix to be used for all output names
prefix_R2=$(basename $j _L001_R2_001.fastq.gz) #to be adjusted according to the name of your fastq file
```

- To run the pipeline, excecute the script main.sh  
- Final result files will be placed in the 'process_results' of the output directory that you have designated


	

		

		

		
		