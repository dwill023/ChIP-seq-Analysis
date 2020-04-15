#!/bin/bash

# Immediately stop on errors
 
set -uex pipefail

# the reference genome path
REF=/BT2_index_files/hg38

module load bowtie2/2.2.9


bowtie2 -p10 -x $REF -U INPUT_I_.fastq.gz -S INPUT_I.sam

bowtie2 -p10 -x $REF -U INPUT_II_.fastq.gz -S INPUT_II.sam

bowtie2 -p10 -x $REF -U Sample_I_.fastq.gz -S Sample_I.sam

bowtie2 -p10 -x $REF -U Sample_II_.fastq.gz -S Sample_II.sam


