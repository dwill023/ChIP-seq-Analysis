#!/bin/bash

# script to run fastQC on multiple files

# Immediately stop on errors
 
set -uex pipefail

module load parallel
module load fastqc


find *fastq.gz | parallel -j 8 fastqc {} -o /fastqc 
