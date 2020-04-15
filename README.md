# ChIP-seq-Analysis

Experimental work flow of ChIP-seq below shows some of the check points that should be done prior to generating libraries. Once the sample has passed the chromatin check and ChIP confirmation. Libraries for the control (Input, un-chip'd) and each ChIP'd sample (2 replicates minimum) can be made. Bioanalyzer results for the libraries should be free of adapter dimers (100-125 bp) and your library should span around 300bp.

![Exp_workflow](https://github.com/dwill023/ChIP-seq-Analysis/blob/master/Pics/exp_workflow.png)

Best practices for sequencing of the libraries:
- Avoid batches or distribute samples evenly over batches
- sequencing depth of 20-40 million reads for standard transcription factors and higher (40-80 million) for broad profiles like histones.

The analysis workflow below:

![analysis_workflow](https://github.com/dwill023/ChIP-seq-Analysis/blob/master/Pics/analysis_workflow.png)

## Quality Control
Usually after sequencing the core facility will automatically trim off the adapters and generate reports for each FASTQ file. If not the report can be generated for all files as shown below, using parallel to run FASTQC over multiple files. Can run [shell script](https://github.com/dwill023/ChIP-seq-Analysis/blob/master/fastQC.sh).

```Shell
module load parallel
module load fastqc

find *fastq.gz | parallel -j 8 fastqc {} -o /fastqc
```

If the quality scores of bases in the reads are good (>20) and there are no N bases in reads, then trimming can be skipped. If not follow below to trim for quality using the [Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) software.

```Shell
module load trimmomatic/0.36
trimmomatic SE -threads 4 sample.fastq.gz sample.trimmed.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30
```

After trimming, rerun FastQC to check the quality.

## Alignment to Genome
Reads are aligned against the genome using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml). An index of the genome must first be generated before alignment. Indexing the genome allows for efficient search and retrieval of matches of the query (sequence read) to the genome.

Indxing:
```Shell
module load bowtie2
bowtie2-build -f hg38.fa /BT2_index_files/hg38 --threads 4
```
Aligning:
By default Bowtie2 performs global end-to-end alignment (no clipping of the bases in the read) and is best for quality trimmed reads.
The default setting is --sensitive setting (-D 15 -R 2 -N 0 -L 22 -i S,1,1.15)
There is no parameter in Bowtie2 to keep uniquely mapped reads so they must be filtered out after alignment.
Run the alignment of the files in parallel.

```Shell
REF=/BT2_index_files/hg38
ls *fastq.gz | parallel bowtie2 -p10 -x $REF -U {} -S {.}.sam
```
## Filter uniquely mapped reads
1. Using [Samtools](http://www.htslib.org/doc/samtools.html) the .sam file will be converted to .bam files for downstream processing.
```Shell
module load samtools
ls *.sam | parallel --eta --verbose "samtools view -h -b {} > {.}.bam"
```
2. Sort the bam files by genomic coordinates using [Sambamba](http://lomereiter.github.io/sambamba/index.html).
```Shell
conda activate bioinfo
sambamba sort -t 2 -o sample.sorted.bam sample.bam
```
3. Filter out uniquely mapped reads
The filter command -F takes out multimappers `[XS]` and not unmapped and not duplicate reads.
```Shell
sambamba view -h -t 2 -f bam -F "[XS] == null and not unmapped  and not duplicate" \
sample.sorted.bam > sample.bam
```

## Peak Calling
Using MACS2 (Model-based Analysis for ChIP-Seq). Merge BAM files of replicates prior to peak calling. Merging the BAM files of technical/biological replicates can improve sensitivity of peak calling, as the number of read depth and coverage increases when the replicates are merged.

1. Merging the bam files: samtools merge – merges multiple sorted input files into a single output.
```Shell
module load samtools
samtools merge -r outfile.bam rep1.bam rep2.bam
```
2. Index the merged bam file
```Shell
ls *.bam | parallel --eta --verbose "samtools index {}"
```
3. Predict the fragment length for each merged bam file and write it down. This will be used for extending the fragment size in the peak calling function.
```Shell
conda activate bioinfo
macs2 predictd -i file.bam -g hs -m 5 20
```
4. Peak calling

--nomodel: while on MACS will bypass building the shifting model. 
--extsize: when nomodel is on, you set this parameter to define the extension of reads in 5'->3' direction to fix-sized fragment. This is used when MACS fails to build a model or when you know the size of the binding region of your protein.
For the option --extsize input the predicted fragment length for the ChIP.bam file.
```Shell
macs2 callpeak -t <ChIP.bam> -c <Control.bam> -f BAM -g hs -n ChIP -B -q 0.05 –nomodel --extsize 200
```



Remove blacklisted regions
These regions within the genome such as repetitive regions (centromeres, telomeres, satellite repeats) tend to have a very high ratio of multi-mapping to unique mapping reads and high variance of mappability. 
