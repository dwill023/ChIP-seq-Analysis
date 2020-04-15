# ChIP-seq-Analysis

Experimental work flow of ChIP-seq below shows some of the check points that should be done prior to generating libraries. Once the sample has passed the chromatin check and ChIP confirmation. Libraries for the control (Input, un-chip'd) and each ChIP'd sample (2 replicates minimum) can be made. Bioanalyzer results for the libraries should be free of adapter dimers (100-125 bp) and your library should span around 300bp.

![Exp_workflow](https://github.com/dwill023/ChIP-seq-Analysis/blob/master/Pics/exp_workflow.png)

Best practices for sequencing of the libraries:
- Avoid batches or distribute samples evenly over batches
- sequencing depth of 20-40 million reads for standard transcription factors and higher (40-80 million) for broad profiles like histones.

The analysis workflow below:

![analysis_workflow](https://github.com/dwill023/ChIP-seq-Analysis/blob/master/Pics/analysis_workflow.png)

## Quality Control
Usually after sequencing the core facility will automatically trim off the adapters and generate reports for each FASTQ file. If not the report can be generated for all files as shown below, using parallel to run FASTQC over multiple files. Can run [shell script](

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
bowtie2-build -f hg38.fa hg38 --threads 4
```
Aligning:
By default Bowtie2 performs global end-to-end alignment (no clipping of the bases in the read) and is best for quality trimmed reads.
The default setting is --sensitive setting (-D 15 -R 2 -N 0 -L 22 -i S,1,1.15)
There is no parameter in Bowtie2 to keep uniquely mapped reads so they must be filtered out after alignment.

```Shell

