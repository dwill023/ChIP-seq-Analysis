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

```powershell
module load parallel
module load fastqc

find *fastq.gz | parallel -j 8 fastqc {} -o /fastqc
```

If the quality scores of bases in the reads are good (>20) and there are no N bases in reads, then trimming can be skipped. If not follow below to trim for quality using the [Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) software.

```powershell
module load trimmomatic/0.36
trimmomatic SE -threads 4 sample.fastq.gz sample.trimmed.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30
```

After trimming, rerun FastQC to check the quality.

## Alignment to Genome
Reads are aligned against the genome using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml). An index of the genome must first be generated before alignment. Indexing the genome allows for efficient search and retrieval of matches of the query (sequence read) to the genome.

Indxing:
```powershell
module load bowtie2
bowtie2-build -f hg38.fa /BT2_index_files/hg38 --threads 4
```
Aligning:
By default Bowtie2 performs global end-to-end alignment (no clipping of the bases in the read) and is best for quality trimmed reads.
The default setting is --sensitive setting (-D 15 -R 2 -N 0 -L 22 -i S,1,1.15)
There is no parameter in Bowtie2 to keep uniquely mapped reads so they must be filtered out after alignment.
Run the alignment of the files in parallel.

```powershell
REF=/BT2_index_files/hg38
ls *fastq.gz | parallel bowtie2 -p10 -x $REF -U {} -S {.}.sam
```
## Filter uniquely mapped reads
1. Using [Samtools](http://www.htslib.org/doc/samtools.html) the .sam file will be converted to .bam files for downstream processing.
```powershell
module load samtools
ls *.sam | parallel --eta --verbose "samtools view -h -b {} > {.}.bam"
```
2. Sort the bam files by genomic coordinates using [Sambamba](http://lomereiter.github.io/sambamba/index.html).
```python
conda activate bioinfo
sambamba sort -t 2 -o sample.sorted.bam sample.bam
```
3. Filter out uniquely mapped reads
The filter command -F takes out multimappers `[XS]` and not unmapped and not duplicate reads.
```python
sambamba view -h -t 2 -f bam -F "[XS] == null and not unmapped  and not duplicate" \
sample.sorted.bam > sample.bam
```

## Peak Calling
Using MACS2 (Model-based Analysis for ChIP-Seq). Merge BAM files of replicates prior to peak calling. Merging the BAM files of technical/biological replicates can improve sensitivity of peak calling, as the number of read depth and coverage increases when the replicates are merged.

1. Merging the bam files: samtools merge – merges multiple sorted input files into a single output.
```powershell
module load samtools
samtools merge -r outfile.bam rep1.bam rep2.bam
```
2. Index the merged bam file
```powershell
ls *.bam | parallel --eta --verbose "samtools index {}"
```
3. Predict the fragment length for each merged bam file and write it down. This will be used for extending the fragment size in the peak calling function.
```python
conda activate bioinfo
macs2 predictd -i file.bam -g hs -m 5 20
```
4. Peak calling

--nomodel: while on MACS will bypass building the shifting model. 
--extsize: when nomodel is on, you set this parameter to define the extension of reads in 5'->3' direction to fix-sized fragment. This is used when MACS fails to build a model or when you know the size of the binding region of your protein.
For the option --extsize input the predicted fragment length for the ChIP.bam file.
```powershell
macs2 callpeak -t <ChIP.bam> -c <Control.bam> -f BAM -g hs -n ChIP -B -q 0.05 –nomodel --extsize 200
```

## Quality Assessment
Evaluate the quality of the chip-seq by generating cross-correlation plots of the chip vs the input. The plot generates two peaks, a peak of enrichment corresponding to the predominant fragment length and a peak corresponding to the read length (“phantom” peak). 
Metrics for assessing signal-to-noise ratios in a ChIP-seq experiment are given by the:
Normalized strand coefficient (NSC) and relative strand correlation (RSC). The higher the value for these metrics, the more enrichment.  The minimum possible value is 1 (no enrichment). Below is an example of what to expect from cross-correlation plots.
![cross-correlation](https://github.com/dwill023/ChIP-seq-Analysis/blob/master/Pics/quality_scores.png)

Other metrics such as the FRiP (fraction of reads in peaks) score measures the fraction of all mapped reads that fall into the called peak regions. i.e. usable reads in significantly enriched peaks divided by all usable reads. The 1% FRiP guideline works well and is useful for comparing results obtained with the same antibody across cell lines or with different antibodies against the same factor.

The complexity of the library can reflect the amount of bias created by a given experimental design, such as amplification preference which shows up as duplicate reads. By calculating PCR Bottlenecking Coefficient (PBC) and the Non-Redundant Fraction (NRF) will show us how complex the sample is. Below are the acceptable values set by [Encode](https://www.encodeproject.org/chip-seq/transcription_factor/#standards) for chip-seq experiments.
![PCB](https://github.com/dwill023/ChIP-seq-Analysis/blob/master/Pics/Image%20%5B4%5D.png)


To obtain these quality metrics use the `ChIC` R package.
```R
library(ChIC)

cc_result <- qualityScores_EM("merged bam file for sample",
                              "merged input bam file",
                              read_length = 75,
                              annotationID = "hg19",
                              mc = 8,
                              savePlotPath = "/metrics")
```

## Remove blacklisted regions
These regions within the genome such as repetitive regions (centromeres, telomeres, satellite repeats) tend to have a very high ratio of multi-mapping to unique mapping reads and high variance of mappability. [See this paper](https://www.nature.com/articles/s41598-019-45839-z).

Blacklist files can be downloaded [here](https://github.com/Boyle-Lab/Blacklist/). Then using the narrowPeak file generated by MACS2, intersect the regions in the blacklist file and remove those regions from the narrowPeak file.

```powershell
module load bedtools
# How many chip peaks overlap with black-listed regions?
bedtools intersect -a peaks.bed -b hg38-blacklist.v2.bed -wa | wc -l

# exclude those peaks
bedtools intersect -a peaks.bed -b hg38-blacklist.v2.bed -v > filtered_peaks.bed
```

## Differential Peak Binding
By performing differential binding we can assess the binding changes a protein may have within different timepoints or conditions. This information can be used to show differences in DNA regulation and infer function of the protein when it is bound to these different regions. 

1. Generate a consensus peak set for two conditions you want to compare i.e. day 0 and day 7. This is compared to the individual bam files to compute the read counts that fall under the consensus peaks. Take the **narrowpeak** files generated from macs2, store only the first three columns (chr, start, end) into a bed file.

```powershell
cut -f 1-3 /macs2/day0_macs_peaks.narrowPeak > day0.bed
cut -f 1-3 /macs2/day7_macs_peaks.narrowPeak > day7.bed
```
2. Merge the bed files using `bedops` to generate a consensus peak set.
```powershell
module load bedops/2.4.24
bedops -m day0.bed day7.bed > consensus_peaks.bed

wc -l consensus_peaks.bed # shows the number of peaks 
```
3. Generate the count matrix file (numpy array .npz file) for input into DESeq2 R package for differential expression.
```powershell
conda activate bioinfo

multiBamSummary BED-file --BED consensus_peaks.bed -bamfiles\                                                                             
d0_FOXO3B_I_.sorted.bam \
d0_FOXO3B_III_.sorted.bam \ 
d7_FOXO3B_I_.sorted.bam \
d7_FOXO3B_II_.sorted.bam \
--outFileName multiBamArray_foxo3B.npz \
--outRawCounts readCounts.tab \
--extendReads=260 \
--blackListFileName /mnt/d/Genomic_Files/Human/hg38-blacklist.v2.bed \
--labels D0FOX3B_1 D0FOX3B_2 D7FOX3B_1 D7FOX3B_2

plotCorrelation --corData multiBamArray_foxo3B.npz \
--plotFile correlation.png \
--whatToPlot heatmap --corMethod pearson --plotNumbers --removeOutliers

plotPCA -in multiBamArray_foxo3B.npz -o pca.png
```
4. Use the readCounts.tab file generated above for input into `DESeq2` for differential expression. 

## Visualization
1. Coverage plots show the peaks across the whole genome.
Use R library `ChIPseeker` to generate plot.
2. Heatmaps generated around transcription start sites (TSS) can show the density of the reads in these regions and which may show if the protein has a function at TSS. 
Use `deepTools` to generate heatmaps
3. Peaks can be visualized in IGV (Integrative Genomics Viewer).
The raw signal track must be normalized to the sequencing depth before visualization and comparison to peak called files. 
Load the peak files (filtered_peaks.bed) with the normalized bigWig (.bw) files to see the peak regions and the peaks.

## Annotation & Functional Enrichment
1. Peak annotation will show the regulatory regions (promoters, TSS, introns, UTRs) these peaks are associated with.
Use `ChIPseeker` to annotate the blacklisted filtered peak files (filtered_peaks.bed).
2. Pathway enrihchment for the annnotated peaks can provide possible regulatory gene ontologies.
Use R library `clusterProfiler`.

## Motif Discovery
A note on generating motifs: they tend to be more subjective as many peak prediction methods still produce a high number of false positives.

1. Take the filtered_peaks.bed file and isolate 500 bp centered around the summit of the peaks. The tool [MEME-ChIP](http://meme-suite.org/tools/meme-chip) requires these centered peaks.
```powershell
cat filtered_peaks.bed | awk 'BEGIN{ OFS="\t";} { midPos=$2+$10; print $1, midPos-250, midPos+250; }' > 500bp.bed
```
2. Extract the sequences from the bed file and store it in a fasta file for uploding into MEME-ChIP. First remove any 'random' or 'alt' chromosome peaks from the bed file.
```powershell
module load bedtools

bedtools getfasta -fi /genomic_files/hg38.fa -bed 500bp.bed -fo 500bp.fa
```
3. Upload the 500bp.fa file into MEME-ChIP website.
4. The motifs that are generated in a "combined.meme" file that can be used to upload to FIMO to scan upstream or promoter sequences for the listed motifs. 
5. [GOMo](http://meme-suite.org/tools/gomo) scans all promoters using nucleotide motifs you provide to determine if any motif is significantly associated with genes linked to one or more Genome Ontology (GO) terms.

## Identifying possible cis-regulatory regions
By analyzing the binding regions outside of proximal promoters, cis-regulatory regions that may be involved with gene function can be identified.Genomic Regions Enrichment of Annotations Tool [GREAT](http://bejerano.stanford.edu/great/public/html/index.php) can be used to identify cis-regulatory regions which then annotates them with nearby genes. From these genes, gene ontologies are generated to show regulatory functions.
