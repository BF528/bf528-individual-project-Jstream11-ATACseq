# RNAseq

## Methods
You will start by performing QC and adapter trimming as you have done for other experiments. I would recommend you align reads using BowTie2 with the -X 2000 flag. After alignment, you will need to first remove any alignments to the mitochondrial chromosome (look at samtools view options).
After this filtering, you will want to shift your reads to account for the bias induced by the “tagmentation” process (there are various tools to accomplish this including deeptools). After this, you will want to perform a quality control analysis by looking at the fragment distribution sizes for your samples using a tool like ATACSeqQC.
You can perform peak calling using MACS3 and their recommended default parameters for ATACseq for each replicate. Similar to ChIPseq, you will need to come up with a bedtool strategy to generate a single set of “reproducible” peaks and then filter any peaks falling into blacklisted regions. After this, you will perform similar analyses to a ChIPseq.
**Data Retrieval and Preprocessing**
Data was downloaded from the European Molecular Biology Laboratory - European Bioinformatics Institute (EMBL-ENA) using the Gene Expression Omnibus (GEO) accession GSE75070. Reads were processed to remove adaptors using Trimommatic [version 0.39] with the following parameters:
```{bash}
  trimmomatic SE -threads {threads} <fastq_input> <trimmed_output> ILLUMINACLIP:<adapters_input>:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
```
The quality of the processed reads was assessed using FastQC [version 0.12.1-0] and default parameters.

**Genome Alignment and Processing**
A genome index was constructed using the 'GRCh38.primary_assembly.genome.fa.gz' file and default parameters. Processed reads were aligned to the human reference genome using Bowtie2 [version 2.5.3] and the Bowtie2 genome index using default parameters. Aligned reads were converted to BAM format with Samtools [version 1.19.2] and subsequently sorted using default parameters. Index files for the sorted BAM files were generated with default parameters to facilitate efficient retrieval of data. MultiQC [version 1.20] and Samtools flagstat were employed with default parameters to assess the quality of each sample.

**Signal Generation and Analysis**
BigWig files were generated from the sorted BAM files using deeptools [version 3.5.4] bamCoverage with default parameters. A single matrix containing information from all BigWig files was generated using deeptools multiBigwig with default parameters. Pearson correlation plots were created to assess the correlation between samples using deeptools plotCorrelation with the following parameters:
```{bash}
plotCorrelation -in <matrix_input> --corMethod spearman --skipZeros --whatToPlot heatmap --plotNumbers -o <plot_output.png>
```
**Peak Calling and Annotation**
Tag directories were constructed using HOMER [version 4.11] makeTagDir on each of the sorted BAM files with default parameters. Peaks were identified using HOMER findPeaks on each experiment replicate, using the following parameters:
```{bash}
findPeaks <path_to_RUNX1_{rep}_tagDir> -style factor -i <path_to_INP_{rep}_tagDir> -o <{rep}_peaks.txt>
```
  with outputs converted to BED format using HOMER pos2bed.pl using default parameters. Reproducible peaks between replicates were identified using bedtools [version 2.31.1] intersect using the following parameters:
```{bash}
bedtools intersect -a <rep1_peaks.bed> -b <rep2_peaks.bed> -f 0.5 -r > <reproducible_peaks.bed>
```
  and filtered against the GENCODE blacklist reference BED file using the following parameters:
```{bash}
bedtools intersect -a <reproducible_peaks.bed> -b <blacklist.bed> -v > <filtered_peaks.bed>
```
Reproducible peaks were annotated to their nearest genomic features using HOMER annotatePeaks.pl with the following parameters:
```{bash}
annotatePeaks.pl <filtered_peaks.bed> hg38 -gtf <gencode.v45.primary_assembly.annotation.gtf> > <annotated_peaks.txt>
```
Motif finding was performed on the list of filtered, reproducible peaks using HOMER findMotifsGenome with the following parameters:
```{bash}
findMotifsGenome.pl <reproducible_peaks.bed> <GRCh38.primary_assembly.genome.fa> <motif_output_directory> -size 200
```
**Functional Analysis and Visualization**
A matrix of signal values across hg38 genes was generated using deeptools computeMatrix, BED file of hg38 from UCSC Genome Browser, and previously generated BigWig files with the following parameters:
```{bash}  
computeMatrix scale-regions -S <RUNX1_{rep}.bigWig> -R <hg38_genes.bed> -b 2000 -a 2000 -o <{rep}_matrix.gz>
```
Signal across hg38 genes was visualized using deeptools plotProfile and matrix outputs with default parameters. Peaks were visualized using the Integrative Genomics Viewer (IGV), incorporating BigWig files, reproducible peaks BED file, and primary assembly GTF file. Specifically, analysis focused on MALAT1 and NEAT1 genes. Annotated peaks were integrated with DESeq2 results from the GEO accession GSE75070, applying filters padj < 0.1 and log2foldchange > 1 on genes.


## Questions to Address
Briefly remark on the quality of the sequencing reads and the alignment statistics, make sure to specifically mention the following:
  - Are there any concerning aspects of the quality control of your sequencing reads?
  - Are there any concerning aspects of the quality control related to alignment?
  - Based on all of your quality control, will you exclude any samples from further analysis?
After alignment, quickly calculate how many alignments were generated from each sample in total and how many alignments were against the mitochondrial chromosome
  - Report the total number of alignments per sample
  - Report the number of alignments against the mitochondrial genome
After performing peak calling analysis, generating a set of reproducible peaks and filtering peaks from blacklisted regions, please answer the following:
  - How many peaks are present in each of the replicates?
  - How many peaks are present in your set of reproducible peaks? What strategy did you use to determine “reproducible” peaks?
  - How many peaks remain after filtering out peaks overlapping blacklisted regions?
After performing motif analysis and gene enrichment on the peak annotations, please answer the following:
  - Briefly discuss the main results of both of these analyses
  - What can chromatin accessibility let us infer biologically?


## Deliverables
1. Produce a fragment length distribution plot for each of the samples.
2. Produce a table of how many alignments for each sample before and after filtering alignments falling on the mitochondrial chromosome.
3. Create a signal coverage plot centered on the TSS (plotProfile) for the nucleosome-free regions (NFR) and the nucleosome-bound regions (NBR).
4. You may consider fragments (<100bp) to be those from the NFR and the rest as the NBR.
5. A table containing the number of peaks called in each replicate, and the number of reproducible peaks.
6. A single BED file containing the reproducible peaks you determined from the experiment.
7. Perform motif finding on your reproducible peaks.
8. Create a single table / figure with the most interesting results.
9. Perform a gene enrichment analysis on the annotated peaks using a well-validated gene enrichment tool.
10. Create a single table / figure with the most interesting results.
11. Produce a figure that displays the proportions of regions that appear to have accessible chromatin called as a peak (Promoter, Intergenic, Intron, Exon, TTS, etc.).

