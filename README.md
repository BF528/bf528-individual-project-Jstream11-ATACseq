# ATACseq

## Methods
**Preprocessing**
Reads were processed to remove adaptors using Trimommatic [version 0.39] with the following parameters:
```{bash}
  trimmomatic SE -threads {threads} <fastq_input> <trimmed_output> ILLUMINACLIP:<adapters_input>:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
```
The quality of the processed reads was assessed using FastQC [version 0.12.1-0] and default parameters.

**Genome Alignment and Processing**
A genome index was constructed using the 'GRCh38.primary_assembly.genome.fa.gz' file and default parameters. Processed reads were aligned to the human reference genome using Bowtie2 [version 2.5.3] and the Bowtie2 genome index using the -X 2000 flag. Aligned reads were processed to remove alignements with mitochondrial chromosome using Samtools [version 1.19.2] view. Filtered reads where shifted to account for biases induced by tagmentation processing using deeptools [version 3.5.4] alignementSieve --ATACshift. Alignment quality control was performed with ATACSeqQC [version 4.3] with default parameters. Aligned reads were converted to BAM format with Samtools [version 1.19.2] and subsequently sorted using default parameters. Index files for the sorted BAM files were generated with default parameters to facilitate efficient retrieval of data. BigWig files were generated from the sorted BAM files using deeptools [version 3.5.4] bamCoverage with default parameters.

**Peak Calling**
Peaks were identified with MACS [version 3.0.0] callpeak using default parameters. Reproducible peaks between replicates were identified using bedtools [version 2.31.1] intersect, using -f 0.5 -r flags, and filtered against the GENCODE blacklist reference BED file, using -v.  Motif finding was performed on the list of filtered, reproducible peaks using HOMER [version 4.11] findMotifsGenome with the -size 200 flag.

**Functional Analysis and Visualization**
A matrix of signal values across hg38 genes was generated using deeptools computeMatrix, BED file of hg38 from UCSC Genome Browser, and previously generated BigWig files. Signal across hg38 genes was visualized using deeptools plotProfile and matrix outputs with default parameters.


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

