---
editor_options: 
  markdown: 
    wrap: 72
---

## Discussion

#### Sample Quality Control

The fastQC reports show warnings or issues for the following statistics:

-   per tile sequence quality for each sample

-   per base sequence content failure

-   sequence duplication levels failure

Per tile sequence quality suggests potential issue during sample
collection with the flow cell, i.e. bubbles. The quality scores for this
statistic are moderate. Per base sequence content failure could be due
to transposase-5 sequence bias. Transposase-5 can have a preference for
cutting certain sequences, leading to an overrepresentation of adenine
and thymine. Sequence duplication levels failure suggests
overamplification during library preparation or PCR artifacts. This
persists in alignment QC results and could potentially increase
alignment ambiguity and inflate signal intensity. Overall, high
duplication levels can introduce noise and bias which should be
considered during downstream analysis

The QC data suggests potential improvements could be made to the
flowcell sample collection and library preparation protocols for future
experiments. While downstream computational filtering can help mitigate
the impact of duplicate reads to some extent, rerunning the library
preparation might be preferable for optimal data quality. However, given
the constraints of this project, I will proceed with the analysis of
these samples using appropriate filtering methods.

#### Alignment Quality

Aligned read counts can be seen in Table 1. After removing ChrM mapped
reads, there are significantly less mapped reads. Some reads might map
ambiguously to multiple locations, including both ChrM and other
chromosomes. Filtering out ChrM might remove the reads mapped to the
ChrM, but keep reads aligned to other chromosome. This would contribute
to a lower ChrM read count than expected by simply subtracting the reads
after filtering form the pre-filtered total.

The fragment size plots (Figure 1) show the large proportion of reads
less than 100bp, representing the nucleosome-free region while the
proportion around 200bp shows the nucleosome-bound region. The fragment
size distribution shows periodicity, in the inset figure, indicating
efficient capture of the nucleosomal fragments. The plots show both
samples consist of quality mapped reads, consistent with expected
ATAC-seq results.

#### Peak Calling Analysis

After peak calling, there are 30193 peaks and 26347 peaks in rep 3 and
rep 4 respectively. 210 of these peaks represent reproducible peaks that
were found in both replicates. After filtering these peaks, removing
blacklisted regions, 195 reproducible peaks were found (Table 2).

Based on the distribution of accessible chromatin region types (Figure
4), There is enrichment of ATAC-seq reads in promoter regions,
indicating open chromatin and potentially high accessibility for
transcription factors to bind and initiate gene expression. The
enrichment of ATAC-seq reads in the intron regions suggests regulatory
elements like enhancers. Enrichment of reads in the intergenic regions
represents regulatory elements between genes, like enhancers and
silencers that influence gene expression from a distance.

The height of the NFR and NBR peaks in the signal coverage plots (Figure
2), along with the relative difference in their intensities, provide
insights into the accessibility of chromatin around the TSS. The NFR
region represents accessibility of DNA to protein binding factors,
suggesting potential regulatory activity. The prominent NFR peak
indicates enrichment of reads in nucleosome-free regions, suggesting
open chromatin facilitating transcription initiation. The NBR regeion
represents areas where nucleosomes are tightly bound to DNA, hindering
protein binding and potentially reducing transcription. NFR and NBR
peaks are at similar intensities, suggesting ambiguous results.
Inflation of the NBR and/or NFR region could be a downstream effect of
the high sequence duplication levels.

#### Motif Analysis

Gene enrichment results (Figure 3) show high correlation to T cell
related pathways. The most significant GO terms are
GSE32901_NAIVE_VS_TH1_CD4_TCELL_DN and
GSE10239_NAIVE_VS_KLRG1HIGH_EFF_CD8_TCELL_DN, representing genes
down-regulated in CD4 T cells and CD8 T cells compared to naive T cells.
This suggests the analyzed cells likely belong to T cell lineage.

Motif analysis results show the top motifs for the 195 total target
sequences (reproducible blacklist filtered peaks). The top 15 motifs can
be seen in Table 3 and largely correspond to transcription factors. The
most significant motifs are from interferon regulatory factors (IRF8 and
IRF3), which have a regulatory role in the immune system as
transcriptional regulators of type1 interferon-dependent immune
responses. These results suggest the cells are poised to activate genes
involved in immune responses.

By identifying accessible regions and enriched transcription factor
motifs cell's identity, potential regulatory mechanisms, and genes that
may be poised for activation can be inferred. Chromatin accessibility
measures the potential for transcription factors to bind. Additional
analysis (i.e. RNA-seq) are needed to confirm gene expression changes.

## Results

Figure 1. Fragment length distribution plots for aligned reads of
samples ATACrep3 and ATACrep4.

![Rep3 Fragment Size
Plot](/results/rep3_frag_sizes.png)
![Rep4 Fragment Size
Plot](/results/rep4_frag_sizes.png)

Table 1. Aligned read counts of samples ATACrep3 and ATACrep4 before
mitochondrial chromosome filtering, after mitochondrial chromosome
filtering, and aligned to the mitochondrial chromosome.

| Sample | Before ChrM Filtering | After ChrM Filtering | Aligned to ChrM |
|--------|-----------------------|----------------------|-----------------|
| Rep3   | 155771633             | 38713076             | 116331965       |
| Rep4   | 114717279             | 29370934             | 84837401        |

Figure 2. Signal coverage plots of samples ATACrep3 and ATACrep4,
centered on the transcription start site (TSS) for the nucleosome-free
regions (NFR, green) and the nucleosome-bound regions (NBR, blue).

![Rep3 Signal Coverage
Plot](/results/ATACrep3_coverage_plot.png)![Rep4
Signal Coverage
Plot](/results/ATACrep4_coverage_plot.png)

Table 2. Number of peaks called in each replicate, the number of
reproducible peaks, and the number of reproducible peaks after blacklist
filtering.

| Rep3 Peaks | Rep4 Peaks | Reproducible Peaks | Reproducible Blacklist filtered peaks |
|---------------|---------------|---------------|-----------------------------|
| 30193      | 26347      | 210                | 195                                   |

Table 3. Motif finding results for reproducible, blacklist filtered
peaks.

![](/results/Known_motif_results.png)

Figure 3. Gene enrichment analysis of reproducible, blacklist filtered
peaks. Top 20 most significant GO terms are colored based on
significance (ln(P)) and plotted against the target gene ratios, the
amount of targets found versus the amount of genes in each GO term
geneset.

![](/results/Gene_enrichment_plot.png)

Figure 4. Distibution of peaks by region biologic feature type:
Promoter, Intergenic, Intron, Exon, and TTS.

![](/results/Region_distribution_plot.png)

## Methods

#### **Preprocessing**

Reads were processed to remove adaptors using Trimommatic [version 0.39]
with the following parameters:

```{bash}
  trimmomatic PE -threads 8 <input.forward.fastq> <input.reverse.fastq> <output.forward_paired.fastq> <output.forward_unpaired.fastq> <output.reverse_paired.fastq> <output.reverse_unpaired.fastq> ILLUMINACLIP:<NexteraPE-PE.fa> :2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
```

The quality of the processed reads was assessed using FastQC [version
0.12.1-0] and default parameters.

#### **Genome Alignment and Processing**

A genome index was constructed using the
'GRCh38.primary_assembly.genome.fa.gz' file and default parameters.
Processed reads were aligned to the human reference genome using Bowtie2
[version 2.5.3] and the Bowtie2 genome index using the -X 2000 flag.

```{bash}
  bowtie2 --threads 32 -X 2000 -x <GRCh38_index> -1 <forward_paired.fastq> -2 <reverse_paired.fastq> | samtools view -h -S -b -o <output.bam> 
```

Aligned reads were sorted using Samtools [version 1.19.2] sort and
default parameters and processed to remove alignments with mitochondrial
chromosome using Samtools view.

```{bash}
  samtools view -h <input.bam> | grep -v chrM | samtools sort -O bam -o <output.bam_noMT> -T .
```

Samtools flagstats and idxstats were used before and after mitochondrial
chromosome removal with default parameters. Filtered reads where indexed
with Samtools index and default parameters and shifted to account for
biases induced by tagmentation processing using deeptools [version
3.5.4] alignementSieve --ATACshift.

```{bash}
  alignmentSieve --ATACshift -b <input.bam_noMT> -o <output.bam_shifted>
```

Alignment quality control was performed with ATACSeqQC [version 4.3]
with fragSizeDist() to assess the fragment size for each sample. Aligned
reads were resorted and index using Samtools sort and index. BigWig
files were generated from the BAM files using deeptools [version 3.5.4]
bamCoverage and sorted by fragment size to separate Nucleosome Binding
Regions (NBR) from Nucleosome Free Regions (NFR).

```{bash}
  bamCoverage -b <input.bam> --minFragmentLength 100 -o <output.nbr.BigWig>
  bamCoverage -b <input.bam> --maxFragmentLength 100 -o >output.nfr.BigWig>
```

MultiQC [version 1.13.1] was performed with default parameters to
aggregate QC reports.

#### **Peak Calling**

Peaks were identified with MACS [version 3.0.0] callpeak.

```{bash}
  macs3 callpeak -f BAMPE -t <input.bam> -g hs -n <params.peaks.bed> -B -q 0.01
```

Reproducible peaks between replicates were identified using bedtools
[version 2.31.1] intersect and filtered against the GENCODE blacklist
reference BED file.

```{bash}
  bedtools intersect -a <bed1> -b <bed2> -f 0.5 -r > <rep_peaks.bed>
  bedtools intersect -a <rep_peaks.bed> -b <blacklist.bed> -v > <filtered_peaks.bed>
```

Motif finding was performed on the list of filtered, reproducible peaks
using HOMER [version 4.11] findMotifsGenome with the -size 200 flag.

```{bash}
  findMotifsGenome.pl <input.peaks.bed> <GRCh38.primary_assembly.genome.fa> <results/MotifOutput> -size 200
```

Peak annotation was performed using HOMER annotatePeaks, with -go used
for gene ontology analysis.

```{bash}
  annotatePeaks.pl <input.peaks.bed> hg38 -gtf <gencode.v45.primary_assembly.annotation.gtf> > <annotated_peaks.txt> -go <results/peaks_GO>
```

#### **Functional Analysis and Visualization**

A matrix of signal values across hg38 genes was generated using
deeptools computeMatrix, BED file of hg38 from UCSC Genome Browser, and
previously generated BigWig files.

```{bash}
  computeMatrix scale-regions -S <nbr.bigWig> <nfr.bigWig> -R <hg38_genes.bed> -b 2000 -a 2000 -o <combined_matrix.gz>
```

Signal across hg38 genes for NBR and NFR regions was visualized using
deeptools plotProfile.

```{bash}
   plotProfile -m <combined_matrix.gz> -out <coverage_plot.png> --perGroup
```

Gene enrichment analysis was performed in R using the GO results from
the geneOntology.html produced through peak annotation.
