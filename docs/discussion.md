## Discussion
Are there any concerning aspects of the quality control of your sequencing reads?
    From fastqc:
        - potential downstream issues → per tile sequence quality for each sample → issues with flowcell
        - per base sequence content failure → 5’ bias due to transposase-5 (Tn5) sequence bias
        - sequence duplication levels failure → normal for ATAC-seq
        - adapter used -> NexteraPE-PE.fa
Are there any concerning aspects of the quality control related to alignment?

Based on all of your quality control, will you exclude any samples from further analysis? After alignment, quickly calculate how many alignments were generated from each sample in total and how many alignments were against the mitochondrial chromosome

Report the total number of alignments per sample

Report the number of alignments against the mitochondrial genome After performing peak calling analysis, generating a set of reproducible peaks and filtering peaks from blacklisted regions, please answer the following:

How many peaks are present in each of the replicates?

How many peaks are present in your set of reproducible peaks? What strategy did you use to determine “reproducible” peaks?

How many peaks remain after filtering out peaks overlapping blacklisted regions? After performing motif analysis and gene enrichment on the peak annotations, please answer the following:

Briefly discuss the main results of both of these analyses.

What can chromatin accessibility let us infer biologically?