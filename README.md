# Using orpheum on spacegraphcats nbhds to predict open reading frames in reads

### Background and motivation

**Motivating Goal**: <u>Determine typical pangenome metrics for spacegraphcats genome query neighborhoods.</u> Traditionally, pangenome metrics are calculated from gene presence/absence. This works well for isolates, where one can be relatively certain that all genes (ORFs) can be annotated in an assembly. In contrast, metagenome-assembled genome bins (MAGs) are potentially missing many genes from sequences that do not assemble or do not bin. Spacegraphcats recovers these sequences, but if a gene is not an existing reference and does not assemble, it will not be counted. 

This question ignores the fact that MAGs are composite genomes from all organisms in a community at the time of sequencing, and genes that are unobserved from low sequencing depth.

*Spacegraphcats multifasta queries* improve annotation of dominating set pieces, but still require a reference with which to perform the cDBG annotation. Using all RefSeq isolates and all metagenome assemblies of a genome query neighborhood, many dominating set pieces are annotated, but I don't have a good way to confirm whether all annotatable nodes end up annotated.

*K-mers* could replace genes as a unit of content in a pangenome, however nucleotide k-mers are brittle to evolutionary distance (e.g. third base pair wobble), so would not be a good content unit for estimating pangenome openness/closedness, and would not hint at functional units present in a metagenome.

*Protein k-mers* may be a nice compromise, where k-mers represent ~all content in the genome query neighborhood, and the protein (or dayhoff, etc.) encoding is less brittle to evolutionary distance. However, protein k-mers from FASTQ files cannot be directly compared against one another because 6-frame translation is used to generate the protein k-mers, meaning ~5/6ths of protein k-mers will be noise. If the open reading frame could be predicted directly from the reads, protein k-mers may be a great solution to the above proble.  

**Problem Statement**: <u>Determine the correct open reading frame for prokaryotic FASTQ sequencing reads to enable protein encodings for k-mers in the sequences.</u>

### Experimental design

**Data**
+ FASTQ files: spacegraphcats genome query neighborhoods (*R. gnavus*). Currently, reads are single-end, but in future iterations that could maybe be paired-end.
+ Databases:
    + PLASS assembly of all *R. gnavus* neighborhoods. Many of these protein sequences will be chimeric and/or fragmented, but the opening reading frame mostly seems to be accurate. 
    + All protein sequences from GTDB 

**Parameters**
+ k = 7 (lots of shared kmers across genomes, best for AAI sort of things)
+ k = 10 (better for distinguishing things)

**Run**
PLASS test:
+ Run PLASS assembly on cat diginorm to generate protein assembly reference
+ Generate orpheum database (`orpheum index`) with PLASS assembly
+ Run orpheum on each of query genome neighborhood (FASTQ)
    + default output stdout to coding peptides > redirect stdout
    + output parquet file of coding for a subset
    + coding nucleotide fasta 

GTDB:
+ sourmash prot signatures of 605 samples
+ run against sourmash gather prot db of gtdb
+ grab genome matches
+ download (or grab from tessa but whatever) & prokka
+ generate orpheum database...

**Evaluation**

+ time all the things!
+ snakemake benchmarks()

**Run on 10 samples**
+ align nucleotides against single sample megahit assembly
    + evaluate num unaligned
+ align nucleotides against all sample combined assembly
    + Evaluate num unalinged
    + extract unaligned
    + BLAST unaligned, make sure they're real
+ align protein sequences against references?
+ compare assemblies
    + num distinct prot kmers in plass
    + num distinct prot kmers in megahit
    + num distinct prot kmers in coding peptides


**Other**
+ back of mind -- annotate ORFs from coding peptides alone?

### data provenance

This repository tests the concepts outlined above using spacegraphcats genome query neighborhoods created by [a different project](https://github.com/dib-lab/2020-ibd/).
The genome query neighborhoods were created by querying 605 gut microbiome metagenomes using a *Rumminococus gnavus* genome (GCA_900036035.1_RGNV35913).
These query neighborhoods were generated in 2020, prior to a bug that was discovered in spacegraphcats in which a large fraction of the CAtlas was accidentally removed (https://github.com/spacegraphcats/spacegraphcats/issues/299).
Even still, they are a good "real" data set to experiment with new methods. 
Note to self: on farm, these data live at `~/github/2020-ibd/sandbox/rgnv_sgc_original_results`. 

Part of this repo also benchmarks orpheus's performance against megahit assemblies of these neighborhoods.
The code used to assemble these fastq files is recorded [here](https://github.com/dib-lab/2020-ibd/tree/master/sandbox/test_megahit_diginorm_nocat). 
Note to self: on farm, these data live at `~/github/2020-ibd/sandbox/test_megahit_diginorm_nocat`).  
