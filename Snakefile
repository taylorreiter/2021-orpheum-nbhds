import pandas as pd

m = pd.read_csv("inputs/working_metadata.tsv", sep = "\t", header = 0)
LIBRARIES = m['library_name'].unique().tolist()
KSIZES = ['7', '10']
	
rule all:
    input:
        #expand("outputs/orpheum/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.faa", library = LIBRARIES, ksize = KSIZES) 
        expand("outputs/nuc_noncoding_bwa/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.flagstat", library = LIBRARIES, ksize = KSIZES), 
        expand("outputs/aa_paladin/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.aa.flagstat", library = LIBRARIES, ksize = KSIZES) 

# mkdir -p outputs/rgnv_sgc_original_results
# cd outputs/rgnv_sgc_original_results
# ln -s ../../../2020-ibd/sandbox/rgnv_sgc_original_results/*gz .

rule cat_sgc_nbhds:
    input: expand("outputs/rgnv_sgc_original_results/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz", library = LIBRARIES)
    output: "outputs/rgnv_sgc_original_results/all_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz"
    resources: mem_mb = 2000
    threads: 1
    shell:'''
    cat {input} > {output}
    '''

# some of the 605 fastq nbhds have shorter ORFs than 75 required for prediction here
# therefore, some nbhds will be better represented than others.
# e.g. qin et al. has 36 bp seqs before trimming.
rule plass_assemble_sgc_nbhds:
    input: "outputs/rgnv_sgc_original_results/all_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz"
    output: "outputs/plass/rgnv_original_sgc_nbhds_plass_assembly.faa"
    conda: "envs/plass.yml"
    benchmark: "benchmarks/plass_assemble_sgc_nbhds.txt"
    resources: mem_mb = 512000
    threads: 8
    shell:'''
    plass assemble --min-length 25 {input} {output} tmp 
    '''

rule orpheum_index_plass_assembly:
    input: "outputs/plass/rgnv_original_sgc_nbhds_plass_assembly.faa"
    output: "outputs/orpheum_index/rgnv_original_sgc_nbhds_plass_assembly_protein_ksize{ksize}.bloomfilter.nodegraph"
    conda: "envs/orpheum.yml"
    benchmark: "benchmarks/orpheum_index_plass_assembly_ksize{ksize}.txt"
    resources: mem_mb = 128000
    threads: 1
    shell:'''
    orpheum index --molecule protein --peptide-ksize {wildcards.ksize} --save-as {output} {input}
    '''

rule orpheum_translate_sgc_nbhds:        
    input: 
        ref="outputs/orpheum_index/rgnv_original_sgc_nbhds_plass_assembly_protein_ksize{ksize}.bloomfilter.nodegraph",
        fastq="outputs/rgnv_sgc_original_results/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz"
    output:
        pep="outputs/orpheum/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.faa", 
        nuc="outputs/orpheum/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_coding.fna",
        nuc_noncoding="outputs/orpheum/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.fna",
        csv="outputs/orpheum/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.coding_scores.csv",
        json="outputs/orpheum/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.summary.json"
    conda: "envs/orpheum.yml"
    benchmark: "benchmarks/orpheum_translate_{library}_plass_assembly_ksize{ksize}.txt"
    resources: mem_mb = 62000
    threads: 1
    shell:'''
    orpheum translate --peptide-ksize {wildcards.ksize}  --peptides-are-bloom-filter --noncoding-nucleotide-fasta {output.nuc_noncoding} --coding-nucleotide-fasta {output.nuc} --csv {output.csv} --json-summary {output.json} {input.ref} {input.fastq} > {output.pep}
    '''

#################################
## Evaluate
#################################

# map noncoding nucleotide sequences against the coding portions of the genome
# approximately zero should align if prediction worked well
# this reference nuc set is separate from the PLASS assembly as was derived from
# megahit assemblies + isolate genomes in RefSeq
# Original file lives here:
# 2020-ibd/sandbox/test_roary/outputs/roary_with_megahit_and_isolates/pan_genome_reference.fa

rule index_ref_nuc_set:
    input: "inputs/pan_genome_reference.fa"
    output: "inputs/pan_genome_reference.fa.bwt"
    conda: "envs/bwa.yml"
    resources: mem_mb = 2000
    threads: 1
    shell:'''
    bwa index {input}
    ''' 

rule map_nuc_noncoding_to_ref_nuc_set:        
    input: 
        ref_nuc_set= "inputs/pan_genome_reference.fa",
        ref_nuc_set_bwt= "inputs/pan_genome_reference.fa.bwt",
        nuc_noncoding="outputs/orpheum/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.fna",
    output:"outputs/nuc_noncoding_bwa/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.bam"
    conda: "envs/bwa.yml"
    resources: mem_mb = 2000
    threads: 1
    shell:'''
    bwa mem -t {threads} {input.ref_nuc_set} {input.nuc_noncoding} | samtools sort -o {output} -
    '''

rule flagstat_map_nuc_noncoding_to_ref_nuc_set:
    input: "outputs/nuc_noncoding_bwa/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.bam"
    output: "outputs/nuc_noncoding_bwa/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.flagstat"
    conda: "envs/bwa.yml"
    resources: mem_mb = 2000
    shell:'''
    samtools flagstat {input} > {output}
    '''

# map proteins against roary ref set aas, and samtools flagstat for % mapping
# % mapping should be high, close to 100%. If low, either PLASS is too
# promiscuous, or orpheum did a bad job and gave a bunch of false postive
# protein sequences
# original prot seqs live here:
# 2020-ibd/sandbox/test_roary/outputs/roary_with_megahit_and_isolates/pan_genome_reference.faa
rule paladin_index:
    input: "inputs/pan_genome_reference.faa",
    output:"inputs/pan_genome_reference.faa.bwt",
    conda: "envs/paladin.yml"
    resources: mem_mb = 8000
    threads: 1
    shell:'''
    paladin index -r3 {input}
    '''

rule paladin_align:
    input: 
        ref="inputs/pan_genome_reference.faa",
        idx="inputs/pan_genome_reference.faa",
        pep="outputs/orpheum/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.faa", 
    output: "outputs/aa_paladin/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.aa.sam"
    conda: "envs/paladin.yml"
    resources: mem_mb = 2000
    threads: 1
    shell:'''
    paladin align -t 1 -p {input.ref} {input.pep} > {output}
    '''

rule samtools_flagstat_paladin:
    input: "outputs/aa_paladin/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.aa.sam"
    output: "outputs/aa_paladin/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.aa.flagstat"
    conda: "envs/paladin.yml"
    resources: mem_mb = 2000
    threads: 1
    shell:'''
    samtools flagstat {input} > {output}
    '''
