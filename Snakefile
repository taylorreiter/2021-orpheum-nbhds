import pandas as pd

m = pd.read_csv("inputs/working_metadata.tsv", sep = "\t", header = 0)
LIBRARIES = m['library_name'].unique().tolist()
KSIZES = ['7', '10']
	
rule all:
    input:
        expand("outputs/orpheum/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.faa", library = LIBRARIES, ksize = KSIZES) 


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
        ref="outputs/orpheum_index/rgnv_original_sgc_nbhds_plass_assembly_protein_ksize7.bloomfilter.nodegraph",
        fastq="outputs/rgnv_sgc_original_results/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz"
    output:
        pep="outputs/orpheum/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.faa", 
        nuc="outputs/orpheum/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_coding.fna",
        nuc_noncoding="outputs/ksize{ksize}/orpheum/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.fna",
        csv="outputs/orpheum/ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.coding_scores.csv",
        json="outputs/orpheum/ksize{ksie}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.summary.json"
    conda: "envs/orpheum.yml"
    benchmark: "benchmarks/orpheum_translate_{library}_plass_assembly_ksize{ksize}.txt"
    resources: mem_mb = 64000
    threads: 1
    shell:'''
    orpheum translate --peptides-are-bloom-filter --noncoding-nucleotide-fasta {output.nuc_noncoding} --coding-nucleotide-fasta {output.nuc} --csv {output.csv} --json-summary {output.json} {input.ref} {input.fastq} > {output.pep}
    '''
