import pandas as pd

TMPDIR = "/scratch/tereiter"

# define samples (libraries) and orpheum databases
m = pd.read_csv("inputs/working_metadata.tsv", sep = "\t", header = 0)
LIBRARIES = m['library_name'].unique().tolist()
#ORPHEUM_DB = ['plass_assembly', "roary_with_megahit_and_isolates"]
ORPHEUM_DB = ["roary_with_megahit_and_isolates", "ruminococcusB", "f__Lachnospiraceae", "p__Firmicutes_A", "gtdb-rs202"]
#ALPHABET = ['protein', 'dayhoff']
#KSIZES = ['7', '10', '15', '17']

# set constrained k sizes
#dayhoff_ksizes = [11, 13, 15, 17]
protein_ksizes = [10]
# Snakemake will use the ALPHA_KSIZE wildcard from rule all to generate output file names
# Then, snakemake will back propagate the strings from the final file names to solve for
# the wildcards "alphabet" and "ksize" throughout the rest of the workflow. 
# The underscore for the chrs in the list ALPHA_KSIZE separates the alphabet string from 
# the ksize string, allowing snakemake to solve {alphabet}_{ksize} wildcard strings. 
# Therefore, the chrs in the ALPHA_KSIZE list also set the alphabet names as "dayhoff" and "protein".
ALPHA_KSIZE = expand('protein_ksize{k}', k=protein_ksizes)
#ALPHA_KSIZE += expand('dayhoff_ksize{k}', k=dayhoff_ksizes)

rule all:
    input:
        "outputs/rgnv_sgc_original_paladin/multiqc_report.html",
        expand("outputs/aa_paladin/{orpheum_db}/{alpha_ksize}/multiqc_report.html", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE),
        expand("outputs/nuc_noncoding_bwa/{orpheum_db}/{alpha_ksize}/multiqc_report.html", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE),
        expand("outputs/nuc_coding_bwa/{orpheum_db}/{alpha_ksize}/multiqc_report.html", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE),
        expand("outputs/nuc_noncoding_bwa/{orpheum_db}/{alpha_ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.stat", alpha_ksize=ALPHA_KSIZE, orpheum_db = ORPHEUM_DB, library = LIBRARIES),
        expand("outputs/nuc_coding_bwa/{orpheum_db}/{alpha_ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_coding.stat", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE, library = LIBRARIES),


# mkdir -p outputs/rgnv_sgc_original_results
# cd outputs/rgnv_sgc_original_results
# ln -s ../../../2020-ibd/sandbox/rgnv_sgc_original_results/*gz .

rule cat_sgc_nbhds:
    input: expand("outputs/rgnv_sgc_original_results/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz", library = LIBRARIES)
    output: "outputs/rgnv_sgc_original_results/all_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz"
    resources: 
        mem_mb = 2000,
        tmpdir = TMPDIR
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
    resources: 
        mem_mb = 512000,
        tmpdir = TMPDIR
    threads: 8
    shell:'''
    plass assemble --min-length 25 {input} {output} tmp 
    '''

rule orpheum_index_plass_assembly:
    input: "outputs/plass/rgnv_original_sgc_nbhds_plass_assembly.faa"
    output: "outputs/orpheum_index/plass_assembly_{alphabet}_ksize{ksize}.bloomfilter.nodegraph"
    conda: "envs/orpheum.yml"
    benchmark: "benchmarks/orpheum_index_plass_assembly_{alphabet}_ksize{ksize}.txt"
    resources: 
        mem_mb = 128000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    orpheum index --molecule {wildcards.alphabet} --peptide-ksize {wildcards.ksize} --save-as {output} {input}
    '''

rule remove_stop_codons:
    input: "inputs/pan_genome_reference.faa"
    output: "inputs/pan_genome_reference.faa.nostop.fa"
    conda: "envs/orpheum.yml"
    threads: 1
    resources: 
        mem_mb = 2000,
        tmpdir = TMPDIR
    shell:'''
    python scripts/remove-stop-plass.py {input}
    '''

rule orpheum_index_roary_with_megahit_and_isolates:
    input: "inputs/pan_genome_reference.faa.nostop.fa"
    output: "outputs/orpheum_index/roary_with_megahit_and_isolates_{alphabet}_ksize{ksize}.bloomfilter.nodegraph"
    conda: "envs/orpheum.yml"
    benchmark: "benchmarks/orpheum_index_pan_genome_reference_{alphabet}_ksize{ksize}.txt"
    resources: 
        mem_mb = 32000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    orpheum index --alphabet {wildcards.alphabet} --peptide-ksize {wildcards.ksize} --save-as {output} {input}
    '''

# GTDB nodegraphs cp'd from tessa's dir, see readme in index folder
rule orpheum_translate_sgc_nbhds:        
    input: 
        ref="outputs/orpheum_index/{orpheum_db}_{alphabet}_ksize{ksize}.bloomfilter.nodegraph",
        fastq="outputs/rgnv_sgc_original_results/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz"
    output:
        pep="outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.faa", 
        nuc="outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_coding.fna",
        nuc_noncoding="outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.fna",
        csv="outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.coding_scores.csv",
        json="outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.summary.json"
    conda: "envs/orpheum.yml"
    benchmark: "benchmarks/orpheum_translate_{library}_{orpheum_db}_{alphabet}_ksize{ksize}.txt"
    #resources: mem_mb = lambda wildcards, attempt: attempt*300000
    resources: 
        mem_mb = 500000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    orpheum translate --alphabet {wildcards.alphabet} --peptide-ksize {wildcards.ksize}  --peptides-are-bloom-filter --noncoding-nucleotide-fasta {output.nuc_noncoding} --coding-nucleotide-fasta {output.nuc} --csv {output.csv} --json-summary {output.json} {input.ref} {input.fastq} > {output.pep}
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
    resources: 
        mem_mb = 2000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    bwa index {input}
    ''' 

rule cut_dedup_nuc_noncoding_read_names:
    input: "outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.fna",
    output:  "outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.cut.dedup.fna.gz",
    resources:
        mem_mb = 8000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    sed '/^>/ s/__.*//g' {input} | awk '/^>/{{f=!d[$1];d[$1]=1}}f' | gzip > {output}
    '''

rule grab_cut_dedup_coding_read_names:
    input: "outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.faa", 
    output: "outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads_aa_names.cut.dedup.txt",  
    resources: 
        mem_mb = 8000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    grep ">" {input} | sed '/^>/ s/__.*//g' | awk '/^>/{{f=!d[$1];d[$1]=1}}f' > {output}
    '''

rule isolate_noncoding_only_reads:
    input: 
        pep = "outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads_aa_names.cut.dedup.txt",  
        noncoding = "outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.cut.dedup.fna.gz" 
    output: "outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.cut.dedup.only.fna.gz",
    resources: 
        mem_mb = 8000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    cat {input.noncoding} | paste - - | grep -v -F -f {input.pep} | tr "\t" "\n" | gzip > {output}
    """

rule map_nuc_noncoding_to_ref_nuc_set:        
    input: 
        ref_nuc_set= "inputs/pan_genome_reference.fa",
        ref_nuc_set_bwt= "inputs/pan_genome_reference.fa.bwt",
        nuc_noncoding="outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.cut.dedup.only.fna.gz",
    output: temp("outputs/nuc_noncoding_bwa/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.bam")
    conda: "envs/bwa.yml"
    resources:
        mem_mb = 2000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    bwa mem -t {threads} {input.ref_nuc_set} {input.nuc_noncoding} | samtools sort -o {output} -
    '''

rule flagstat_map_nuc_noncoding_to_ref_nuc_set:
    input: "outputs/nuc_noncoding_bwa/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.bam"
    output: "outputs/nuc_noncoding_bwa/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.flagstat"
    conda: "envs/bwa.yml"
    resources:
        mem_mb = 2000,
        tmpdir = TMPDIR
    shell:'''
    samtools flagstat {input} > {output}
    '''

rule multiqc_flagstat_map_nuc_noncoding_to_ref_nuc_set:
    input: expand("outputs/nuc_noncoding_bwa/{{orpheum_db}}/{{alphabet}}_ksize{{ksize}}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.flagstat", library = LIBRARIES), 
    output: "outputs/nuc_noncoding_bwa/{orpheum_db}/{alphabet}_ksize{ksize}/multiqc_report.html"
    params: 
        iodir = lambda wildcards: "outputs/nuc_noncoding_bwa/" + wildcards.orpheum_db + "/" + wildcards.alphabet + "_ksize" + wildcards.ksize,
    conda: "envs/multiqc.yml"
    resources: 
        mem_mb = 8000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    multiqc {params.iodir} -o {params.iodir} 
    '''

rule stat_map_nuc_noncoding_to_ref_nuc_set:
    input: "outputs/nuc_noncoding_bwa/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.bam"
    output: "outputs/nuc_noncoding_bwa/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_noncoding.stat"
    conda: "envs/bwa.yml"
    resources:
        mem_mb = 2000,
        tmpdir = TMPDIR
    shell:'''
    samtools view -h {input} | samtools stats | grep '^SN' | cut -f 2- > {output}
    '''

# map proteins against roary ref set aas, and samtools flagstat for % mapping
# % mapping should be high, close to 100%. If low, either PLASS is too
# promiscuous, or orpheum did a bad job and gave a bunch of false postive
# protein sequences
# original prot seqs live here:
# 2020-ibd/sandbox/test_roary/outputs/roary_with_megahit_and_isolates/pan_genome_reference.faa
rule paladin_index:
    input: "inputs/pan_genome_reference.faa",
    output:"inputs/pan_genome_reference.faa.pro",
    conda: "envs/paladin.yml"
    resources: 
        mem_mb = 8000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    paladin index -r3 {input}
    '''

rule paladin_align:
    input: 
        ref="inputs/pan_genome_reference.faa",
        idx="inputs/pan_genome_reference.faa.pro",
        pep="outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.faa", 
    output: temp("outputs/aa_paladin/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.aa.sam")
    conda: "envs/paladin.yml"
    resources: 
        mem_mb = 2000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    paladin align -t 1 -p {input.ref} {input.pep} > {output}
    '''

rule samtools_flagstat_paladin:
    input: "outputs/aa_paladin/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.aa.sam"
    output: "outputs/aa_paladin/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.aa.flagstat"
    conda: "envs/paladin.yml"
    resources: 
        mem_mb = 2000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    samtools flagstat {input} > {output}
    '''

rule multiqc_samtools_flagstat_paladin:
    input: expand("outputs/aa_paladin/{{orpheum_db}}/{{alphabet}}_ksize{{ksize}}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.aa.flagstat", library = LIBRARIES), 
    output: "outputs/aa_paladin/{orpheum_db}/{alphabet}_ksize{ksize}/multiqc_report.html"
    resources: mem_mb = 8000
    threads: 1
    params: 
        iodir = lambda wildcards: "outputs/aa_paladin/" + wildcards.orpheum_db + "/" + wildcards.alphabet + "_ksize" + wildcards.ksize,
    conda: "envs/multiqc.yml"
    resources: 
        mem_mb = 8000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    multiqc {params.iodir} -o {params.iodir} 
    '''

# map (deduplicated) nucleotide coding sequences back to nucleotide reference
# This number allows us to see how many more reads we discover as protein coding 
# using orpheum than we would discover from mapping against reference. 
# Ideally, it will be (much) lower than average AA map against AA reference 
# pangenome be predicted by orpheum. 


rule cut_dedup_nuc_coding_read_names:
    input: "outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_coding.fna",
    output:  "outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_coding.cut.dedup.fna.gz",
    resources: 
        mem_mb = 8000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    sed '/^>/ s/__.*//g' {input} | awk '/^>/{{f=!d[$1];d[$1]=1}}f' | gzip > {output}
    '''

rule map_nuc_coding_to_ref_nuc_set:        
    input: 
        ref_nuc_set= "inputs/pan_genome_reference.fa",
        ref_nuc_set_bwt= "inputs/pan_genome_reference.fa.bwt",
        nuc_noncoding="outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_coding.cut.dedup.fna.gz",
    output:temp("outputs/nuc_coding_bwa/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_coding.bam")
    conda: "envs/bwa.yml"
    resources: 
        mem_mb = 2000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    bwa mem -t {threads} {input.ref_nuc_set} {input.nuc_noncoding} | samtools sort -o {output} -
    '''

rule flagstat_map_nuc_coding_to_ref_nuc_set:
    input: "outputs/nuc_coding_bwa/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_coding.bam"
    output: "outputs/nuc_coding_bwa/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_coding.flagstat"
    conda: "envs/bwa.yml"
    resources: 
        mem_mb = 2000,
        tmpdir = TMPDIR
    shell:'''
    samtools flagstat {input} > {output}
    '''

rule multiqc_flagstat_map_nuc_coding_to_ref_nuc_set:
    input: expand("outputs/nuc_coding_bwa/{{orpheum_db}}/{{alphabet}}_ksize{{ksize}}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_coding.flagstat", library = LIBRARIES), 
    output: "outputs/nuc_coding_bwa/{orpheum_db}/{alphabet}_ksize{ksize}/multiqc_report.html"
    params: 
        iodir = lambda wildcards: "outputs/nuc_coding_bwa/" + wildcards.orpheum_db + "/" + wildcards.alphabet + "_ksize" + wildcards.ksize,
    conda: "envs/multiqc.yml"
    resources: 
        mem_mb = 8000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    multiqc {params.iodir} -o {params.iodir} 
    '''

rule stat_map_nuc_coding_to_ref_nuc_set:
    input: "outputs/nuc_coding_bwa/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_coding.bam"
    output: "outputs/nuc_coding_bwa/{orpheum_db}/{alphabet}_ksize{ksize}/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.nuc_coding.stat"
    conda: "envs/bwa.yml"
    resources:
        mem_mb = 2000,
        tmpdir = TMPDIR
    shell:'''
    samtools view -h {input} | samtools stats | grep '^SN' | cut -f 2- > {output}
    '''

#########################################
## Controls
#########################################

# Average all reads map against AA reference pangenome:
#    This number estimates the lower limit of reads that should be protein coding; 
#    at least this many reads should be predicted by orpheum

rule paladin_og_fastq_seqs:
    input: 
        ref="inputs/pan_genome_reference.faa",
        idx="inputs/pan_genome_reference.faa.pro",
        fastq="outputs/rgnv_sgc_original_results/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz"
    output: temp("outputs/rgnv_sgc_original_paladin/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.sam")
    conda: "envs/paladin.yml"
    resources: 
        mem_mb = 2000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    paladin align -t 1 {input.ref} {input.fastq} > {output}
    '''

rule flagstat_paladin_og_fastq_seqs:
    input: "outputs/rgnv_sgc_original_paladin/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.sam"
    output: "outputs/rgnv_sgc_original_paladin/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.flagstat"
    conda: "envs/paladin.yml"
    resources:
        mem_mb = 2000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    samtools flagstat {input} > {output}
    '''

rule multiqc_flagstat_paladin_og_fastq_seqs:
    input: expand("outputs/rgnv_sgc_original_paladin/{library}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.flagstat", library = LIBRARIES)
    output: "outputs/rgnv_sgc_original_paladin/multiqc_report.html"
    params: 
        indir = "outputs/rgnv_sgc_original_paladin",
        outdir = "outputs/rgnv_sgc_original_paladin"
    conda: "envs/multiqc.yml"
    resources: 
        mem_mb = 8000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    multiqc {params.indir} -o {params.outdir} 
    '''

