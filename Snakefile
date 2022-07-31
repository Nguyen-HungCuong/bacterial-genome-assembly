import pandas as pd
import os, glob
from snakemake.io import expand

sample_sheet  = config['sample_sheet']
sample_dict   = sample_info.set_index('Sample',drop=False).to_dict()['name']
outdir        = config['outdir']

output_prefix = [expand(outdir + "/{Sample}_{name}", Sample=key, name=value) for key, value in sample_dict.items()]
output_suffix = [expand(outdir)]
# ----------------------------------------------------------------------------

def fastq1(wildcards):
    return sample_info['fq1'][sample_info['Sample'] == wildcards.Sample]
def fastq2(wildcards):
    return sample_info['fq2'][sample_info['Sample'] == wildcards.Sample]
def fastq_np(wildcards):
    return sample_info['fq3'][sample_info['Sample'] == wildcards.Sample]
def reference_fasta(wildcards):
    return sample_info['reference'][sample_info['Sample'] == wildcards.Sample]
def species(wildcards):
    return sample_info['species'][sample_info['Sample']== wildcards.Sample]
def gram(wildcards):
    return sample_info['gram'][sample_info['Sample']== wildcards.Sample]

# ------------------------------------------------------------------------------

rule all:


rule mkdir_sample:
    input: outdir=lambda wildcards: expand("{outdir}/{Sample}-{{name}}", Sample=sample_dict[wildcards.Sample])
    output: "{outdir}/{Sample}-{{name}}"
    shell: 
        "mkdir -p {input.outdir}"

# Check the quality of Illumina short reads
rule fastqc_initial:
    input:
        read1 = fastq1,
        read2 = fastq2
    output:
        html_read1 = "{ouput}/{Sample}-{name}_R1_fastqc.html",
        html_read2 = "{output}/{Sample}-{name}_R2_fastqc.html" 
    threads: 24
    params: outputpath = "{outdir}/{Sample}-{{dot_giai}}"
    run:
        for file in input:
            cmd = "fastqc -t " + str(threads) + " " + file + " -o " + str(params.outputpath)

# Trimming short reads
rule fastp:
    input:
        read1 = fastq1,
        read2 = fastq2
    output:
        read1_out = "{outdir}/{Sample}-{name}_read1.trimmed.fastq.gz",
        read2_out = "{outdir}/{Sample}-{name}_read2.trimmed.fastq.gz",
        html      = "{outdir}/{Sample}-{name}.trimmed.html" 
    threads: 24
    log:
        "{outdir}/{sample}.fastp.log"
    params:
        trim_front1 = 15,
        trim_tail1 = 10,
        trim_front2 = 15,
        trim_tail1 = 10
    shell:
        "fastp --in1 {input.read1} --in2 {input.read2} --thread {threads} "
        "--trim_front1 {trim_front1} --trim_tail1 {trim_tail1} "
        "--trim_front2 {trim_front2} --trim_tail2 {trim_tail2} "
        "--out1 {read1_out} --out2 {read2_out}"

# Check the quality of long read
rule nanoplot:
    input: fastq_np
    output: "{outdir}/{Sample}-{name}"
    threads: 24
    shell:
        "NanoPlot -t {threads} --fastq {input} -o {output} --maxlength 40000 --plots dot --legacy hex"

# Remove lambda phage in file fastqa
rule nanofilt_nanolyse:
    input: fastq_np               
    output:
        fastq_nanofilt = "{outdir}/{Sample}-{name}_np.trimmed.fastq.gz",
        logfile = "{outdir}/{sample}-{name}_np.trimmed.log"
    params:
        quality = 10,
        length = 500, 
        headcrop = 50
    shell:
        "gunzip -c {input} | NanoLyse | "
        "NanoFilt -q {params.quality} -l {params.length} --headcrop {params.headcrop} --logfile {output.logfile} | "
        "gzip | {output.fastq_nanofilt}"

# Using Unicycler to assembly the whole genome with 3 modes: conservative, normal and bold
rule mkdir_unicycler:
    input: outdir=lambda wildcards: expand("outdir/{Sample}-{{name}}/unicycler", Sample=sample_dict[wildcards.Sample])
    output: "{outdir}/{Sample}-{{name}}/unicycler"
    shell: 
        "mkdir -p {input.outdir}"

rule assembly_unicycler:
    input: 
        short_read1  = "{outdir}/{Sample}-{name}_read1.trimmed.fastq.gz",
        short_read2  = "{outdir}/{Sample}-{name}_read2.trimmed.fastq.gz",
        long_read    = "{outdir}/{Sample}-{name}_np.trimmed.fastq.gz"
    output: 
        conservative = directory("{outdir}/{Sample}-{name}/unicycler/{Sample}-{name}-conservative")
        normal       = directory("{outdir}/{Sample}-{name}/unicycler/{Sample}-{name}-normal")
        bold         = directory("{outdir}/{Sample}-{name}/unicycler/{Sample}-{name}-bold")
    threads: 30
    params: min_len  = 200
    log: "{outdir}/{Sample}/{Sample}-{name}.unicycler.log"
    shell:
        "Unicycler --short1 {input.short_read1} --short2 {input.short_read2} --long{input.long_read} --threads {str(threads)} --min_fasta_length {params.min_len} --mode {output.conservative}"

        cmd_normal         =  f"Unicycler --short1 {input.short_read1} --short2 {input.short_read2} --long{input.long_read} --threads {str(threads)} --min_fasta_length {params.min_len} --mode {output.normal}

        cmd_bold           = f"Unicycler --short1 {input.short_read1} --short2 {input.short_read2} --long{input.long_read} --threads {str(threads)} --min_fasta_length {params.min_len} --mode {output.bold}"

rule copy_genome_contigs:
    input:
        conservative = "{outdir}/{Sample}/unicycler/{Sample}-{name}-conservative",
        normal       = "{outdir}/{Sample}/unicycler/{Sample}-{name}-normal",
        bold         = "{outdir}/{Sample}/unicycler/{Sample}-{name}-bold"
    output:
        conservative = "{outdir}/{Sample}/{Sample}-{name}-conservative.assembly.fasta",
        normal       = "{outdir}/{Sample}/{Sample}-{name}-normal.assembly.fasta",
        bold         = "{outdir}/{Sample}/{Sample}-{name}-bold.assembly.fasta"
    shell:
        "cp {input.conservative} {output.conservative} ; cp {input.normal} {output.normal} ; cp {input.bold} {output.bold}"

# rule assess_assembly:
#     input:
#         conservative = "{outdir}/{Sample}/{Sample}-{name}-conservative.assembly.fasta",
#         normal       = "{outdir}/{Sample}/{Sample}-{name}-normal.assembly.fasta",
#         bold         = "{outdir}/{Sample}/{Sample}-{name}-bold.assembly.fasta"
#     params:
#         reference    = reference,

#     output: directory('{output}/')
#     shell:
#         "quast.py {input} -r {}"

rule mkdir_annotation:
    input: 
        outdir=lambda wildcards: expand("outdir/{Sample}-{{name}}/annotation", Sample=sample_dict[wildcards.Sample])
    output: "{outdir}/{Sample}-{{name}}/annotation"
    shell: 
        "mkdir -p {input.outdir}"

rule annotation:
    input:
        conservative = "{outdir}/{Sample}/{Sample}-{name}-conservative.assembly.fasta",
        normal       = "{outdir}/{Sample}/{Sample}-{name}-normal.assembly.fasta",
        bold         = "{outdir}/{Sample}/{Sample}-{name}-bold.assembly.fasta"
    output:

    params:
        species = species,
        gram    = gram,
        prefix  = 

    threads: 20
    shell:
        "prokka --cpus {threads} --species {params.species} --gram {params.gram} --prefix  "






 




