rule get_genome:
    input:
        fasta=lambda wildcards: (
            config["get_genome"]["fasta"]
            if config["get_genome"]["database"] == "manual"
            else []
        ),
        gff=lambda wildcards: (
            config["get_genome"]["gff"]
            if config["get_genome"]["database"] == "manual"
            else []
        ),
    output:
        fasta="results/get_genome/genome.fasta",
        gff="results/get_genome/genome.gff",
        fai="results/get_genome/genome.fasta.fai",
    params:
        database=config["get_genome"]["database"],
        assembly=config["get_genome"]["assembly"],
        gff_source_types=config["get_genome"]["gff_source_type"],
    message:
        "parsing genome GFF and FASTA files"
    log:
        path="results/get_genome/log/get_genome.log",
    wrapper:
        "https://raw.githubusercontent.com/MPUSP/mpusp-snakemake-wrappers/refs/heads/main/get_genome"


rule get_fastq:
    input:
        get_fastq,
    output:
        fastq="results/get_fastq/{sample}_{read}.fastq.gz",
    conda:
        "../envs/basic.yml"
    message:
        "obtaining fastq files"
    log:
        "results/get_fastq/{sample}_{read}.log",
    shell:
        "ln -s {input} {output.fastq};"
        "echo 'made symbolic link from {input} to {output.fastq}' > {log}"


rule fastp:
    input:
        sample=get_fastq_pairs,
    output:
        html="results/fastp/{sample}.html",
        json="results/fastp/{sample}.json",
        trimmed=expand(
            "results/fastp/{{sample}}_{read}.fastq.gz",
            read=["read1", "read2"] if is_paired_end() else ["read1"],
        ),
    log:
        "results/fastp/{sample}.log",
    message:
        "trimming and QC filtering reads using fastp"
    params:
        extra=config["processing"]["fastp"]["extra"],
    threads: 2
    resources:
        mem_mb=4096,
    wrapper:
        "v7.0.0/bio/fastp"

rule trim_galore:
    input:
        get_fastq_pairs,
    output:
        fasta_fwd="results/trim_galore/{sample}_read1.fastq.gz",
        report_fwd="results/trim_galore/{sample}_read1.fastq.gz_trimming_report.txt",
        fasta_rev="results/trim_galore/{sample}_read2.fastq.gz",
        report_rev="results/trim_galore/{sample}_read2.fastq.gz_trimming_report.txt",
    message:
        "trimming and QC filtering reads using fastp"
    threads: 2
    log:
        "results/trim_galore/{sample}.log",
    params:
        extra=config["processing"]["trim_galore"]["extra"],
    wrapper:
       "v3.14.1/bio/trim_galore/pe"

rule umi_tools_extract:
    input:
        fastq=expand(
            "results/{folder}/{{sample}}_{read}.fastq.gz",
            read=["read1", "read2"] if is_paired_end() else ["read1"],
            folder = config["processing"]["tool"]
        ),
    output:
        fastq=expand(
            "results/umi_tools/extract/{{sample}}_{read}.fastq.gz",
            read=["read1", "read2"] if is_paired_end() else ["read1"],
        ),
    log:
        "results/umi_tools/extract/{sample}.log",
    message:
        "extracting UMIs using umi_tools"
    params:
        extra=config["processing"]["umi_tools_extract"]["extra"],
        input_args=lambda wildcards, input: (
            f"-I {input.fastq[0]} -S {wildcards.sample}_read1.fastq.gz "
            f"--read2-in={input.fastq[1]} --read2-out={wildcards.sample}_read2.fastq.gz"
            if is_paired_end()
            else f"-I {input.fastq[0]} -S {wildcards.sample}_read1.fastq.gz"
        ),
    conda:
        "../envs/umi_tools.yml"
    threads: 1
    shell:
        """
        umi_tools extract {params.extra} {params.input_args} -L {log}
        mv {wildcards.sample}_*.fastq.gz results/umi_tools/extract/
        """