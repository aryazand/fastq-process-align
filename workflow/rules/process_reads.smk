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


rule umi_tools_extract:
    input:
        fastq=expand(
            "results/get_fastq/{{sample}}_{read}.fastq.gz",
            read=["read1", "read2"] if is_paired_end() else ["read1"],
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
        input_args=lambda wildcards, input, output: (
            f"-I {input.fastq[0]} -S {output.fastq[0]} "
            f"--read2-in={input.fastq[1]} --read2-out={output.fastq[1]}"
            if is_paired_end()
            else f"-I {input.fastq[0]} -S {output.fastq[0]}"
        ),
    container:
        "docker://quay.io/biocontainers/umi_tools:1.1.6--py310h1fe012e_0"
    threads: 1
    shell:
        """
        mkdir -p $(dirname {output.fastq[0]})
        umi_tools extract {params.extra} {params.input_args} -L {log}
        """


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
