# import basic packages
import pandas as pd
from snakemake.utils import validate
from pathlib import Path
import re

# read sample sheet
samples = (
    pd.read_csv(config["samplesheet"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)


# validate sample sheet and config file
validate(samples, schema="../../config/schemas/samples.schema.yml")
validate(config, schema="../../config/schemas/config.schema.yml")


# determine input type
def is_paired_end():
    if samples["read2"].isna().all():
        return False
    elif samples["read2"].notna().all():
        return True
    else:
        raise ValueError(
            f"Some samples seem to have a read2 fastq file, while others have only a "
            + "read1 fastq file. \nYou may not mix single-end and paired-end samples."
        )


# get fastq files
def get_fastq(wildcards):
    file = Path(samples.loc[wildcards["sample"]][wildcards["read"]])
    if file.is_absolute():
        return file
    else:
        input_dir = Path.absolute(Path.cwd())
        return input_dir / file


# get pairs of fastq files for fastp
def get_fastq_pairs(wildcards):
    return expand(
        "results/{folder}/{sample}_{read}.fastq.gz",
        folder=(
            "umi_tools/extract"
            if config["processing"]["umi_tools_extract"]["enabled"]
            else "get_fastq"
        ),
        sample=wildcards.sample,
        read=["read1", "read2"] if is_paired_end() else ["read1"],
    )


# get processed fastq files (after fastp or umi_tools)
def get_processed_fastq(wildcards, regex=None):

    if config["processing"]["umi_tools_extract"]["enabled"]:
        processed_fastq = expand(
            "results/umi_tools/extract/{{sample}}_{read}.fastq.gz",
            read=["read1", "read2"] if is_paired_end() else ["read1"],
        )
    else:
        processed_fastq = expand(
            "results/{tool}/{{sample}}_{read}.fastq.gz",
            read=["read1", "read2"] if is_paired_end() else ["read1"],
            tool=config["processing"]["tool"],
        )
    if regex is None:
        return processed_fastq
    else:
        return [s for s in processed_fastq if re.search(regex, s)]


# determine processing tool output directory
def get_processing_dir():
    return f"results/{config['processing']['tool']}"


# determine version of genome to get
def get_genome_for_mapping(wildcards):
    if len(config["get_genome"]["structure"]["circular"]) > 0:
        return rules.add_overhang_for_circular_chromosomes.output
    else:
        return rules.get_genome.output.fasta


def get_fasta_index(wildcards):
    if len(config["get_genome"]["structure"]["circular"]) > 0:
        return rules.index_genome_with_overhang_chromosomes.output
    else:
        return rules.get_genome.output.fai


# get bam files
def get_bam(wildcards):
    return expand(
        "results/{tool}/align/{sample}/mapped.bam",
        sample=wildcards.sample,
        tool=config["mapping"]["tool"],
    )


def get_bam_2(wildcards):
    if config["mapping"]["umi_tools_dedup"]["enabled"]:
        return f"results/umi_tools/dedup/{wildcards.sample}.bam"
    else:
        return f"results/samtools/sort/{wildcards.sample}.bam"


def get_bai(wildcards):
    if config["mapping"]["umi_tools_dedup"]["enabled"]:
        return f"results/umi_tools/dedup/{wildcards.sample}.bai"
    else:
        return f"results/samtools/sort/{wildcards.sample}.bai"


# get input for multiqc
def get_multiqc_input(wildcards):
    result = []
    result += expand(
        "results/fastqc/{sample}_{read}_fastqc.{ext}",
        sample=samples.index,
        read=["read1", "read2"] if is_paired_end() else ["read1"],
        ext=["html", "zip"],
    )
    if config["processing"]["tool"] == "fastp":
        result += expand(
            "results/fastp/{sample}.json",
            sample=samples.index,
        )
    elif config["processing"]["tool"] == "trim_galore":
        result += expand(
            "results/trim_galore/{sample}_{read}.fastq.gz_trimming_report.txt",
            sample=samples.index,
            read=["read1", "read2"] if is_paired_end() else ["read1"],
        )
    result += expand(
        "results/{tool}/align/{sample}/mapped.bam",
        sample=samples.index,
        tool=config["mapping"]["tool"],
    )
    result += expand(
        "results/rseqc/{tool}/{sample}.txt",
        sample=samples.index,
        tool=["infer_experiment", "bam_stat"],
    )
    result += expand(
        "results/deeptools/coverage/{sample}.bw",
        sample=samples.index,
    )
    return result
