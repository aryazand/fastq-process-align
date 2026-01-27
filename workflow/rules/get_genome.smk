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


rule add_overhang_for_circular_chromosomes:
    input:
        "results/get_genome/genome.fasta",
    output:
        "results/get_genome/genome_with_overhang.fasta",
    params:
        overhang=config["get_genome"]["structure"]["circular"],
    container:
        "docker://quay.io/biocontainers/bioconductor-biostrings:2.74.0--r44h3df3fcb_1"
    script:
        "../scripts/add_overhang.R"


rule index_genome_with_overhang_chromosomes:
    input:
        "results/get_genome/genome_with_overhang.fasta",
    output:
        "results/get_genome/genome_with_overhang.fasta.fai",
    log:
        "results/get_genome/genome_with_overhang.faidx.log",
    params:
        extra="",
    wrapper:
        "v8.1.1/bio/samtools/faidx"
