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
