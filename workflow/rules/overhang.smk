rule add_overhang_for_circular_chromosomes:
    input:
        "results/get_genome/genome.fasta",
    output:
        fasta="results/get_genome/genome_with_overhang.fasta",
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
