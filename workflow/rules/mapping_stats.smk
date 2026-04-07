rule gffread_gff:
    input:
        fasta=get_genome_for_mapping,
        annotation=rules.get_genome.output.gff,
    output:
        records="results/get_genome/genome.bed",
    threads: 1
    log:
        "results/get_genome/gffread.log",
    message:
        "convert genome annotation from GFF to BED format"
    params:
        extra=config["mapping_stats"]["gffread"]["extra"],
    wrapper:
        "v7.0.0/bio/gffread"


rule rseqc_infer_experiment:
    input:
        aln=get_cram,
        refgene="results/get_genome/genome.bed",
    output:
        "results/rseqc/infer_experiment/{sample}.txt",
    log:
        "results/rseqc/infer_experiment/{sample}.log",
    message:
        "infer experiment type from mapping to features"
    params:
        extra="--sample-size 10000",
    wrapper:
        "v7.0.0/bio/rseqc/infer_experiment"


rule rseqc_bam_stat:
    input:
        get_cram,
    output:
        "results/rseqc/bam_stat/{sample}.txt",
    threads: 2
    params:
        extra="--mapq 5",
    log:
        "results/rseqc/bam_stat/{sample}.log",
    message:
        "collect mapping statistics using RSeQC"
    wrapper:
        "v7.0.0/bio/rseqc/bam_stat"


rule deeptools_coverage:
    input:
        bam=get_cram,
        bai=get_crai,
    output:
        "results/deeptools/coverage/{sample}_{direction}.bw",
    wildcard_constraints:
        direction="for|forward|plus|rev|reverse|minus",
    threads: 4
    params:
        extra=deeptools_extra,
    log:
        "results/deeptools/coverage/{sample}_{direction}.log",
    message:
        "generate normalized coverage using deeptools"
    wrapper:
        "v7.0.0/bio/deeptools/bamcoverage"
