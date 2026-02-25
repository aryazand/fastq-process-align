rule bwa_mem2_index:
    input:
        ref=get_genome_for_mapping,
    output:
        multiext(
            "results/bwa_mem2/index/genome",
            ".0123",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
        ),
    log:
        "results/bwa_mem2/index/genome.log",
    message:
        "build bwa_mem2 index"
    wrapper:
        "v6.2.0/bio/bwa-mem2/index"


rule bwa_mem2:
    input:
        reads=get_processed_fastq,
        idx=rules.bwa_mem2_index.output,
    output:
        "results/bwa_mem2/align/{sample}/mapped.bam",
    log:
        "results/bwa_mem2/align/{sample}/mapped.log",
    message:
        "make bwa_mem2 alignment"
    params:
        extra=config["mapping"]["bwa_mem2"]["extra"],
        sort=config["mapping"]["bwa_mem2"]["sort"],
        sort_order=config["mapping"]["bwa_mem2"]["sort_order"],
        sort_extra=config["mapping"]["bwa_mem2"]["sort_extra"],
    threads: 8
    wrapper:
        "v7.6.0/bio/bwa-mem2/mem"
