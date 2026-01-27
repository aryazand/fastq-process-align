rule star_index:
    input:
        fasta=get_genome_for_mapping,
    output:
        directory("results/star/index/"),
    threads: 1
    params:
        extra=config["mapping"]["star"]["index"],
    log:
        "results/star/index/index.log",
    message:
        "build star index"
    wrapper:
        "v7.2.0/bio/star/index"


rule star_align:
    input:
        fq1=lambda wildcards: get_processed_fastq(wildcards.sample, regex="read1"),
        fq2=lambda wildcards: (
            get_processed_fastq(wildcards.sample, regex="read2")
            if is_paired_end()
            else []
        ),
        idx=rules.star_index.output,
    output:
        aln="results/star/align/{sample}/mapped.bam",
        log_final="results/star/align/{sample}/Log.final.out",
        sj="results/star/align/{sample}/SJ.out.tab",
    log:
        "results/star/align/{sample}/mapped.log",
    message:
        "make star alignment"
    params:
        extra=config["mapping"]["star"]["extra"],
    threads: 8
    wrapper:
        "v7.2.0/bio/star/align"
