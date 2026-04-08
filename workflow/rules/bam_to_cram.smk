rule samtools_sort:
    input:
        get_bam,
    output:
        temp("results/alignment/{sample}_sorted.bam"),
    log:
        "results/alignment/{sample}_samtools_sort.log",
    message:
        "re-sort reads after mapping regardless if mapper did"
    params:
        extra=config["mapping"]["samtools_sort"]["extra"],
    threads: 2
    wrapper:
        "v7.0.0/bio/samtools/sort"


rule bam_to_cram:
    input:
        bam=rules.samtools_sort.output,
        fa=get_genome_for_mapping,
    output:
        "results/alignment/{sample}.cram",
    log:
        "results/alignment/{sample}_cram.log",
    params:
        extra=lambda wildcards, input: f"-C -T {input.fa}",  # optional params string
        region="",  # optional region string
    threads: 2
    wrapper:
        "v9.4.1/bio/samtools/view"


rule index_cram:
    input:
        rules.bam_to_cram.output,
    output:
        "results/alignment/{sample}.crai",
    log:
        "results/alignment/{sample}.crai.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v8.1.1/bio/samtools/index"
