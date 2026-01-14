rule samtools_sort:
    input:
        get_bam,
    output:
        "results/samtools/sort/{sample}.bam",
    log:
        "results/samtools/sort/{sample}.log",
    message:
        "re-sort reads after mapping regardless if mapper did"
    params:
        extra=config["mapping"]["samtools_sort"]["extra"],
    threads: 2
    wrapper:
        "v7.0.0/bio/samtools/sort"


rule samtools_index:
    input:
        "results/samtools/sort/{sample}.bam",
    output:
        "results/samtools/sort/{sample}.bai",
    log:
        "results/samtools/sort/{sample}_index.log",
    message:
        "index reads"
    params:
        extra=config["mapping"]["samtools_index"]["extra"],
    threads: 2
    wrapper:
        "v7.0.0/bio/samtools/index"


rule gffread_gff:
    input:
        fasta=rules.get_genome.output.fasta,
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


rule umi_tools_dedup:
    input:
        bam="results/samtools/sort/{sample}.bam",
        bai="results/samtools/sort/{sample}.bai",
    output:
        bam="results/umi_tools/dedup/{sample}.bam",
    log:
        "results/umi_tools/dedup/{sample}.log",
    message:
        "deduplicate reads using umi_tools"
    params:
        extra=config["mapping"]["umi_tools_dedup"]["extra"],
        paired="--paired" if is_paired_end() else "",
    conda:
        "../envs/umi_tools.yml"
    threads: 1
    shell:
        """
        umi_tools dedup \
            -I {input.bam} \
            -S {output.bam} \
            --log={log} \
            {params.paired} \
            {params.extra} 
        """


rule samtools_index_dedup:
    input:
        "results/umi_tools/dedup/{sample}.bam",
    output:
        "results/umi_tools/dedup/{sample}.bai",
    log:
        "results/umi_tools/dedup/{sample}_index.log",
    message:
        "index deduplicated reads"
    params:
        extra=config["mapping"]["samtools_index"]["extra"],
    threads: 2
    wrapper:
        "v7.0.0/bio/samtools/index"


rule rseqc_infer_experiment:
    input:
        aln="results/samtools/sort/{sample}.bam",
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
        "results/samtools/sort/{sample}.bam",
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
        bam=get_bam_2,
        bai=get_bai,
    output:
        "results/deeptools/coverage/{sample}.bw",
    threads: 4
    params:
        effective_genome_size=config["mapping_stats"]["deeptools_coverage"][
            "genome_size"
        ],
        extra=config["mapping_stats"]["deeptools_coverage"]["extra"],
    log:
        "results/deeptools/coverage/{sample}.log",
    message:
        "generate normalized coverage files using deeptools"
    wrapper:
        "v7.0.0/bio/deeptools/bamcoverage"
