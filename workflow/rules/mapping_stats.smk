rule remove_overhang:
    input:
        bam=get_bam,
        fai=get_fasta_index,
        awk=workflow.source_path("../scripts/remove_overhang.awk")
    output:
        temp("results/samtools/{sample}.remove_overhang.bam"),
    log:
        "results/samtools/{sample}.remove_overhang.log",
    message:
        "process overhangs in reads using custom awk script"
    params:
        circular_chroms=",".join(config["get_genome"]["structure"]["circular"].keys()),
        overhang_lengths=",".join(
            f"{chrom}:{overhang}"
            for chrom, overhang in config["get_genome"]["structure"][
                "circular"
            ].items()
        ),
        chromosome_lengths=lambda wildcards, input: get_chrom_lengths_from_fai(
            input.fai[0]
        ),
    shell:
        """
        samtools view -h {input.bam} \
            | gawk -f {input.awk} \
                  -v circular_chroms="{params.circular_chroms}" \
                  -v chromosome_lengths="{params.chromosome_lengths}" \
                  -v overhang_lengths="{params.overhang_lengths}" \
            | samtools view -b -o {output} - 2> {log}
        """


rule samtools_sort:
    input:
        rules.remove_overhang.output,
    output:
        temp("results/samtools/sort/{sample}.bam"),
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
        rules.samtools_sort.output,
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


rule bam_to_cram:
    input:
        bam=rules.samtools_sort.output,
        fa=get_genome_for_mapping,
    output:
        "results/samtools/cram/{sample}.cram",
    log:
        "results/samtools/cram/{sample}.log",
    params:
        extra=lambda input:"-C -T {input.fa}",  # optional params string
        region="",  # optional region string
    threads: 2
    wrapper:
        "v8.1.1/bio/samtools/view"


rule index_cram:
    input:
        rules.bam_to_cram.output,
    output:
        "results/samtools/cram/{sample}.crai",
    log:
        "results/samtools/cram/{sample}.crai.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v8.1.1/bio/samtools/index"


rule umi_tools_dedup:
    input:
        bam=rules.samtools_sort.output,
        bai=rules.samtools_index.output,
    output:
        temp("results/umi_tools/dedup/{sample}.bam"),
    log:
        "results/umi_tools/dedup/{sample}.log",
    message:
        "deduplicate reads using umi_tools"
    params:
        extra=config["mapping"]["umi_tools_dedup"]["extra"],
        paired="--paired" if is_paired_end() else "",
    container:
        "docker://quay.io/biocontainers/umi_tools:1.1.6--py310h1fe012e_0"
    threads: 5
    shell:
        """
        umi_tools dedup \
            -I {input} \
            -S {output} \
            --log={log} \
            {params.paired} \
            {params.extra} 
        """


rule bam_to_cram_dedup:
    input:
        bam=rules.umi_tools_dedup.output,
        ref="results/get_genome/genome.fasta",
    output:
        "results/umi_tools/dedup/{sample}.cram",
    log:
        "results/umi_tools/dedup/{sample}.cram.log",
    params:
        extra=lambda wildcards, input:"-C",  # optional params string
        region="",  # optional region string
    threads: 2
    wrapper:
        "v8.1.1/bio/samtools/view"


rule index_cram_dedup:
    input:
        rules.bam_to_cram_dedup.output,
    output:
        "results/umi_tools/dedup/{sample}.crai",
    log:
        "results/umi_tools/dedup/{sample}.crai.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v8.1.1/bio/samtools/index"


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
        "results/deeptools/coverage/{sample}.bw",
    wildcard_constraints:
        sample="|".join(map(re.escape, samples.index)),
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
