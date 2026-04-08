#################################################################
# PROCESS ALIGNEMNTS
# 1. Filter alignments per user specifications (e.g. MAPQ, flags, etc.)
# 2. Process overhangs in reads (e.g. for circular genomes)
# 3. Deduplicate reads using umi_tools (if UMIs are present)
# 4. Convert to CRAM and index
#################################################################


rule samtools_filter:
    input:
        cram=rules.bam_to_cram.output,
        crai=rules.index_cram.output,
        fa=get_genome_for_mapping,
    output:
        bam="results/processed_alignment/samtools_filter/{sample}.bam",
    log:
        view="results/processed_alignment/samtools_filter/{sample}_view.log",
        reaheader="results/processed_alignment/samtools_filter/{sample}_reheader.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"
    params:
        extra=config["mapping"]["postprocessing"]["filter"],  # optional params string
        region=config["mapping"]["postprocessing"]["region"],  # optional region string
    threads: 2
    shell:
        """
        # 1. samtools view extracts the target chromosome and includes the header (-h)
        # 2. awk processes the SAM stream line by line
        # 3. The output is piped back into samtools to compress into the final BAM
        
        samtools view -h {input.cram} {params.region} | \
        awk -F'\\t' -v OFS='\\t' -v chrom="{params.region}" '
        /^@SQ/ {{
            # Split the second column (e.g., SN:chr1) by colon
            split($2, a, ":");
            # Only print the header line if the chromosome matches our target
            if (a[2] == chrom) print;
            next;
        }}
        /^@/ {{ 
            # Print all other header lines (e.g., @HD, @RG, @PG)
            print; 
            next; 
        }}
        {{
            # For alignments: column 7 ($7) holds the mate reference sequence name.
            # "=" means the mate is on the exact same chromosome.
            # "*" means the read is single-end (or mate is unmapped).
            # chrom catches edge cases where the aligner writes the explicit name instead of "=".
            if ($7 == "=" || $7 == chrom || $7 == "*") {{
                print;
            }}
        }}' | \
        samtools view -@ {threads} -O BAM -o {output.bam}

        # samtools view -bh {params.extra} -T {input.fa} -@ {threads} -o {output.bam}.temp {input.cram} {params.region} 2> {log.view}
        # samtools view -h {output.bam}.temp | awk '$1 == "@SQ" && $2 != "SN:{params.region}" {{ next }}{{ print }}' | samtools reheader - {output.bam}.temp > {output.bam} 2> {log.reaheader} 
        # rm {output.bam}.temp
        """


rule remove_overhang:
    # Needs to be updated with clearer code
    # Needs a container
    input:
        bam=get_removed_overhang_input,
        fai=get_fasta_index,
        adjust_overhang_awk=workflow.source_path("../scripts/adjust_overhang.awk"),
    output:
        "results/processed_alignment/remove_overhang/{sample}.bam",
    log:
        sam="results/processed_alignment/remove_overhang/{sample}.log",
        stats="results/processed_alignment/remove_overhang/{sample}_overhang_stats.txt",
    message:
        "process overhangs in reads using custom awk script"
    params:
        overhang_lengths=config["get_genome"]["structure"]["circular"][
            "overhang_length"
        ],
        chromosome_lengths=lambda wildcards, input: get_chrom_lengths_from_fai(
            input.fai[0]
        ),
    shell:
        """
        samtools view -h {input.bam} | \
            awk -f {input.adjust_overhang_awk} \
                -v chrom_lengths="{params.chromosome_lengths}" \
                2> {log.stats} | \
            samtools view -bh -o {output} - 2> {log.sam}
        """


rule index_umitools_input:
    input:
        bam=get_umi_tools_dedup_input,
    output:
        "results/processed_alignment/index_prior_to_dedup/{sample}.bam.bai",
    log:
        "results/processed_alignment/index_prior_to_dedup/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v9.4.1/bio/samtools/index"


rule umi_tools_dedup:
    input:
        bam=get_umi_tools_dedup_input,
        bai=rules.index_umitools_input.output,
    output:
        temp("results/processed_alignment/dedup/{sample}.bam"),
    log:
        "results/processed_alignment/dedup/{sample}.log",
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


rule bam_to_cram_post_processing:
    input:
        bam=bam_to_cram_post_processing_input,
        fa=get_genome_for_mapping,
    output:
        "results/processed_alignment/cram/{sample}.cram",
    log:
        "results/processed_alignment/cram/{sample}.cram.log",
    params:
        extra=lambda wildcards, input: f"-C -T {input.fa}",  # optional params string
        region="",  # optional region string
    threads: 2
    wrapper:
        "v8.1.1/bio/samtools/view"


rule index_cram_post_processing:
    input:
        rules.bam_to_cram_post_processing.output,
    output:
        "results/processed_alignment/cram/{sample}.crai",
    log:
        "results/processed_alignment/cram/{sample}.crai.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v8.1.1/bio/samtools/index"
