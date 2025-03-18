rule align_reads_pe:
    input:
        fasta="resources/genome.fasta",
        reads1="resources/{SRA}/{SRA}_R1.fastq",
        reads2="resources/{SRA}/{SRA}_R2.fastq",
    output:
        "resources/{SRA}/aligned-reads.bam",
    conda:
        "../envs/bwa-mem2.yaml"
    log:
        "logs/align_reads{SRA}.log",
    resources:
        mem_mb=1024,
    shell:
        """
        bwa-mem2 index {input.fasta}
        bwa-mem2 mem -t 8 {input.fasta} {input.reads1} {input.reads2} | samtools view -S -b - > {output}
        """


rule sort_aligned_reads:
    input:
        "resources/{SRA}/aligned-reads.bam",
    output:
        "resources/{SRA}/aligned-reads-sorted.bam",
    log:
        "logs/sort_aligned_reads{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """
        samtools sort -@ {threads}  {input} -o {output}    
        """


# Focus on chromosome
rule filter_aligned_reads:
    input:
        "resources/{SRA}/aligned-reads-sorted.bam",
    output:
        bam="resources/{SRA}/aligned-reads-focused.bam",
    log:
        "logs/filter_aligned_reads{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
        chromosome=chromosome,
    threads: 10
    shell:
        """ 
        samtools index -@ {threads} {params.pipeline_path}{input}
        if [ "{params.chromosome}" != "all" ]; then
            samtools view -b -o {output.bam} {input} {params.chromosome}
        else
            samtools view -b -o {output.bam} {input}
        fi
        """


####################################


rule filter_mapping_quality:
    input:
        "resources/{SRA}/aligned-reads-focused.bam",
    output:
        "resources/{SRA}/aligned-reads-focused_filtered.bam",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
        min_quality=10,
    threads: 10
    shell:
        """
        samtools view -q {params.min_quality} -b -o {output} {input}
        """


rule markduplicates_bam:
    input:
        bams="resources/{SRA}/aligned-reads-focused_filtered.bam",
    output:
        bam="resources/{SRA}/aligned-reads-focused_dedup.bam",
        metrics="resources/{SRA}/alignment_focused_dedup.metrics.txt",
    log:
        "logs/markduplicates_bam_{SRA}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.6.0/bio/picard/markduplicates"


rule index_markduplicated_bam:
    input:
        bam="resources/{SRA}/aligned-reads-focused_dedup.bam",
    output:
        bam="resources/{SRA}/aligned-reads-focused_dedup.bam.bai",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"


rule aligned_reads_candidates_region:
    input:
        alignment="resources/{SRA}/aligned-reads-focused_dedup.bam",
        index="resources/{SRA}/aligned-reads-focused_dedup.bam.bai",
        candidate=expand(
            "resources/{chromosome}/bed_avg_{{scatteritem}}_filtered.bedGraph",
            chromosome=chromosome,
        ),
        # candidate="resources/bed_avg_{chromosome}_{scatteritem}_filtered.bedGraph",
    output:
        "resources/{SRA}/candidate_specific/alignment_{scatteritem}.bam",
    conda:
        "../envs/samtools.yaml"
    params:
        window_size=config["max_read_length"],
        chromosome=chromosome,
    shell:
        """
        samtools view -b -L {input.candidate} {input.alignment} > {output}

        """
        # set +o pipefail;
        # start=$(awk 'NR==1 {{print $2}}' {input.candidate})
        # end=$(awk 'END {{print $2}}' {input.candidate})
        # end=$((end + {params.window_size}))
        # samtools view -b {input.alignment} "{params.chromosome}:$start-$end" > {output}


rule aligned_reads_candidates_region_to_sam:
    input:
        "resources/{SRA}/candidate_specific/alignment_{scatteritem}.bam",
    output:
        "resources/{SRA}/candidate_specific/alignment_{scatteritem}.sam",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -o {output} {input}
        """


rule aligned_reads_candidates_region_index:
    input:
        "resources/{SRA}/candidate_specific/alignment_{scatteritem}.bam",
    output:
        "resources/{SRA}/candidate_specific/alignment_{scatteritem}.bam.bai",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {input}
        """
