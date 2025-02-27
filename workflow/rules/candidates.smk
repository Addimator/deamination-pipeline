

rule find_candidates:
    input:
        "resources/chromosome.fasta",
    output:
        "resources/candidates.bcf",
    log:
        "logs/find_candidates.log",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        varlo_path=config["varlo_path"],
        pipeline_path=config["pipeline_path"],
    shell:
        """ 
        cd {params.varlo_path}
        cargo run -- methylation-candidates {params.pipeline_path}{input} {params.pipeline_path}{output}
        """


rule candidates_to_vcf:
    input:
        "resources/candidates.bcf",
    output:
        "resources/candidates.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/convert_to_vcf.log",
    shell:
        """
        bcftools view {input} > {output}
        """


rule split_candidates:
    input:
        "resources/candidates.bcf",
    output:
        scatter.split_candidates("resources/candidates_{scatteritem}.bcf"),
    log:
        "logs/split_candidates.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output}"


rule candidate_splits_to_vcf:
    input:
        "resources/candidates_{scatteritem}.bcf",
    output:
        "resources/candidates_{scatteritem}.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/convert_splits_to_vcf_{scatteritem}.log",
    shell:
        """
        bcftools view {input} > {output}
        """


# Just for debugging reasons
rule get_candidate_neighborhood:
    input:
        genome="resources/genome.fasta",
    output:
        "resources/candidate_neighborhood.txt",
    params:
        pipeline_path=config["pipeline_path"],
        chrom=21,
        candidate=candidate,
        left_border=candidate - 150,
        right_border=candidate + 150,
        candidate_left=candidate - 1,
        candidate_right=candidate + 1,
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools faidx resources/genome.fasta "{params.chrom}:{params.left_border}-{params.candidate_left}" && samtools faidx resources/genome.fasta "{params.chrom}:{params.candidate}-{params.candidate}" && samtools faidx resources/genome.fasta "{params.chrom}:{params.candidate_right}-{params.right_border}"
        """


rule filter_candidates:
    input:
        candidates="resources/candidates_{scatteritem}.vcf",
        mutations="resources/HG002_GRCh38_1_22_v4.2.1_benchmark_filtered.vcf",
        consistent="resources/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent_filtered.bed",
    output:
        "resources/candidates_filtered_{scatteritem}.bcf",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        pipeline_path=config["pipeline_path"],
        deamination_path=config["deamination_path"],
    log:
        "logs/filter_candidates{scatteritem}.log",
    shell:
        """
        cd {params.deamination_path}   
        cargo run  -- filter-candidates {params.pipeline_path}{input.candidates} {params.pipeline_path}{input.mutations} {params.pipeline_path}{input.consistent} {params.pipeline_path}{output}
        """


rule candidate_filtered_to_vcf:
    input:
        "resources/candidates_filtered_{scatteritem}.bcf",
    output:
        "resources/candidates_filtered_{scatteritem}.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/candidate_filtered_to_vcf{scatteritem}.log",
    shell:
        """
        bcftools view {input} > {output}
        """
