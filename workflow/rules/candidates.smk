

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


rule split_candidates_bed_avg:
    input:
        "resources/{chromosome}/bed_avg.bedGraph",
    output:
        scatter.split_candidates(
            "resources/{{chromosome}}/bed_avg_{scatteritem}.bedGraph"
        ),
    log:
        "logs/split_candidates_{chromosome}.log",
    conda:
        "../envs/rbt.yaml"
    params:
        num_parts=config["num_candidates"],
    shell:
        """
        total_lines=$(wc -l < {input})
        lines_per_file=$((total_lines / {params.num_parts} + 1))
        split -l $lines_per_file -d -a 2 {input} resources/{wildcards.chromosome}/bed_avg_
        total_parts={params.num_parts}
        for i in $(seq -f "%02g" 0 $((total_parts - 1))); do
            new_i=$(printf "%d" $((10#$i + 1)))
            mv resources/{wildcards.chromosome}/bed_avg_$i resources/{wildcards.chromosome}/bed_avg_${{new_i}}-of-${{total_parts}}.bedGraph
        done
        """


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
        candidates="resources/{chromosome}/bed_avg_{scatteritem}.bedGraph",
        # candidates="resources/candidates_{scatteritem}.vcf",
        mutations="resources/HG002_GRCh38_1_22_v4.2.1_benchmark_filtered.vcf",
        consistent="resources/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent_filtered.bed",
    output:
        "resources/{chromosome}/bed_avg_{scatteritem}_filtered.bedGraph",
    conda:
        "../envs/bedtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
        deamination_path=config["deamination_path"],
        filtered=config["filter"],
    log:
        "logs/filter_candidates_{chromosome}_{scatteritem}.log",
    shell:
        """
        echo {params.filtered}
        if [ "{params.filtered}" = "True" ]; then
            echo "Filtering candidates"
            bedtools intersect -a {input.candidates} -b {input.consistent} > {output}
        else
            cp {input.candidates} {output}
        fi
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
