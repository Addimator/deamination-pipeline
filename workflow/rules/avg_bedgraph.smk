chr_chromosome = "chr" + str(config["chromosome"])


rule download_bedGraphs:
    output:
        "resources/HG002/{bedGraph}.bedGraph.gz",
    params:
        pipeline_path=config["pipeline_path"],
        bedGraphs=config["bedGraphs_HG002"],
        bedgraph_path="resources/HG002",
    script:
        "../scripts/get_bedGraph_data.py"


ruleorder: filter_bedGraphs > extract_data


rule extract_data:
    input:
        "resources/HG002/{bedGraph}.bedGraph.gz",
    output:
        "resources/HG002/{bedGraph}.bedGraph",
    shell:
        "gunzip -c {input} > {output}"


rule filter_bedGraphs:
    input:
        "resources/HG002/{bedGraph}.bedGraph",
    output:
        "resources/HG002/{bedGraph}-{chromosome}.bedGraph",
    wildcard_constraints:
        chromosome="[^_]+",
    shell:
        """
        awk '$1 == "{wildcards.chromosome}" {{print}}' {input} > {output}
        """


rule compute_avg_bedGraph:
    input:
        expand(
            "resources/HG002/{bedGraph}-{chromosome}.bedGraph",
            bedGraph=config["bedGraphs_HG002"],
            chromosome=chr_chromosome,
        ),
    output:
        "resources/bed_avg_{chromosome}.bedGraph",
    wildcard_constraints:
        chromosome="[^_]+",
    params:
        chromosome=lambda wildcards: "chr" + wildcards.chromosome,
    script:
        "../scripts/compute_avg_bedGraph.py"


# rule plot_avg_bedGraph:
#     input:
#         candidates=expand(
#             "resources/{chrom}/candidates.vcf", chrom=[chrom for chrom in chromosomes]
#         ),
#         bedgraphs=expand(
#             "resources/HG002/{bedGraph}-{chromosome}.bedGraph",
#             bedGraph=config["bedGraphs_HG002"],
#             chromosome=[chrom for chrom in chr_chromosome],
#         ),
#     output:
#         cov="resources/bed_avg_cov.{plot_type}",
#         meth="resources/bed_avg_meth.{plot_type}",
#     conda:
#         "../envs/plot.yaml"
#     wildcard_constraints:
#         chromosome="[^_]+",
#     log:
#         "logs/plot_avg_bedGraph.log",
#     params:
#         plot_type=config["plot_type"],
#     script:
#         "../scripts/plot_avg_bedGraph.py"
