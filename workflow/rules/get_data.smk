### Get genome and filter to chromosome
rule get_genome:
    output:
        "resources/genome.fasta",
    params:
        species=genome_config["species"],
        datatype=genome_config["datatype"],
        build=genome_config["build"],
        release=genome_config["release"],
    log:
        "logs/get_genome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


rule genome_index:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome_index.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """ 
        samtools faidx {params.pipeline_path}{input}
        """


rule filter_genome:
    input:
        "resources/genome.fasta",
    output:
        "resources/chromosome.fasta",
    log:
        "logs/filter_genome.log",
    conda:
        "../envs/samtools.yaml"
    params:
        chromosome=chromosome_conf["chromosome"],
    threads: 10
    shell:
        """ 
        samtools faidx {input} {params.chromosome} > {output}
        """


### Giab annotations for filtering regions
rule download_giab_annotations:
    output:
        consistent="resources/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
        mutations="resources/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
    params:
        consistent="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
        mutations="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
    shell:
        """
        wget -O {output.consistent} {params.consistent}
        wget -O {output.mutations} {params.mutations}
        """


rule filter_giab_annotations:
    input:
        consistent="resources/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
        mutations="resources/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
    output:
        consistent="resources/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent_filtered.bed",
        mutations="resources/HG002_GRCh38_1_22_v4.2.1_benchmark_filtered.vcf",
    params:
        chromosome="chr" + chromosome_conf["chromosome"],
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bcftools index {input.mutations}
        bcftools view -r {params.chromosome} {input.mutations} -o {output.mutations}
        awk '$1 == "{params.chromosome}"' {input.consistent} > {output.consistent}
        """


### Get reads, difference between single end, paired end and pcr_free reads

# rule get_fastq_pe:
#     output:
#         # the wildcard name must be accession, pointing to an SRA number
#         "resources/{SRA}/{accession}_1.fastq",
#         "resources/{SRA}/{accession}_2.fastq",
#     log:
#         "logs/pe/get_fastq_pe{SRA}_{accession}.log",
#     params:
#         extra="--skip-technical",
#     threads: 6  # defaults to 6
#     conda:
#         "envs/fastq-wrapper.yaml"
#     resources:
#         mem_mb=512,
#     wrapper:
#         "v2.6.0/bio/sra-tools/fasterq-dump"


# rule align_reads_se:
#     input:
#         fasta="resources/genome.fasta",
#         reads1="resources/{SRA}/{accession}.fastq",
#     output:
#         "resources/{SRA}{accession}/aligned-reads.bam",
#     conda:
#         "envs/bwa-mem2.yaml"
#     log:
#         "logs/align_reads_{SRA}.log",
#     shell:
#         """
#         bwa-mem2 index {input.fasta}
#         bwa-mem2 mem -t 8 {input.fasta} {input.reads1} | samtools view -S -b - > {output}
#         """


rule get_pcr_free_data_forward:
    output:
        "resources/{SRA}/{accession}_R1.fastq.gz",
    log:
        "logs/get_pcr_free_data_forward{SRA}/{accession}.log",
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        "wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R1.fastq.gz -O {params.pipeline_path}{output}"


rule get_pcr_free_data_reverse:
    output:
        "resources/{SRA}/{accession}_R2.fastq.gz",
    log:
        "logs/get_pcr_free_data_reverse{SRA}/{accession}.log",
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        "wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R2.fastq.gz -O {params.pipeline_path}{output}"


rule decompress_fastq_reverse:
    input:
        "resources/{SRA}/{accession}_{number}.fastq.gz",
    output:
        "resources/{SRA}/{accession}_{number}.fastq",
    shell:
        "gunzip -c {input} > {output}"
