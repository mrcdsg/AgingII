configfile: "config1.yaml"

SAMPLES = ["sample"]["single"]
GENOME = "reference_genome/mm10"

rule all:
    input:
        expand("results/fastqc/{sample}_fastqc.html", sample=SAMPLES),
        expand("results/fastqc/{sample}_fastqc.zip", sample=SAMPLES),
        expand("results/trimmed/{sample}.fastq.gz", sample=SAMPLES),
        expand("results/aligned/{sample}.bam", sample=SAMPLES),
        expand("results/aligned/{sample}.bam.bai", sample=SAMPLES),
        expand("results/filtered/{sample}_dedup.bam", sample=SAMPLES),
        "results/multiqc_report.html"

rule fastqc:
    input:
        "data/{sample}.fastq.gz"
    output:
        html="results/fastqc/{sample}_fastqc.html",
        zip="results/fastqc/{sample}_fastqc.zip"
    threads: 2
    shell:
        "mkdir -p results/fastqc && fastqc -t {threads} {input} --outdir results/fastqc"

rule trim:
    input:
        "data/{sample}.fastq.gz"
    output:
        "results/trimmed/{sample}.fastq.gz"
    log:
        "logs/snakemake/trimmomatic/{sample}.log"
    shell:
        """
        mkdir -p results/trimmed logs/snakemake/trimmomatic && \
        trimmomatic SE -phred33 {input} {output} SLIDINGWINDOW:4:20 MINLEN:25 > {log} 2>&1
        """
rule align:
    input:
        "results/trimmed/{sample}.fastq.gz"
    output:
        bam="results/aligned/{sample}.bam"
        bai="results/aligned/{sample}.bam.bai"
    params:
        index=GENOME
    log:
        "logs/snakemake/bowtie2/{sample}.log"
    threads: 8
    shell:
        """
        mkdir -p results/aligned logs/snakemake/bowtie2 && \bowtie2 -x {params.index} -U {input} -p {threads} 2>> {log} | \samtools view -bS - | \samtools sort -o {output.bam} - && \
        samtools index {output.bam}
        """

rule dedup:
    input:
        "results/aligned/{sample}.bam"
    output:
        "results/filtered/{sample}_dedup.bam"
    log:
        "logs/snakemake/picard/{sample}.log"
    shell:
        """
        mkdir -p results/filtered logs/snakemake/picard && \ picard MarkDuplicates I={input} O={output} REMOVE_DUPLICATES=true M=logs/snakemake/picard/{wildcards.sample}_metrics.txt > {log} 2>&1
        """
rule call_peaks:
    input:
        "results/filtered/{sample}_dedup.bam"
    output:
        "results/peaks/{sample}_peaks.narrowPeak"
    log:
        "logs/snakemake/macs2/{sample}.log"
    shell:
        """
        mkdir -p results/peaks logs/snakemake/macs2 && \macs2 callpeak -t {input} -f BAM -g mm -n {wildcards.sample} --outdir results/peaks > {log} 2>&1
        """
rule multiqc:
    input:
        expand("results/fastqc/{sample}_fastqc.zip", sample=SAMPLES),
        expand("logs/snakemake/trimmomatic/{sample}.log", sample=SAMPLES),
        expand("logs/snakemake/bowtie2/{sample}.log", sample=SAMPLES),
        expand("logs/snakemake/macs2/{sample}.log", sample=SAMPLES),
        expand("logs/snakemake/picard/{sample}.log", sample=SAMPLES)
    output:
        "results/multiqc_report.html"
    shell:
        """
        mkdir -p results/ && multiqc results/ -o results/
        """



