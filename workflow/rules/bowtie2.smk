rule bowtie2_index:
    input:
        "data/genome/genome.fa",
    output:
        expand("data/genome/genome{ext}", ext = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]),
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        bowtie2-build --threads {threads} {input} data/genome/genome
        """

rule bowtie2_align:
    input:
        trimmed1=lambda wildcards: expand("results/trimmed/{sample}_1_val_1.fq.gz", sample = SampleTable["sample_name"][SampleTable["sample_ID"] == wildcards.ID]),
        trimmed2=lambda wildcards: expand("results/trimmed/{sample}_2_val_2.fq.gz", sample = SampleTable["sample_name"][SampleTable["sample_ID"] == wildcards.ID]),
        index=rules.bowtie2_index.output,
    output:
        bam="results/BAM/align/{ID}.bam",
    params:
        bowtie2="--end-to-end --very-sensitive --no-unal --no-mixed --no-discordant --dovetail -I 10 -X 700",
    conda:
        "../envs/env_preprocessing.yaml",
    log:
        "results/logs/bowtie2/{ID}.log",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        bowtie2 -x data/genome/genome --threads {threads} {params.bowtie2} -1 {input.trimmed1} -2 {input.trimmed2} 2> {log} | samtools view -Sbh -o {output.bam}
        """