rule bowtie2_index:
    input:
        "data/genome/combined.fa",
    output:
        expand("data/genome/combined{ext}", ext = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]),
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        bowtie2-build --threads {threads} {input} data/genome/combined
        """

rule bowtie2_align:
    input:
        trimmed1="results/trimmed/{ID}_1_val_1.fq.gz",
        trimmed2="results/trimmed/{ID}_2_val_2.fq.gz",
        index=rules.bowtie2_index.output,
    output:
        bam="results/BAM/{ID}.bam",
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
        bowtie2 -x data/genome/combined --threads {threads} {params.bowtie2} -1 {input.trimmed1} -2 {input.trimmed2} 2> {log} | samtools view -Sbh -o {output.bam}
        """
