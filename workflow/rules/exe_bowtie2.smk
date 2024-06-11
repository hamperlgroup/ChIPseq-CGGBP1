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

def get_list(my_id, my_columns):
    out = [i.replace('.fastq.gz','') for i in SampleTable[SampleTable['id']==my_id][my_columns].values.flatten().tolist()]
    return out

rule bowtie2_align:
    input:
        trimmed1=lambda wildcards: expand("results/trimmed/{LID}_val_1.fq.gz", LID=get_list(wildcards.ID, ['fastq_1','fastq_3'])),
        trimmed2=lambda wildcards: expand("results/trimmed/{LID}_val_2.fq.gz", LID=get_list(wildcards.ID, ['fastq_2','fastq_4'])),
        index=rules.bowtie2_index.output,
    output:
        bam="results/BAM/raw/{ID}.bam",
    params:
        trimmed1=lambda wildcards: ",".join(expand("results/trimmed/{LID}_val_1.fq.gz", LID=get_list(wildcards.ID, ['fastq_1','fastq_3']))),
        trimmed2=lambda wildcards: ",".join(expand("results/trimmed/{LID}_val_2.fq.gz", LID=get_list(wildcards.ID, ['fastq_2','fastq_4']))),
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
        bowtie2 -x data/genome/combined --threads {threads} {params.bowtie2} -1 {params.trimmed1} -2 {params.trimmed2} 2> {log} | samtools view -Sbh -o {output.bam}
        """
