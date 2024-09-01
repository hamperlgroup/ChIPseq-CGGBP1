rule trim_galore:
    input:
        fastq1="data/raw/{sample}_1.fastq.gz",
        fastq2="data/raw/{sample}_2.fastq.gz"
    output:
        trimmed1="results/trimmed/{sample}_1_val_1.fq.gz",
        trimmed2="results/trimmed/{sample}_2_val_2.fq.gz",
    params:
        dir="results/trimmed/"
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 8000
    threads: 8
    shell:
        """
        trim_galore -j {threads} --quality 28 --paired {input.fastq1} {input.fastq2} -o {params.dir}
        """