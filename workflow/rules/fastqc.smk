rule fastqc:
    input:
        fastq = 'data/raw/{fastq}.fastq.gz',
    output:
        fastqc = "results/fastqc/{fastq}_fastqc.html",
    params:
        outdir = "results/fastqc/",
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 8000
    threads: 2
    shell:
        """
            fastqc {input.fastq} -o {params.outdir} \
        """
