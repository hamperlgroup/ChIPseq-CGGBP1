rule samtools:
    input:
        bam="results/BAM/align/{ID}.bam",
    output:
        flt="results/BAM/flt/{ID}.flt.bam",
        bai="results/BAM/flt/{ID}.flt.bam.bai",
    params:
        mapqc = 12
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        samtools view -bS -@ {threads} -q {params.mapqc} {input.bam} | \
        samtools sort -@ {threads} - | tee {output.flt} | \
        samtools index - {output.bai}
        """


rule picard:
    input:
        flt="results/BAM/flt/{ID}.flt.bam",
        bai="results/BAM/flt/{ID}.flt.bam.bai",
    output:
        rmdup="results/BAM/rmdup/{ID}.rmdup.bam",
        bai="results/BAM/rmdup/{ID}.rmdup.bam.bai",
        idx="results/logs/picard/{ID}_idx.txt",
    params:
        rm = "TRUE"
    conda:
        "../envs/env_preprocessing.yaml",
    log:
        "results/logs/picard/{ID}.log",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        picard MarkDuplicates -I {input.flt} -O {output.rmdup} \
        -CREATE_INDEX FALSE -REMOVE_DUPLICATES {params.rm} \
        -METRICS_FILE {log}

        samtools index {output.rmdup}
        samtools idxstats {output.rmdup} | grep -v "K" | grep -v "G" | grep -v "\*" > {output.idx}
        """
