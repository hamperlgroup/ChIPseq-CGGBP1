rule deseq2:
    input:
        table="data/Design_file.csv",
        counts="results/counts/counts_per_bin.tab",
    output:
        bed="results/deseq2_genomewide/{TID}/sessionInfo.txt",
        session="results/deseq2_genomewide/{TID}/final_sites.bed",
    params:
    conda:
        "../envs/env_deseq2.yaml",
    resources:
        mem_mb = 24000
    threads: 4
    shell:
        """
        Rscript --vanilla workflow/scripts/deseq2.R {input} {params} {output}
        """

rule deseq2_insertions:
    input:
        table="data/Design_file.csv",
        counts="results/counts/counts_per_bin.tab",
    output:
        session="results/deseq2_insertions/{TID}/sessionInfo.txt",
    params:
    conda:
        "../envs/env_deseq2.yaml",
    resources:
        mem_mb = 24000
    threads: 4
    shell:
        """
        Rscript --vanilla workflow/scripts/deseq2_insertions.R {input} {params} {output}
        """

rule deseq2_comparison:
    input:
        table="data/Design_file.csv",
        counts="results/counts/counts_per_bin.tab",
        res1 = "results/deseq2_genomewide/DMSO/sessionInfo.txt",
        res2 = "results/deseq2_insertions/DMSO/sessionInfo.txt",
    output:
        session="results/deseq2_comparison/sessionInfo.txt",
    params:
    conda:
        "../envs/env_deseq2.yaml",
    resources:
        mem_mb = 24000
    threads: 4
    shell:
        """
        Rscript --vanilla workflow/scripts/deseq2_comparison.R {input} {params} {output}
        """
