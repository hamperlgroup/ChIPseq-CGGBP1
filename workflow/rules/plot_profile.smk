rule scaled_matrix:
    input:
        bw="results/coverage/{ID}.coverage.bw",
        gtf="data/genome/combined.gtf",
    output:
        mat="results/matrix_scaled/{ID}.scaled_matrix.gz",
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        computeMatrix scale-regions \
        --numberOfProcessors {threads} \
        -S {input.bw} -R {input.gtf} \
        --binSize 20 --skipZeros \
        --downstream 2000 --upstream 2000 \
        --regionBodyLength 4000 \
        -o {output.mat}
        """

rule scaled_plot:
    input:
        mat="results/matrix_scaled/{ID}.scaled_matrix.gz",
    output:
        plot="results/metaplots/{ID}.scaled_plot.png",
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 8000
    threads: 4
    shell:
        """
        plotProfile -m {input.mat} \
        --plotHeight 10 --plotWidth 14 \
        --yMin 0.05 --yMax 1.0 \
        -o {output.plot}
        """
