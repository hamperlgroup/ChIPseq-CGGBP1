rule site_matmean:
    input:
        bw="results/covmean/{GID}.coverage.bw",
        bed="results/deseq2_genomewide/DMSO/final_sites.bed",
    output:
        mat="results/matmean_site/{GID}.site_matrix.gz",
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        computeMatrix reference-point \
        --numberOfProcessors {threads} \
        -S {input.bw} -R {input.bed} \
        --referencePoint 'center' \
        --binSize 100 --skipZeros \
        --downstream 10000 --upstream 10000 \
        -o {output.mat}
        """


rule site_heatmapmean:
    input:
        mat="results/matmean_site/{GID}.site_matrix.gz",
    output:
        plot="results/heatmapmean/{GID}.site.png",
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 8000
    threads: 4
    shell:
        """
        plotHeatmap -m {input.mat} \
        --heatmapHeight 8 --heatmapWidth 8 \
        --whatToShow "heatmap and colorbar" \
        --refPointLabel 'site' \
        --regionsLabel 'sites' \
        --colorMap 'YlGnBu' \
        --zMin 0 --zMax 1.5 \
        --xAxisLabel 'distance' \
        --sortRegions 'keep' \
        -o {output.plot}
        """
