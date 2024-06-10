rule site_matdiff:
    input:
        bw="results/covdiff/{GID1}_vs_{GID2}.coverage.bw",
        bed="results/deseq2_genomewide/DMSO/final_sites.bed",
    output:
        mat="results/matdiff_site/{GID1}_vs_{GID2}.site_matrix.gz",
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
        --binSize 100  \
        --downstream 5000 --upstream 5000 \
        -o {output.mat}
        """


rule site_heatmapdiff:
    input:
        mat="results/matdiff_site/{GID1}_vs_{GID2}.site_matrix.gz",
    output:
        plot="results/heatdiff/{GID1}_vs_{GID2}.site.pdf",
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
        --colorMap 'RdBu_r' \
        --zMin -1.5 --zMax 1.5 \
        --xAxisLabel 'distance' \
        --sortRegions 'keep' \
        -o {output.plot}
        """
