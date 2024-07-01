rule coding_genes_overlap:
    input:
        peaks = "results/ChIP-peaks/narrowPeak/peak_calling/{IP}_peaks.narrowPeak",
        genes = "data/genome/protein_coding_genes.bed"
    output:
        "results/overlap_coding_genes/overlapping_coding_genes_{IP}.bed",
        "results/overlap_coding_genes/notOverlapping_coding_genes_{IP}.bed",
        "results/overlap_coding_genes/overlappingRegions_{IP}.pdf"
    params:
        sampleName = "{IP}"
    conda:
        'overlap_OK-DRIPc'
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        Rscript workflow/scripts/protein_coding_genes_overlap.R --sampleName {params.sampleName} --peaksFile {input.peaks} --genesFile {input.genes}
        """


rule codingGenes_scaled_matrix:
    input:
        bw="results/coverage/{IP}.coverage.bw",
        genesOverlap="results/overlap_coding_genes/overlapping_coding_genes_{IP}.bed",
        genesNotOverlap="results/overlap_coding_genes/notOverlapping_coding_genes_{IP}.bed"
    output:
        mat="results/overlap_coding_genes/matrix_scaled/replicates/{IP}.scaled_matrix.gz"
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        computeMatrix scale-regions \
        --numberOfProcessors {threads} \
        -S {input.bw} -R {input.genesOverlap} {input.genesNotOverlap} \
        --binSize 20 --skipZeros \
        --downstream 2000 --upstream 2000 \
        --regionBodyLength 4000 \
        --missingDataAsZero \
        -o {output.mat}
        """

rule codingGenes_scaled_plot:
    input:
        mat="results/overlap_coding_genes/matrix_scaled/replicates/{IP}.scaled_matrix.gz"
    output:
        plot="results/overlap_coding_genes/metaplots/replicates/{IP}.scaled_plot.png",
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

rule codingGenes_heatmap_TSS:
    input:
        mat="results/overlap_coding_genes/matrix_scaled/replicates/{IP}.scaled_matrix.gz"
    output:
        plot="results/overlap_coding_genes/metaplots/replicates/{IP}.heatmap_TSS.png",
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 8000
    threads: 4
    shell:
        """
        plotHeatmap -m {input.mat} \
        --colorMap magma --heatmapHeight 50 --heatmapWidth 10 \
        -o {output.plot}
        """