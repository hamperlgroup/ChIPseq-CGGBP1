rule merging_replicates:
    input:
        peaksFiles = lambda wildcards: expand("results/ChIP-peaks/narrowPeak/peak_calling/{ID}_peaks.narrowPeak", ID = SampleTable['id'][SampleTable['group'] == wildcards.IPGID])
    output:
        "results/ChIP-peaks/narrowPeak/merged/merged_{IPGID}.bed"
    params:
        groupName = "{IPGID}",
        peaks = lambda wildcards: ",".join(expand("results/ChIP-peaks/narrowPeak/peak_calling/{ID}_peaks.narrowPeak", ID = SampleTable['id'][SampleTable['group'] == wildcards.IPGID])),
    conda:
        'overlap_OK-DRIPc'
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        Rscript workflow/scripts/merge_peaks.R --groupName {params.groupName} --peaksFiles {params.peaks}
        """


rule CGGBP1_overlap_all_genes:
    input:
        mergedCGGBP1 = "results/ChIP-peaks/narrowPeak/merged/merged_{CGGBP1}.bed",
        genes = "data/genome/genes_genome.gtf"
    output:
        "results/overlap_CGGBP1_all_genes/overlapping_CGGBP1_all_genes_{CGGBP1}.bed",
        "results/overlap_CGGBP1_all_genes/nonOverlapping_CGGBP1_all_genes_{CGGBP1}.bed"
    params:
        sampleName = "{CGGBP1}"
    conda:
        'overlap_OK-DRIPc'
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        Rscript workflow/scripts/CGGBP1_overlap_allGenes.R --sampleName {params.sampleName} --mergedCGGBP1 {input.mergedCGGBP1} --genesFile {input.genes}
        """


# rule CGGBP1allGenes_scaled_matrix_group:
#     input:
#         bw="/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1_insertion/results/covmean/{PS2_5}.coverage.bw",
#         genesOverlap=expand("results/overlap_CGGBP1_all_genes/overlapping_CGGBP1_all_genes_{CGGBP1}.bed", CGGBP1 = ['ChIP-seq_CGGBP1_ENCODE']),
#         genesNotOverlap=expand("results/overlap_CGGBP1_all_genes/notOverlapping_CGGBP1_all_genes_{CGGBP1}.bed", CGGBP1 = ['ChIP-seq_CGGBP1_ENCODE'])
#     output:
#         mat="results/overlap_CGGBP1_all_genes/matrix_scaled/{PS2_5}.scaled_matrix.gz",
#     conda:
#         "../envs/env_preprocessing.yaml",
#     resources:
#         mem_mb = 24000
#     threads: 24
#     shell:
#         """
#         computeMatrix scale-regions \
#         --numberOfProcessors {threads} \
#         -S {input.bw} -R {input.genesOverlap} {input.genesNotOverlap} \
#         --binSize 20 --skipZeros \
#         --downstream 2000 --upstream 2000 \
#         --regionBodyLength 4000 \
#         --missingDataAsZero \
#         -o {output.mat}
#         """

# rule CGGBP1allGenes_scaled_plot_group:
#     input:
#         mat="results/overlap_CGGBP1_all_genes/matrix_scaled/{PS2_5}.scaled_matrix.gz",
#     output:
#         plot="results/overlap_CGGBP1_all_genes/metaplots/{PS2_5}.scaled_plot.png",
#     conda:
#         "../envs/env_preprocessing.yaml",
#     resources:
#         mem_mb = 8000
#     threads: 4
#     shell:
#         """
#         plotProfile -m {input.mat} \
#         --plotHeight 10 --plotWidth 14 \
#         --yMin 0.05 --yMax 1.0 \
#         -o {output.plot}
#         """

# rule CGGBP1allGenes_heatmap_TSS_group:
#     input:
#         mat="results/overlap_CGGBP1_all_genes/matrix_scaled/{PS2_5}.scaled_matrix.gz",
#     output:
#         plot="results/overlap_CGGBP1_all_genes/metaplots/{PS2_5}.heatmap_TSS.png",
#     conda:
#         "../envs/env_preprocessing.yaml",
#     resources:
#         mem_mb = 8000
#     threads: 4
#     shell:
#         """
#         plotHeatmap -m {input.mat} \
#         --colorMap magma --heatmapHeight 50 --heatmapWidth 10 \
#         -o {output.plot}
#         """