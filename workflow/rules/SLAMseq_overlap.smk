rule SLAMseq_overlap:
    input:
        peaks = "results/ChIP-peaks/{PEAKS}/peak_annotation/annotation_{IP}.tsv"
    output:
        "results/ChIP-peaks/{PEAKS}/SLAM-seq_overlap/replicates/overlap_{IP}-SLAM_bulk_unchanged.tsv",
        "results/ChIP-peaks/{PEAKS}/SLAM-seq_overlap/replicates/overlap_{IP}-SLAM_nascent_unchanged.tsv"
    params:
        sample = "{IP}",
        fc = 1,
        pval = 1,
        peakMode = "{PEAKS}",
        sampleMode = "replicates"
    conda:
        # "../envs/env_chipSeeker.yaml" <-- DONT FORGET TO UPDATE THE ENVIRONMENT BASED ON THE LOCAL ONE
        "overlap_OK-DRIPc"
    resources:
        mem_mb = 50000
    threads: 24
    shell:
        """
            Rscript workflow/scripts/slam-seq_gene_overlap.R --sampleName {params.sample} --inputFile {input.peaks} --absCutoffFC {params.fc} --cutoffPval {params.pval} --peakMode {params.peakMode} --sampleMode {params.sampleMode}
        """

rule SLAMseq_consensus:
    input:
        bulk = lambda wildcards: expand("results/ChIP-peaks/{PEAKS}/SLAM-seq_overlap/replicates/overlap_{IPG}-SLAM_bulk_unchanged.tsv", IPG = SampleTable['sample_ID'][SampleTable['group'] == wildcards.IPGID], PEAKS = wildcards.PEAKS),
        nascent = lambda wildcards: expand("results/ChIP-peaks/{PEAKS}/SLAM-seq_overlap/replicates/overlap_{IPG}-SLAM_nascent_unchanged.tsv", IPG = SampleTable['sample_ID'][SampleTable['group'] == wildcards.IPGID], PEAKS = wildcards.PEAKS)
    output:
        bulk_group = "results/ChIP-peaks/{PEAKS}/SLAM-seq_overlap/replicates/grouped_overlap_{IPGID}-SLAM_bulk_unchanged.tsv",
        nascent_group = "results/ChIP-peaks/{PEAKS}/SLAM-seq_overlap/replicates/grouped_overlap_{IPGID}-SLAM_nascent_unchanged.tsv"
    conda:
        "bioinfo"
    params:
        group = "{IPGID}",
        peaks = "{PEAKS}",
        path = "results/ChIP-peaks/{PEAKS}/SLAM-seq_overlap/replicates",
        files_bulk = lambda wildcards: " ".join(expand("results/ChIP-peaks/{PEAKS}/SLAM-seq_overlap/replicates/overlap_{IPG}-SLAM_bulk_unchanged.tsv", IPG = SampleTable['sample_ID'][SampleTable['group'] == wildcards.IPGID], PEAKS = wildcards.PEAKS)),
        files_nascent = lambda wildcards: " ".join(expand("results/ChIP-peaks/{PEAKS}/SLAM-seq_overlap/replicates/overlap_{IPG}-SLAM_nascent_unchanged.tsv", IPG = SampleTable['sample_ID'][SampleTable['group'] == wildcards.IPGID], PEAKS = wildcards.PEAKS))
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        for file_b in {params.files_bulk}; do cut -f1 $file_b >> ${params.path}/temp_{params.group}_b.tsv; done;
        sort ${params.path}/temp_{params.group}_b.tsv | uniq > {output.bulk_group}
        rm ${params.path}/temp_{params.group}_b.tsv

        for file_n in {params.files_nascent}; do cut -f1 $file_n >> ${params.path}/temp_{params.group}_n.tsv; done;
        sort ${params.path}/temp_{params.group}_n.tsv | uniq > {output.nascent_group}
        rm ${params.path}/temp_{params.group}_n.tsv
        """

## In merged replicates
rule SLAMseq_overlap_merged:
    input:
        peaks = "results/ChIP-peaks/{PEAKS}/merged/peak_annotation/annotation_{IPGID}.tsv"
    output:
        "results/ChIP-peaks/{PEAKS}/SLAM-seq_overlap/merged/overlap_{IPGID}-SLAM_bulk_unchanged.tsv",
        "results/ChIP-peaks/{PEAKS}/SLAM-seq_overlap/merged/overlap_{IPGID}-SLAM_nascent_unchanged.tsv"
    params:
        sample = "{IPGID}",
        fc = 1,
        pval = 1,
        peakMode = "{PEAKS}",
        sampleMode = "merged"
    conda:
        # "../envs/env_chipSeeker.yaml" <-- DONT FORGET TO UPDATE THE ENVIRONMENT BASED ON THE LOCAL ONE
        "overlap_OK-DRIPc"
    resources:
        mem_mb = 50000
    threads: 24
    shell:
        """
            Rscript workflow/scripts/slam-seq_gene_overlap.R --sampleName {params.sample} --inputFile {input.peaks} --absCutoffFC {params.fc} --cutoffPval {params.pval} --peakMode {params.peakMode} --sampleMode {params.sampleMode}
        """