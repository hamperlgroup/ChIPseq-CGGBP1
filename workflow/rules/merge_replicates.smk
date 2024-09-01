rule merging_replicates_bedtools:
    input:
        peaksFiles = lambda wildcards: expand("results/ChIP-peaks/narrowPeak/peak_calling/{IPG}_peaks.narrowPeak", IPG = SampleTable['sample_ID'][SampleTable['group'] == wildcards.IPGID])
    output:
        "results/ChIP-peaks/narrowPeak/merged/{IPGID}.bed"
    params:
        peaks = lambda wildcards: " ".join(expand("results/ChIP-peaks/narrowPeak/peak_calling/{IPG}_peaks.narrowPeak", IPG = SampleTable['sample_ID'][SampleTable['group'] == wildcards.IPGID]))
    conda:
        'bedtools'
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        cat {params.peaks} | awk 'OFS="\\t" {{$1="chr"$1; print}}' | sortBed | bedtools merge -c 1 -o count | awk '$4>=2' > {output}
        """

rule peak_annot_merged:
    input:
        peaks = "results/ChIP-peaks/narrowPeak/merged/{IPGID}.bed"
    output:
        "results/ChIP-peaks/narrowPeak/merged/peak_annotation/annotation_{IPGID}.tsv"
    params:
        sample = "{IPGID}",
        peakMode = "narrowPeak"
    conda:
        # "../envs/env_chipSeeker.yaml" <-- DONT FORGET TO UPDATE THE ENVIRONMENT BASED ON THE LOCAL ONE
        "overlap_OK-DRIPc"
    resources:
        mem_mb = 50000
    threads: 24
    shell:
        """
            Rscript workflow/scripts/peak_annot_merged.R --sampleName {params.sample} --inputFile {input.peaks} --peakMode {params.peakMode}
        """
