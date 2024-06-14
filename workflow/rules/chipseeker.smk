rule peak_calling:
    input:
        bam_sample = "results/BAM/{IP}.rmdup.bam",
        bam_ctrl = 
        lambda wildcards: expand("results/BAM/{ctrl}.rmdup.bam", ctrl = SampleTable['control'][SampleTable['fastq_1'] == f'{wildcards.IP}_1.fastq.gz'])
    output:
        peaks = "results/ChIP-peaks/peak_calling/{IP}_peaks_chr.narrowPeak"
    params:
        outdir = "results/ChIP-peaks/peak_calling",
        sample = "{IP}",
        genomeSize = "hs",
        narrowPeak = "results/ChIP-peaks/peak_calling/{IP}_peaks.narrowPeak"
    log:
        "results/logs/MACS3/{IP}-macs3.log"
    conda:
        "macs3"
    resources:
        mem_mb = 50000
    threads: 24
    shell:
        """
            macs3 callpeak -t {input.bam_sample} \
                -c {input.bam_ctrl} \
                -f BAM -g {params.genomeSize} \
                -n {params.sample} \
                --outdir {params.outdir} 2> {log}

            awk 'OFS="\\t" {{$1="chr"$1; print}}' {params.narrowPeak} > {output.peaks}
        """

rule peak_annotation:
    input:
        peaks = "results/ChIP-peaks/peak_calling/{IP}_peaks_chr.narrowPeak"
    output:
        regions = "results/ChIP-peaks/peak_annotation/peak_regions_{IP}.bed"
    params:
        sample = "{IP}"
    conda:
        # "../envs/env_chipSeeker.yaml"
        "overlap_OK-DRIPc"
    resources:
        mem_mb = 50000
    threads: 24
    shell:
        """
            Rscript workflow/scripts/ChIPseeker.R --sampleName {params.sample} --inputFile {input.peaks}

            sed -i 's/"//g' {output.regions}
        """