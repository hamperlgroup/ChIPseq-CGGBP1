rule peak_calling_narrow:
    input:
        bam_sample = "results/BAM/rmdup/{IP}.rmdup.bam",
        bam_ctrl = 
        lambda wildcards: expand("results/BAM/rmdup/{ctrl}.rmdup.bam", ctrl = SampleTable['control'][SampleTable['fastq_1'] == f'{wildcards.IP}_1.fastq.gz'])
    output:
        peaks = "results/ChIP-peaks/narrowPeak/peak_calling/{IP}_peaks.narrowPeak"
    params:
        outdir = "results/ChIP-peaks/narrowPeak/peak_calling",
        sample = "{IP}",
        genomeSize = "hs"
    log:
        "results/logs/MACS3/{IP}-narrowPeak_macs3.log"
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
        """


# rule peak_calling_broad:
#     input:
#         bam_sample = "results/BAM/{IP}.rmdup.bam",
#         bam_ctrl = 
#         lambda wildcards: expand("results/BAM/{ctrl}.rmdup.bam", ctrl = SampleTable['control'][SampleTable['fastq_1'] == f'{wildcards.IP}_1.fastq.gz'])
#     output:
#         peaks = "results/ChIP-peaks/broadPeak/peak_calling/{IP}_peaks.broadPeak"
#     params:
#         outdir = "results/ChIP-peaks/broadPeak/peak_calling",
#         sample = "{IP}",
#         genomeSize = "hs",
#     log:
#         "results/logs/MACS3/{IP}-broadPeak_macs3.log"
#     conda:
#         "macs3"
#     resources:
#         mem_mb = 50000
#     threads: 24
#     shell:
#         """
#             macs3 callpeak -t {input.bam_sample} \
#                 -c {input.bam_ctrl} \
#                 -f BAM -g {params.genomeSize} \
#                 -n {params.sample} \
#                 --broad \
#                 --outdir {params.outdir} 2> {log}
#         """


rule peak_annotation:
    input:
        peaks = "results/ChIP-peaks/{PEAKS}/peak_calling/{IP}_peaks.{PEAKS}"
    output:
        "results/ChIP-peaks/{PEAKS}/peak_annotation/annotation_{IP}.tsv"
    params:
        sample = "{IP}",
        peakMode = "{PEAKS}"
    conda:
        # "../envs/env_chipSeeker.yaml" <-- DONT FORGET TO UPDATE THE ENVIRONMENT BASED ON THE LOCAL ONE
        "overlap_OK-DRIPc"
    resources:
        mem_mb = 50000
    threads: 24
    shell:
        """
            Rscript workflow/scripts/ChIPseeker.R --sampleName {params.sample} --inputFile {input.peaks} --peakMode {params.peakMode}
        """

# onsuccess:
#     print("Workflow finished, no error")
# rule MultiQC:
#     input:
#         peaks = "results/ChIP-peaks/{PEAKS}/peak_annotation/annotation_{IP}.tsv"
#     output:
#         "results/MultiQC/{PEAKS}/multiqc_report.html"
#     params:
#         outdir = "results/MultiQC/{PEAKS}",
#         peakMode = "{PEAKS}"
#     conda:
#         "multiQC"
#     resources:
#         mem_mb = 50000
#     threads: 24
#     shell:
#         """
#             multiqc results/fastqc \
#             results/logs/bowtie2 \
#             results/logs/fragsize \
#             results/logs/MACS3/*{params.peakMode}* \
#             results/logs/picard \
#             results/trimmed \
#             results/ChIP-peaks/{params.peakMode}/peak_calling \
#             --outdir {params.outdir}
#         """


# multiqc results/fastqc results/logs/bowtie2 results/logs/fragsize results/logs/MACS3/*broadPeak* results/logs/picard results/trimmed results/ChIP-peaks/broadPeak/peak_calling --outdir results/MultiQC/broadPeak

# multiqc results/fastqc results/logs/bowtie2 results/logs/fragsize results/logs/MACS3/*narrowPeak* results/logs/picard results/trimmed results/ChIP-peaks/narrowPeak/peak_calling --outdir results/MultiQC/narrowPeak