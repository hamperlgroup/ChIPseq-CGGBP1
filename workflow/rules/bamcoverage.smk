rule fragsize:
    input:
        rmdup="results/BAM/rmdup/{ID}.rmdup.bam",
    output:
        png="results/fragsize/{ID}.fragmentsize.png",
    conda:
        "../envs/env_preprocessing.yaml",
    log:
        "results/logs/fragsize/{ID}.log",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        bamPEFragmentSize -hist {output.png} -p {threads} \
        -T "Fragment size" --maxFragmentLength 800 \
        --samplesLabel {wildcards.ID} \
        -b {input.rmdup} > {log}
        """

rule bamcoverage2:
    input:
        rmdup=lambda wildcards: expand("results/BAM/rmdup/{ID}.rmdup.bam", ID = SampleTable['sample_ID'][SampleTable['nid'] == wildcards.NID]),
        bed="data/genome/hg38-blacklist.v2.bed",
    output:
        bw="results/cov2/{NID}.coverage.bw",
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        bamCoverage --bam {input.rmdup} -o {output.bw} \
        --blackListFileName {input.bed} \
        --ignoreForNormalization MT \
        --binSize 20 --smoothLength 60 \
        --extendReads  \
        --maxFragmentLength 700 \
        --normalizeUsing "CPM" \
        --numberOfProcessors {threads}
        """

rule bamcoverage:
    input:
        rmdup="results/BAM/rmdup/{ID}.rmdup.bam",
        bed="data/genome/hg38-blacklist.v2.bed",
    output:
        bw="results/coverage/{ID}.coverage.bw",
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        bamCoverage --bam {input.rmdup} -o {output.bw} \
        --blackListFileName {input.bed} \
        --ignoreForNormalization MT \
        --binSize 20 --smoothLength 60 \
        --extendReads  \
        --maxFragmentLength 700 \
        --normalizeUsing "CPM" \
        --numberOfProcessors {threads}
        """


rule make_windows:
    input:
        fai="data/genome/genome.fa.fai",
    output:
        bed="data/genome/bins.bed",
    params:
        binsize = 1000
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 8000
    threads: 4
    shell:
        """
        bedtools makewindows -g {input.fai} -w {params.binsize} | grep -v "KI\|GL\|MT" | sort -k1,1 -k2,2n > {output.bed}
        """


rule multisummary:
    input:
        bam=expand("results/BAM/rmdup/{ID}.rmdup.bam", ID = sample_id),
        bins="data/genome/bins.bed",
        bed="data/genome/hg38-blacklist.v2.bed",
    output:
        npz="results/counts/counts_per_bin.npz",
        tab="results/counts/counts_per_bin.tab",
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        multiBamSummary BED-file --BED {input.bins}  \
        --numberOfProcessors {threads} \
        --bamfiles {input.bam} \
        --blackListFileName {input.bed} \
        --smartLabels \
        --extendReads \
        --centerReads \
        --samFlagInclude 64 \
        --outFileName  {output.npz} \
        --outRawCounts {output.tab}
        """


rule meancov_bw:
    input:
        bw= lambda wildcards: expand("results/coverage/{ID}.coverage.bw", ID = SampleTable['sample_ID'][SampleTable['group'] == wildcards.GID]),
    output:
        bw="results/covmean/{GID}.coverage.bw",
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        bigwigAverage --bigwigs {input.bw} -o {output.bw} \
        --binSize 20 --numberOfProcessors {threads}
        """

rule meancov_bg:
    input:
        bw= lambda wildcards: expand("results/coverage/{ID}.coverage.bw", ID = SampleTable['sample_ID'][SampleTable['group'] == wildcards.GID]),
    output:
        bg="results/covmean/{GID}.coverage.bedgraph",
    conda:
        "../envs/env_preprocessing.yaml",
    params:
        tmp="results/covmean/tmp_{GID}.coverage.bedgraph",
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        bigwigAverage --bigwigs {input.bw} -o {params.tmp} \
        --binSize 20 --numberOfProcessors {threads} \
        --outFileFormat bedgraph

        awk 'OFS="\\t" {{$1="chr"$1; print}}' {params.tmp} > {output.bg}
        """


# rule diffmeancov:
#     input:
#         bw1="results/covmean/{GID1}.coverage.bw",
#         bw2="results/covmean/{GID2}.coverage.bw",
#     output:
#         diff="results/covdiff/{GID1}_vs_{GID2}.coverage.bw",
#     conda:
#         "../envs/env_preprocessing.yaml",
#     resources:
#         mem_mb = 24000
#     threads: 24
#     shell:
#         """
#         bigwigCompare --bigwig1 {input.bw1} --bigwig2 {input.bw2} -o {output.diff} \
#         --binSize 100 --fixedStep --numberOfProcessors {threads}
#         """
