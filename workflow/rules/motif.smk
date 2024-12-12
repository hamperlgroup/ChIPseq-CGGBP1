rule homer:
    input:
        peaks = "results/ChIP-peaks/narrowPeak/merged/{IPGID}.bed"
    output:
        "results/motifs/narrowPeak/homer/{IPGID}/homerResults.html"
    params:
        genomeVersion = "hg38",
        outDir = "results/motifs/narrowPeak/homer/{IPGID}/",
        cores = 8
    conda:
        'motif'
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        findMotifsGenome.pl {input.peaks} {params.genomeVersion} {params.outDir} -size 200 -mask -p {params.cores}
        """

# bedtools intersect -a results/ChIP-peaks/narrowPeak/merged/ChIP-seq_CGGBP1_ENCODE_noCHR.bed -b data/genome/protein_coding_genes.bed -wa | uniq > results/ChIP-peaks/narrowPeak/merged/ChIP-seq_CGGBP1_ENCODE_proteinCodingGeneRegions.bed
# tail -n +2 results/ChIP-peaks/narrowPeak/merged/peak_annotation/annotation_ChIP-seq_CGGBP1_ENCODE.tsv | cut -f 1,2,3,7 | grep 'Promoter' | sed 's/chr//g' > results/ChIP-peaks/narrowPeak/merged/ChIP-seq_CGGBP1_ENCODE_promoterRegions.bed
# scale the FASTA sequences to a unique size e.g. 100 bp for a transcription factor or 500 bp for histone marks
rule meme_input:
    input:
        # peaks = "results/ChIP-peaks/narrowPeak/merged/{IPGID}.bed"
        # peaks = "results/ChIP-peaks/narrowPeak/merged/{IPGID}_proteinCodingGeneRegions.bed"
        peaks = "results/ChIP-peaks/narrowPeak/merged/{IPGID}_promoterRegions.bed"
    output:
        fasta = "results/motifs/narrowPeak/meme/fasta/{IPGID}.fa",
        shuffled = "results/motifs/narrowPeak/meme/fasta/{IPGID}_shuffled.fa"
    params:
        outDir = "results/motifs/narrowPeak/meme/fasta",
        sample = "{IPGID}",
        genome = "data/genome/genome.fa",
        resizeFasta = 50
    conda:
        'motif'
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        mkdir -p {params.outDir}
        
        awk -F '\\t' '{{X={params.resizeFasta}; mid=(int($2)+int($3))/2;printf("%s\\t%d\\t%d\\n",$1,(mid-X<0?0:mid-X),mid+X);}}' {input.peaks} | sed 's/chr//g' > {params.outDir}/{params.sample}_{params.resizeFasta}-{params.resizeFasta}bp.bed

        bedtools getfasta -fi {params.genome} -bed {params.outDir}/{params.sample}_{params.resizeFasta}-{params.resizeFasta}bp.bed -fo {output.fasta}

        fasta-shuffle-letters -k 2 -dna {output.fasta} {output.shuffled}
        """


rule meme_backgroud:
    input:
        fasta = "results/motifs/narrowPeak/meme/fasta/{IPGID}.fa",
        shuffled = "results/motifs/narrowPeak/meme/fasta/{IPGID}_shuffled.fa"
    output:
        "results/motifs/narrowPeak/meme/{IPGID}_Overshuffle/meme-chip.html"
    params:
        outDir = "results/motifs/narrowPeak/meme/{IPGID}_Overshuffle",
        motifDB = "data/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"
    conda:
        'motif'
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        meme-chip -oc {params.outDir} -db {params.motifDB} -neg {input.shuffled} -dna -meme-nmotifs 25 -minw 8 -maxw 12 -meme-mod zoops {input.fasta}
        """

rule meme_no_backgroud:
    input:
        fasta = "results/motifs/narrowPeak/meme/fasta/{IPGID}.fa"
    output:
        "results/motifs/narrowPeak/meme/{IPGID}_noBck/meme-chip.html"
    params:
        outDir = "results/motifs/narrowPeak/meme/{IPGID}_noBck",
        motifDB = "data/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"
    conda:
        'motif'
    resources:
        mem_mb = 24000
    threads: 24
    shell:
        """
        meme-chip -oc {params.outDir} -db {params.motifDB} -dna -meme-nmotifs 25 -minw 8 -maxw 12 -meme-mod zoops {input.fasta}
        """
