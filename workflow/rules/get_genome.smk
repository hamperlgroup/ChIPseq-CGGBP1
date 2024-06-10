rule get_genome_files:
    params:
        genomeLocation=config['genomeFTP'],
        gtfLocation=config['gtfFTP'],
    output:
        fasta1="data/genome/genome.fa",
        gtf1=  "data/genome/genome.gtf",
        fai1=  "data/genome/genome.fa.fai",
        fasta2="data/genome/combined.fa",
        gtf2=  "data/genome/combined.gtf",
        fai2=  "data/genome/combined.fa.fai",
    conda:
        "../envs/env_preprocessing.yaml",
    resources:
        mem_mb = 8000
    threads: 2
    shell:
        """
        wget -O {output.fasta1}.gz {params.genomeLocation}
        gunzip {output.fasta1}.gz
        samtools faidx {output.fasta1}

        wget -O {output.gtf1}.gz {params.gtfLocation}
        gunzip {output.gtf1}.gz

        # cat {output.fasta1} data/external/reporters.fa > {output.fasta2}
        cp {output.fasta1} {output.fasta2}
        samtools faidx {output.fasta2}

        # cat {output.gtf1} data/external/reporters.gtf > {output.gtf2}
        cp {output.gtf1} {output.gtf2}
        """
