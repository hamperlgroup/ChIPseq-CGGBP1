rule get_genome_files:
    params:
        genomeLocation=config['genomeFTP'],
        gtfLocation=config['gtfFTP'],
    output:
        fasta1="data/genome/genome.fa",
        gtf1=  "data/genome/genome.gtf",
        fai1=  "data/genome/genome.fa.fai",
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

        samtools faidx {output.fasta1}
        """
