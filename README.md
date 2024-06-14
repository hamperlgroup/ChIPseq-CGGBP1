# Pipeline developed by Tamas Schauer for project 202205_chipseq_polii_ser2p

/lustre/groups/ies/projects/metp_lab/tamas.schauer/Projects/Marcel/202205_chipseq_polii_ser2p

# Environment

```
mamba env create --name env_chipseq_snakemake --file workflow/envs/env_snakemake.yaml
```

# Workflow

```
snakemake --rulegraph | dot -Tpng > images/snakemake_rulegraph.png
```

![](images/snakemake_rulegraph.png)

# Run

```
conda activate env_chipseq_snakemake

export TMPDIR="/localscratch"

sbatch run.sbatch
```
