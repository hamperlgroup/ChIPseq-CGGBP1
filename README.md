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
