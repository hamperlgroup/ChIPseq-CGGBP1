#!/bin/bash

#SBATCH -J Snake_ChIP-seq_pub
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=8G
#SBATCH -t 03-00
#SBATCH --output=/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/logs/submit_snake.out
#SBATCH --error=/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/logs/submit_snake.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elizabeth.marquezgomez@helmholtz-munich.de

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

echo "------------------------------------------------------"
echo "------------------------------------------------------"
echo SLURM: qsub is running on $SLURM_SUBMIT_HOST
echo SLURM: originating queue is $SLURM_JOB_PARTITION
echo SLURM: executing queue is $SLURM_JOB_PARTITION
echo SLURM: working directory is $SLURM_SUBMIT_DIR
echo SLURM: job identifier is $SLURM_JOB_ID
echo "------------------------------------------------------"

# The working directory for the job is inside the scratch directory

echo "------------------------------------------------------"
echo "------------------------------------------------------"
echo " "
echo " "

###############################################################
#                                                             #
#                      Main Program                           #
#                                                             #
###############################################################

# snakemake --profile slurm_profile --snakefile workflow/Snakefile --configfile config/config.yaml --dry-run --unlock
# snakemake --profile slurm_profile --snakefile workflow/Snakefile --configfile config/config.yaml --dry-run > snakemake_log.txt

# snakemake --dag | dot -Tpdf > dag_all.pdf
# snakemake --rulegraph | dot -Tpdf > rule_graph.pdf

WORKDIR=/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1
cd $WORKDIR
log_path=${WORKDIR}/logs/

mamba activate env_chipseq_pipeline

start=`date +%s`
date
echo 'Running snakemake'

# --rerun-triggers mtime
snakemake --profile slurm_profile --snakefile workflow/Snakefile --configfile config/config.yaml --stats $log_path/snakemake.stats >& $log_path/snakemake.log  --rerun-incomplete --use-conda 


end=`date +%s`
runtime=$((end-start))
echo 'Running time ' ${runtime} ' seconds'