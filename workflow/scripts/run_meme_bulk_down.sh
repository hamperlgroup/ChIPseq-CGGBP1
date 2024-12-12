#!/bin/bash

#SBATCH -J bulk_down_MEME
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=24G
#SBATCH -t 03:00:00
#SBATCH --output=/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/logs/meme_bound_genes_manual/submit_bulk_down.out
#SBATCH --error=/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/logs/meme_bound_genes_manual/submit_bulk_down.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elizabeth.marquezgomez@helmholtz-munich.de

start=`date +%s`
##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

echo "------------------------------------------------------"
echo -n "Job is running on node "; $SLURM_JOB_NODELIST
echo "------------------------------------------------------"
echo SLURM: qsub is running on $SLURM_SUBMIT_HOST
echo SLURM: originating queue is $SLURM_JOB_PARTITION
echo SLURM: executing queue is $SLURM_JOB_PARTITION
echo SLURM: working directory is $SLURM_SUBMIT_DIR
echo SLURM: job identifier is $SLURM_JOB_ID
echo SLURM: job name is $SLURM_JOB_NAME
echo SLURM: node file is $SLURM_JOB_NODELIS
echo "------------------------------------------------------"

# The working directory for the job is inside the scratch directory

echo "------------------------------------------------------"
echo -n "Job is running on node "; $SLURM_JOB_NODELIST
echo "------------------------------------------------------"
echo " "
echo " "

###############################################################
#                                                             #
#                      Main Program                           #
#                                                             #
###############################################################

mamba activate motif
cd /lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1

## Scale bed files
awk -F '\t' '{X=50; mid=(int($2)+int($3))/2;printf("%s\t%d\t%d\n",$1,(mid-X<0?0:mid-X),mid+X);}' results/ChIP-peaks/narrowPeak/SLAM-seq_overlap/merged/overlap_ChIP-seq_CGGBP1_ENCODE-SLAM_bulk_down_genes.bed | sed 's/chr//g' > results/motifs/narrowPeak/meme/fasta/bound_genes_bulk_down_50-50bp.bed

## Convert bed to fasta
bedtools getfasta -fi data/genome/genome.fa -bed results/motifs/narrowPeak/meme/fasta/bound_genes_bulk_down_50-50bp.bed -fo results/motifs/narrowPeak/meme/fasta/bound_genes_bulk_down.fa

## Generate shuffled fasta
fasta-shuffle-letters -k 2 -dna results/motifs/narrowPeak/meme/fasta/bound_genes_bulk_down.fa results/motifs/narrowPeak/meme/fasta/bound_genes_bulk_down_shuffled.fa

## Run MEME
meme-chip -oc results/motifs/narrowPeak/meme/bound_genes_bulk_down_Overshuffle -db data/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme -neg results/motifs/narrowPeak/meme/fasta/bound_genes_bulk_down_shuffled.fa -dna -meme-nmotifs 25 -minw 8 -maxw 12 -meme-mod zoops results/motifs/narrowPeak/meme/fasta/bound_genes_bulk_down.fa

