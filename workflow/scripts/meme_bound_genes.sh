#!/bin/bash

############    -  Motif MEME CGGBP1 bound genes  -    ############

# ---
# Author: Elizabeth Marquez-Gomez
# Date: 2024-12
# Version: 1
# Subversion: 0
# Environment: motif
# ---

### -------------------------- Description -------------------------- ###
# Overlap


### RUNS ##
## sh /lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/workflow/scripts/meme_bound_genes.sh --rnas bulk,nascent --states up,down --fastaLength 50 --execute 0

####################################### ----------------------------------- #######################################
if [ "$#" == "0" ]; then
echo "Error: No arguments" >&2; exit 1;
fi

while (( "$#" )); do
  case "$1" in
    --rnas)
      rnas=$2
      shift 2
      ;;
    --states)
      states=$2
      shift 2
      ;;
    --fastaLength)
      length=$2
      shift 2
      ;;
    --execute)
      execute=$2
      shift 2
      ;;
    -*|--*=|*) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
  esac
done


root_path=/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1

echo "------------------------------------------------------"
echo -n "MEME Analysis"
echo "------------------------------------------------------"
echo "Running RNAs: ${rnas}"
echo "States: ${states}"
echo "Execution mode: ${execute}"
echo "------------------------------------------------------"


for rna in $(echo $rnas | sed "s/,/ /g"); do
echo "Working on: ${rna}"

for state in $(echo $states | sed "s/,/ /g"); do
echo "Working on: ${state}"

file=${root_path}/workflow/scripts/run_meme_${rna}_${state}.sh
if test -f "${file}"; then
rm "${file}"
fi

## MEME array script --------------------------
cat <<EOT >> "${file}"
#!/bin/bash

#SBATCH -J ${rna}_${state}_MEME
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=24G
#SBATCH -t 03:00:00
#SBATCH --output=${root_path}/logs/meme_bound_genes_manual/submit_${rna}_${state}.out
#SBATCH --error=${root_path}/logs/meme_bound_genes_manual/submit_${rna}_${state}.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elizabeth.marquezgomez@helmholtz-munich.de

start=\`date +%s\`
##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

echo "------------------------------------------------------"
echo -n "Job is running on node "; \$SLURM_JOB_NODELIST
echo "------------------------------------------------------"
echo SLURM: qsub is running on \$SLURM_SUBMIT_HOST
echo SLURM: originating queue is \$SLURM_JOB_PARTITION
echo SLURM: executing queue is \$SLURM_JOB_PARTITION
echo SLURM: working directory is \$SLURM_SUBMIT_DIR
echo SLURM: job identifier is \$SLURM_JOB_ID
echo SLURM: job name is \$SLURM_JOB_NAME
echo SLURM: node file is \$SLURM_JOB_NODELIS
echo "------------------------------------------------------"

# The working directory for the job is inside the scratch directory

echo "------------------------------------------------------"
echo -n "Job is running on node "; \$SLURM_JOB_NODELIST
echo "------------------------------------------------------"
echo " "
echo " "

###############################################################
#                                                             #
#                      Main Program                           #
#                                                             #
###############################################################

mamba activate motif
cd ${root_path}

## Scale bed files
awk -F '\t' '{X=${length}; mid=(int(\$2)+int(\$3))/2;printf("%s\t%d\t%d\n",\$1,(mid-X<0?0:mid-X),mid+X);}' results/ChIP-peaks/narrowPeak/SLAM-seq_overlap/merged/overlap_ChIP-seq_CGGBP1_ENCODE-SLAM_${rna}_${state}_genes.bed | sed 's/chr//g' > results/motifs/narrowPeak/meme/fasta/bound_genes_${rna}_${state}_${length}-${length}bp.bed

## Convert bed to fasta
bedtools getfasta -fi data/genome/genome.fa -bed results/motifs/narrowPeak/meme/fasta/bound_genes_${rna}_${state}_${length}-${length}bp.bed -fo results/motifs/narrowPeak/meme/fasta/bound_genes_${rna}_${state}.fa

## Generate shuffled fasta
fasta-shuffle-letters -k 2 -dna results/motifs/narrowPeak/meme/fasta/bound_genes_${rna}_${state}.fa results/motifs/narrowPeak/meme/fasta/bound_genes_${rna}_${state}_shuffled.fa

## Run MEME
meme-chip -oc results/motifs/narrowPeak/meme/bound_genes_${rna}_${state}_Overshuffle -db data/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme -neg results/motifs/narrowPeak/meme/fasta/bound_genes_${rna}_${state}_shuffled.fa -dna -meme-nmotifs 25 -minw 8 -maxw 12 -meme-mod zoops results/motifs/narrowPeak/meme/fasta/bound_genes_${rna}_${state}.fa

EOT

if [ "${execute}" == 1 ]
then
echo "Excuting: ${file}"
sbatch ${file}
fi

done;
done;
