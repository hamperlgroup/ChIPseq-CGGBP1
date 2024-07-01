## One liner to extract protein coding genes from genome gtf

# /lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/data/genome
# /lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1_insertion/data/genome

cat genome.gtf | awk '$3 == "gene"' | grep "protein_coding" | awk -v OFS='\t' '{print $1,$4,$5,$10,0,$7}' | sed 's/"//g' | sed 's/;//' > protein_coding_genes.bed