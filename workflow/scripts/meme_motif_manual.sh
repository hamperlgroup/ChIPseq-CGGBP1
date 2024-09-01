mamba activate motif
cd /lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1_insertion/test/meme

sed 's/chr//g' /lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/results/ChIP-peaks/narrowPeak/merged/ChIP-seq_CGGBP1_ENCODE.bed > ChIP-seq_CGGBP1_ENCODE_noCHR.bed

bedtools getfasta -fi /lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/data/genome/genome.fa -bed ChIP-seq_CGGBP1_ENCODE_noCHR.bed -fo ChIP-seq_CGGBP1_ENCODE.fa

fasta-shuffle-letters -k 2 -dna ChIP-seq_CGGBP1_ENCODE.fa ChIP-seq_CGGBP1_ENCODE_shuffled.fa

#with background sequences
meme-chip -o ChIP-seq_CGGBP1_ENCODE_Overshuffle -db /lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1_insertion/data/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme -neg ChIP-seq_CGGBP1_ENCODE_shuffled.fa -dna -meme-nmotifs 25 -minw 8 -maxw 12 -meme-mod zoops ChIP-seq_CGGBP1_ENCODE.fa

#without background data
meme-chip -o ChIP-seq_CGGBP1_ENCODE_noBck -db /lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1_insertion/data/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme -dna -meme-nmotifs 25 -minw 8 -maxw 12 -meme-mod zoops ChIP-seq_CGGBP1_ENCODE.fa
