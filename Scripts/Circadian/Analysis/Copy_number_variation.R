### copy number variation - do circadian unbalanced triads show higher copy number variation 
### for our enriched motifs???

library(data.table)
library(ggplot2)

balanced_phase_copy_number_analysis <- function(balanced_data, filtered_fimo, dataset){
  with_balance <- balanced_data[is.na(balance_raw), balance := 'Balanced']
  with_balance <- balanced_data[!is.na(balance_raw), balance := 'Unbalanced']
  fimo_with_balance <- filtered_fimo[with_balance, on='gene']
  fimo_with_balance$chr <- substring(fimo_with_balance$gene, 9, 9)
  just_triads <- na.omit(fimo_with_balance, 'triad')
#  pos_strand_hits <- just_triads[strand == '+',]
  pos_strand_hits <- just_triads
  copy_number_per_gene <- unique(pos_strand_hits[,.(n = .N, triad = triad,balance=balance), by=.(gene,motif_id)])
  grouped_per_triad <- copy_number_per_gene[,.(all_genes_contain = .N,n=n, balance = balance),by=.(triad,motif_id)][all_genes_contain == 3,.(homeologs=.N, balance=balance),by=.(triad,motif_id,n)]
  variation_per_triad <- unique(grouped_per_triad)[,.(copy_number_variation = .N,balance=balance),by=.(triad,motif_id)][copy_number_variation == 1,variation := 'No_variation'][copy_number_variation != 1, variation := 'Variation']
  proportion_variation <- unique(variation_per_triad[,.(n = .N,variation,balance),by=triad][,.(proportion = .N/n,balance),by=.(triad,variation)])
  no_variation <- proportion_variation[variation == 'No_variation',]
  variation_plot <- ggplot(no_variation, aes(x = balance, y = proportion, fill=balance)) +
    geom_boxplot(alpha = 0.5) +
    geom_segment(aes(x = 1, xend = 1, y = 1.05, yend = 1.06)) +
    geom_segment(aes(x = 2, xend = 2, y = 1.05, yend = 1.06)) +
    geom_segment(aes(x = 1, xend = 2, y = 1.06, yend = 1.06)) +
    annotate('text', x = 1.5, y = 1.1, label = paste('p = ',round(kruskal.test(proportion~balance, no_variation)$p.value,3),sep=''),size = 5) +
    theme_classic()+
    ylab('Proportion of motifs in triad showing conserved copy number') +
    xlab('') +
    theme(legend.position='none',
          text = element_text(size = 20)) +
    scale_fill_brewer(palette = 'Set1')
  if (dataset == 'Enriched_motifs'){
    ggsave('Figures/Circadian/Copy_number_variation_enriched_motifs.png', variation_plot, height = 10, width = 7.5)
  }
  if (dataset == 'Other_motifs'){
    ggsave('Figures/Circadian/Copy_number_variation_other_motifs.png', variation_plot, height = 10, width = 7.5)
    
  }
  print(kruskal.test(proportion~balance, no_variation))
  print(median(no_variation[balance == 'Balanced']$proportion))
  print(median(no_variation[balance == 'Unbalanced']$proportion))
  numbers <- no_variation[,.N,by=balance]
  print(numbers)
}

circadian_dat <- fread('data/circadian/datasets/circadian_data.txt')
circadian_dat_filtered <- circadian_dat[Ex2mrg_cluster_id != '#N/A',]

fimo_results <- fread('data/circadian/produced/all_module_genes/fimo_v1.1_all_module_genes.tsv')
enriched_motifs <- fread('data/circadian/produced/all_module_genes/enriched_motifs_grouped.txt')
fimo_filtered <- fimo_results[motif_id %in% enriched_motifs$Motif_id,]
fimo_filtered$gene <- substring(fimo_filtered$sequence_name, 1, 18)

simple_circadian <- circadian_dat_filtered[,.(gene, triad, Unbalanced_phase)]
setnames(simple_circadian, 'Unbalanced_phase', 'balance_raw')

balanced_phase_copy_number_analysis(simple_circadian, fimo_filtered, 'Enriched_motifs')

fimo_other_motifs <- fimo_results[!motif_id %in% enriched_motifs$Motif_id,]
fimo_other_motifs$gene <- substring(fimo_other_motifs$sequence_name, 1, 18)
balanced_phase_copy_number_analysis(simple_circadian, fimo_other_motifs, 'Other_motifs')
