### compare the variation in the enriched motifs between triads with balanced phase and those with unbalanced phase (>4h lag)

library(data.table)

balanced_phase_analysis <- function(balanced_data, filtered_fimo, dataset){
  with_balance <- balanced_data[is.na(balance_raw), balance := 'Balanced']
  with_balance <- balanced_data[!is.na(balance_raw), balance := 'Unbalanced']
  fimo_with_balance <- filtered_fimo[with_balance, on='gene']
  fimo_with_balance$chr <- substring(fimo_with_balance$gene, 9, 9)
  just_triads <- na.omit(fimo_with_balance, 'triad')
  by_triad <- unique(just_triads, by=c('gene', 'triad', 'motif_id', 'balance', 'chr'))[,.(category = sort(paste(chr, collapse = ''))),by=.(triad,motif_id,balance)]
  
  categories_list <- c('ABD','ADB','BAD','BDA','DAB','DBA','AB','BA','AD','DA','BD','DB','A','B','D')
  simple_categories <- c('No_variation','No_variation','No_variation','No_variation','No_variation','No_variation', 
                         'Variation','Variation','Variation','Variation','Variation','Variation','Variation','Variation','Variation')
  for (i in seq(1,length(categories_list))){
    by_triad <- by_triad[category == categories_list[i], simple_category := simple_categories[i]]
  }

  proportion_variation <- unique(by_triad[,.(n = .N,simple_category,balance),by=triad][,.(proportion = .N/n,balance),by=.(triad,simple_category)])
  no_variation <- proportion_variation[simple_category == 'No_variation',]
  balance_bplot <- ggplot(no_variation, aes(x = balance, y = proportion, fill = balance)) +
    geom_boxplot(alpha = 0.5) +
    geom_segment(aes(x = 1, xend = 1, y = 0.85, yend = 0.86)) +
    geom_segment(aes(x = 2, xend = 2, y = 0.85, yend = 0.86)) +
    geom_segment(aes(x = 1, xend = 2, y = 0.86, yend = 0.86)) +
    annotate('text', x = 1.5, y = 0.9, label = paste('p = ',round(kruskal.test(proportion~balance, no_variation)$p.value,3),sep=''),size = 5) +
    theme_classic() +
    theme(text = element_text(size = 20),
          legend.position = 'none',
          axis.title.x = element_blank()) +
    ylab('Proportion of enriched motifs conserved across homoeologs') +
    scale_fill_brewer(palette = 'Set1')
  if (dataset == 'Enriched_motifs'){
    ggsave('Figures/Circadian/Balanced_phase_boxplot.png', height = 10, width = 7.5)
  }
  if (dataset == 'Other_motifs'){
    ggsave('Figures/Circadian/Balanced_phase_boxplot_other_motifs.png', height = 10, width = 7.5)
    
  }
  print(kruskal.test(proportion ~ balance, data = no_variation))
  numbers <- no_variation[,.N,by=balance]
  numbers
}

circadian_dat <- fread('data/circadian/datasets/circadian_data.txt')
circadian_dat_filtered <- circadian_dat[Ex2mrg_cluster_id != '#N/A',]

fimo_results <- fread('data/circadian/produced/all_module_genes/fimo_v1.1_all_module_genes.tsv')
enriched_motifs <- fread('data/circadian/produced/all_module_genes/enriched_motifs_grouped.txt')
fimo_filtered <- fimo_results[motif_id %in% enriched_motifs$Motif_id,]
fimo_filtered$gene <- substring(fimo_filtered$sequence_name, 1, 18)

simple_circadian <- circadian_dat_filtered[,.(gene, triad, Unbalanced_phase)]
setnames(simple_circadian, 'Unbalanced_phase', 'balance_raw')

balanced_phase_analysis(simple_circadian, fimo_filtered, 'Enriched_motifs')

fimo_other_motifs <- fimo_results[!motif_id %in% enriched_motifs$Motif_id,]
fimo_other_motifs$gene <- substring(fimo_other_motifs$sequence_name, 1, 18)
balanced_phase_analysis(simple_circadian, fimo_other_motifs, 'Other_motifs')
