###looks at the specificity of the conservation difference

library(data.table)
library(ggplot2)
library(FSA)
circadian_data <- fread('data/circadian/datasets/circadian_data.txt')
two_one_modules_dat <- fread('data/circadian/produced/two_one_module_genes/two_one_modules_phase_categories.csv')

generate_conserved_modules_list <- function(circadian_dat){
  circadian_filtered <- circadian_dat[Ex2mrg_cluster_id != '#N/A',]
  by_triad <- circadian_filtered[,.(n=.N,Ex2mrg_cluster_id),by=.(triad)][n==3,]
  only_conserved <- unique(by_triad,by=c('triad','Ex2mrg_cluster_id'))[,.(modules = .N),by=.(triad)][modules == 1,]$triad
  genes_list <- circadian_dat[triad %in% only_conserved,.(gene,triad)]
  genes_list
}

balanced_phase_analysis <- function(module_conservation_data, filtered_fimo){
  fimo_with_conservation <- filtered_fimo[module_conservation_data, on='gene']
  fimo_with_conservation$chr <- substring(fimo_with_conservation$gene, 9, 9)
  just_triads <- na.omit(fimo_with_conservation, 'triad')
  by_triad <- unique(just_triads, by=c('gene', 'triad', 'motif_id', 'module_cat', 'chr'))[,.(category = sort(paste(chr, collapse = ''))),by=.(triad,motif_id,module_cat)]
  
  categories_list <- c('ABD','ADB','BAD','BDA','DAB','DBA','AB','BA','AD','DA','BD','DB','A','B','D')
  simple_categories <- c('No_variation','No_variation','No_variation','No_variation','No_variation','No_variation', 
                         'Variation','Variation','Variation','Variation','Variation','Variation','Variation','Variation','Variation')
  for (i in seq(1,length(categories_list))){
    by_triad <- by_triad[category == categories_list[i], simple_category := simple_categories[i]]
  }
  
  proportion_variation <- unique(by_triad[,.(n = .N,simple_category,module_cat),by=triad][,.(proportion = .N/n,module_cat),by=.(triad,simple_category)])
  no_variation <- proportion_variation[simple_category == 'No_variation',]

  print(kruskal.test(proportion~module_cat,no_variation))
  print(dunnTest(proportion~module_cat,no_variation,method='bonferroni'))
  multiple_tests <- dunnTest(proportion~module_cat,no_variation,method='bonferroni')$res
  no_variation$module_cat <- factor(no_variation$module_cat, levels = c('All_same','A','B','D'))
  balance_bplot <- ggplot(no_variation, aes(x = module_cat, y = proportion, fill = module_cat)) +
    geom_boxplot(linetype = 'dashed',alpha = 0.5) +
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..,alpha=0.5)) +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..,alpha=0.5),width=0.5) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..,alpha=0.5),width=0.5) +
    geom_segment(aes(x=1,xend=1,y=0.89,yend=0.9)) +
    geom_segment(aes(x=2,xend=2,y=0.89,yend=0.9)) +
    geom_segment(aes(x=1,xend=2,y=0.9,yend=0.9)) +
    geom_segment(aes(x=1,xend=1,y=0.99,yend=1.0)) +
    geom_segment(aes(x=3,xend=3,y=0.99,yend=1.0)) +
    geom_segment(aes(x=1,xend=3,y=1.0,yend=1.0)) +
    geom_segment(aes(x=1,xend=1,y=1.09,yend=1.1)) +
    geom_segment(aes(x=4,xend=4,y=1.09,yend=1.1)) +
    geom_segment(aes(x=1,xend=4,y=1.1,yend=1.1)) +
    annotate('text',x=1.5,y=0.92,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'A - All_same',]$P.adj,digits=3)),size=7.5) +
    annotate('text',x=2,y=1.02,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'All_same - B',]$P.adj,digits=3)),size=7.5) +
    annotate('text',x=2.5,y=1.12,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'All_same - D',]$P.adj,digits=3)),size=7.5) +
    theme_classic() +
    theme(legend.position = 'none',
          text = element_text(size = 20)) +
    ylab('Proportion of enriched motifs conserved across homoeologs') +
    xlab('Subgenome gaining/losing module') +
    scale_fill_brewer(palette='Set1') 
#  ggsave('Figures/Circadian/subgenome_variation_boxplot.png', balance_bplot, height = 10, width = 10)
  ggsave('figures/circadian/random_subgenome_variation_boxplot.png', balance_bplot, height = 10, width = 10)
}

conserved_genes <- generate_conserved_modules_list(circadian_data)
conserved_genes$module_cat <- 'All_same'

two_one_modules_genes <- circadian_data[two_one_modules_dat,.(gene,triad,phase_category),on=.(triad)]
two_one_modules_genes$module_cat <- paste(substring(two_one_modules_genes$phase_category, 1,1), sep='')

two_one_modules_genes <- two_one_modules_genes[,.(gene,triad,module_cat)]
conserved_modules_dat <- rbind(two_one_modules_genes, conserved_genes)

fimo_results_all <- fread('data/circadian/produced/all_module_genes/fimo_results_all_module_genes.tsv')
fimo_results <- fread('data/circadian/produced/two_one_module_genes/fimo_results_two_one_module_genes.tsv')
enriched_motifs <- fread('data/circadian/produced/all_module_genes/enriched_motifs_grouped.txt')

fimo_combined <- unique(rbind(fimo_results, fimo_results_all))
fimo_filtered <- fimo_combined[motif_id %in% enriched_motifs$Motif_id,]
fimo_filtered$gene <- substring(fimo_filtered$sequence_name, 1, 18)

balanced_phase_analysis(conserved_modules_dat, fimo_filtered)


fimo_random <- fread('data/circadian/produced/all_module_genes/random_seq_fimo.tsv')
fimo_random_two_one <- fread('data/circadian/produced/two_one_module_genes/random_seq_fimo.tsv')
fimo_random_combined <- unique(rbind(fimo_random,fimo_random_two_one))
fimo_random_combined$gene <- substring(fimo_random_combined$sequence_name, 1, 18)
balanced_phase_analysis(conserved_modules_dat, fimo_random_combined)
