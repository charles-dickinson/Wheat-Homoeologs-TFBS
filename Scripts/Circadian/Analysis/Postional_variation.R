library(data.table)
library(FSA)
combine_datasets <- function(motif_group_data, fimo_hits, circadian_dat){
  fimo_filtered <- fimo_hits[motif_group_data, on = 'motif_id']
  fimo_filtered$gene <- substring(fimo_filtered$sequence_name, 1, 18)
  with_module <- na.omit(circadian_data[fimo_filtered, on='gene'],cols=c('Ex2mrg_cluster_id'))
}

positional_variation <- function(combined_datasets, group_data){
  for (i in unique(groups$group)){
    fimo_group <- combined_datasets[group == i,]
    fimo_group <- setorder(fimo_group, gene, -start)
    first_occurrence <- unique(fimo_group, by=c('gene'))
    stats <- kruskal.test(start ~ Ex2mrg_cluster_id, data = first_occurrence)
    if (stats$p.value < 0.01){
      print(i)
      print(stats)
      dunn <- dunnTest(first_occurrence$start, first_occurrence$Ex2mrg_cluster_id,
                  method = 'bonferroni')$res
      signif_diffs <- dunn[dunn$P.adj < 0.05,]$Comparison
      print(signif_diffs)
    }      
    if (stats$p.value < 0.001){
      caption <- labs(title = paste(i,', Kruskal-Wallis p < 0.001', sep=''))
    }
    if (stats$p.value > 0.001){
      caption <- labs(title = paste(i,', Kruskal-Wallis p = ', signif(stats$p.value, digits=3),sep=''))
    }
    if (stats$p.value == 0.001){
      caption <- labs(title = paste(i,', Kruskall-Wallis p = ', signif(stats$p.value, digits=3),sep=''))
    }
    plot <- ggplot(first_occurrence, aes(y = Ex2mrg_cluster_id, x = start)) +
      geom_violin(alpha = 0.25, adjust = 1/4,aes(fill = Ex2mrg_cluster_id)) +
      geom_boxplot(fill = 'white',alpha=0.1,width=0.75,linetype='dashed') +
      annotate('text', x=1600,y=1, label = 'TSS',size=5) +
 #     geom_point(position = 'jitter', alpha = 0.1) +
      theme_classic() +
      theme(text = element_text(size = 20),
            axis.title.y = element_blank(),
            axis.text = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            legend.position = 'none') +
      xlab("Start position of motif along promoter (5' to 3')") +
      caption +
      scale_fill_brewer(palette = 'BrBG')
    ggsave(paste('Figures/Circadian/Groups/',i,'.png',sep=''), height = 7.5, width = 10)
  }
}

fimo <- fread('data/circadian/produced/all_module_genes/fimo_v1.1_all_module_genes.tsv')

circadian_data <- fread('data/circadian/produced/all_module_genes/all_module_genes_list.tsv')

motif_groups <- fread('data/circadian/produced/all_module_genes/enriched_motifs_grouped.txt')

groups <- unique(motif_groups)
setnames(groups, 'Name', 'group')
setnames(groups, 'Motif_id', 'motif_id')

combined <- combine_datasets(groups, fimo, circadian_data)
positional_variation(combined, groups)
