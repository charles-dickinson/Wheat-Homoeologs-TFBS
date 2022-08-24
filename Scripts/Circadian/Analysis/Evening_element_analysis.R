library(data.table)
library(ggplot2)

phase_analysis <- function(fimo_results_filtered, phase_data){
  presence_evening <- unique(fimo_results_filtered, by='gene', with = F)
  phase_contains_evening <- phase_data[gene %in% presence_evening$gene, ]
  phase_no_evening <- phase_data[!gene %in% presence_evening$gene,]
  
  phase_no_evening$dataset <- "Genes not containing EE"
  phase_contains_evening$dataset <- "Genes containing EE"
  phase_bplot <- ggplot() +
    geom_violin(data = phase_no_evening, aes(y = dataset, x = CT, fill = dataset), alpha = 0.5) +
    geom_violin(data = phase_contains_evening, aes(y = dataset, x = CT, fill= dataset), alpha = 0.5) +
    geom_boxplot(data = phase_contains_evening, aes(y = dataset, x = CT),
                 width = 0.25) +
    geom_boxplot(data = phase_no_evening, aes(y = dataset, x = CT),
                 width = 0.25) +
    theme_classic() +
    theme(text = element_text(size = 20),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      legend.position = 'top',
    ) +
    xlab('Phase (hours after dawn)') +
    scale_fill_brewer(palette = 'Set1')
  phase_bplot
  ggsave('Figures/Circadian/EE_phase_violin.png', height = 5, width = 7.5)
  
  phase_stats_result <- wilcox.test(phase_no_evening$CT, 
                                    phase_contains_evening$CT,
                                    alternative = 'greater')
  phase_stats_result
}

module_enrichment <- function(circadian_dat, fimo_just_EE){
  just_module <- circadian_dat[Ex2mrg_cluster_id != '#N/A',.(gene,Ex2mrg_cluster_id)]
  presence_evening <- unique(fimo_just_EE, by='gene', with = F)
  fimo_0461_with_module <- just_module[gene %in% presence_evening$gene,]
  all_genes_module_proportions <- just_module[, .(all_genes_prop = .N/length(just_module$gene)), by=Ex2mrg_cluster_id]
  all_genes_module_proportions$dataset <- 'all_genes'
  id_0461_module_proportions <- fimo_0461_with_module[, .(all_genes_prop = .N/length(fimo_0461_with_module$gene)), by=Ex2mrg_cluster_id]
  id_0461_module_proportions$dataset <- 'genes_containing_EE'
  total_for_plot <- rbind(all_genes_module_proportions, id_0461_module_proportions)
  total_for_plot$dataset <- factor(total_for_plot$dataset, levels = c('genes_containing_EE','all_genes'))
  module_plot <- ggplot(total_for_plot, aes(x = Ex2mrg_cluster_id, y = all_genes_prop, fill = dataset)) +
    geom_col(position='dodge', alpha = 0.5) +
    annotate('text',x=4, y = 0.32, label='*', size = 10) +
    annotate('text',x=5, y = 0.32, label='*', size =10) +
    theme_classic() +
    theme(text = element_text(size = 20),
          legend.position = 'top') +
    ylab('Proportion of genes') +
    xlab('Phase Module') +
    scale_fill_brewer(palette = 'Set1', name = '', labels = c('Genes containing EE','All genes'))
  ggsave('Figures/Circadian/EE_module_proportions.png',module_plot, height = 5, width = 7.5)
  
  expected_cont_row <- data.frame(length(fimo_0461_with_module$gene)/9, 
                                  (length(just_module$gene) - length(fimo_0461_with_module$gene))/9)
  colnames(expected_cont_row) <- c('Motif_in', 'Motif_not_in')
  module_gene_counts <- just_module[,.(count = .N),by='Ex2mrg_cluster_id']
  module_numbers_0461 <- fimo_0461_with_module[,.(motif_in = .N),by='Ex2mrg_cluster_id']
  for (i in unique(fimo_0461_with_module$Ex2mrg_cluster_id)){
    motif_in <- module_numbers_0461[Ex2mrg_cluster_id == i, .(motif_in)]
    motif_not_in <- module_gene_counts[Ex2mrg_cluster_id == i, .(count)] - motif_in
    actual_cont_row <- data.frame(motif_in, motif_not_in)
    colnames(actual_cont_row) <- c('Motif_in', 'Motif_not_in')
    cont_table <- rbind(actual_cont_row, expected_cont_row)
    fisher_result <- fisher.test(cont_table,alternative = 'greater')
    # cont_table <- as.table(as.matrix(cont_table))
    print(i)
    print(fisher_result)
  }
}

triad_analysis <- function(circadian_dat, fimo_just_EE){
  triad_17859_circadian <- circadian_dat[triad == 17859,]
  triad_17859_circadian_simple <- triad_17859_circadian[,.(gene,X0,X4,X8,X12,X16,
                                                         X20,X24,X28,X32,X36,
                                                         X40,X44,X48,X52,X56,
                                                         X60,X64,X68)]
  longer <- melt(triad_17859_circadian_simple, id.vars = 'gene', measure.vars = c('X0','X4','X8','X12','X16','X20','X24','X28','X32','X36','X40','X44','X48','X52','X56','X60','X64','X68'),
                variable.name = 'Time', value.name = 'TPM')
  times <- c('X0','X4','X8','X12','X16','X20','X24','X28','X32','X36','X40','X44','X48','X52','X56','X60','X64','X68')
  numeric_times <- c(0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68)
  for (i in seq(1, length(times))){
    longer <- longer[Time == times[i], ZT := numeric_times[i]]
  }
  darkdf <- data.frame(c(18, 42, 66), c(5, 5, 5))
  colnames(darkdf) <- c('time','height')
  triad_plot <- ggplot(longer, aes(x = ZT, y = TPM)) +
    geom_col(data = darkdf, aes(x = time, y = height), width = 12, fill = 'darkgrey', 
             colour = 'white', alpha = 0.5) +
    geom_line(aes(colour = gene, size = 20)) +
    theme_classic() +
    theme(text = element_text(size = 30),
          legend.position = 'none') +
    scale_colour_brewer(palette = 'Set1') +
    xlab('Time (hours after transfer into L:L)')
  ggsave('Figures/Circadian/triad_17859_circadian.png',triad_plot,
         height = 10, width = 15)
}

circadian_data <- fread('Data/Circadian/Used/Datasets/Circadian_data.txt')
circadian_data$CT <- (circadian_data$meta2d_phase_ZT24.68 * 24)/circadian_data$meta2d_period_ZT24.68
for (i in c('0', '4', '8','12','16',
           '20','24','28','32','36',
           '40','44','48','52','56',
           '60','64', '68' )){
  setnames(circadian_data, i, paste('X',i,sep=''))
}
fimo_results_all_module_genes <- fread('Data/Circadian/Used/Produced/All_module_genes/fimo_results_all_module_genes.tsv')

just_phase <- circadian_data[ZT24.68_q.01 == 'Y',.(gene,CT)]

fimo_evening_element <- fimo_results_all_module_genes[motif_id == 'TF_motif_seq_0461',]
fimo_evening_element$gene <- substring(fimo_evening_element$sequence_name, 1, 18)
phase_analysis(fimo_evening_element, just_phase)
module_enrichment(circadian_data, fimo_evening_element)
triad_analysis(circadian_data, fimo_evening_element)
