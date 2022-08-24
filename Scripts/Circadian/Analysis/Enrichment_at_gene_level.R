### this script represents the enrichment at gene level analysis
library(data.table)
library(tibble)
# only select genes that have modules
just_modules <- function(circadian){
  simple_circadian <- circadian[Ex2mrg_cluster_id != '#N/A',.(gene, Ex2mrg_cluster_id)]
  simple_circadian
}

# merge the motif hits and the module data
merge_fimo_modules <- function(fimo, modules){
  merged <- modules[fimo, on = .(gene),nomatch = NULL]
  merged
}

# work out the overall proportion, proportion per module and percentage enrichment
calculate_enrichment_motifs <- function(motif_and_phase_per_gene, circadian_dat_simple){
  number_of_genes <- length(unique(circadian_dat_simple$gene))
  overall_proportions <- motif_and_phase_per_gene[,.(proportion_of_genes_containing_motif=.N/number_of_genes), by=motif_id]
  n_genes_by_module <- circadian_dat_simple[,.(n_genes_per_module = .N), by = Ex2mrg_cluster_id]
  module_numbers <- motif_and_phase_per_gene[,.(motif_hits_per_phase = .N), by=.(motif_id, Ex2mrg_cluster_id)]
  module_proportions <- module_numbers[n_genes_by_module, on=.(Ex2mrg_cluster_id), 
                                       nomatch = NULL][,proportion_per_module := motif_hits_per_phase/n_genes_per_module]
  enrichment <- overall_proportions[module_proportions, on =.(motif_id)][,enrichment := 100*proportion_per_module/proportion_of_genes_containing_motif]
  enrichment
}

# most enriched stats
generate_list_of_sig_enriched_motifs <- function(motif_hits_by_phase, genes_module_list){
  total_genes <- length(unique(genes_module_list$gene))
  motif_counts_all_genes <- motif_hits_by_phase[, .(motif_in = .N, motif_not_in = total_genes - .N), by=motif_id]
  n_genes <- genes_module_list[,.(genes_in_module = .N),by=Ex2mrg_cluster_id]
  module_counts_per_motif <- unique(motif_hits_by_phase[n_genes, on = .(Ex2mrg_cluster_id)][,.(motif_in = .N, motif_not_in = genes_in_module - .N),
                                                                                            by = .(motif_id, Ex2mrg_cluster_id)])
  motif_list <- unique(module_counts_per_motif$motif_id)
  df <- data.frame("", "", 0)
  names(df) <- c('motif_id', 'module', 'enrichment') ### delete enrichment
  counter <- 0
  for (i in motif_list){
    for (j in unique(module_counts_per_motif$Ex2mrg_cluster_id)){
      contingency <- module_counts_per_motif[motif_id == i & Ex2mrg_cluster_id == j, .(motif_in,motif_not_in)]
      total <- rbind(contingency, motif_counts_all_genes[motif_id == i,.(motif_in,motif_not_in)])
      cont_table <- as.table(as.matrix(total))
      fisher_result <- fisher.test(cont_table,alternative = 'greater')
      if (fisher_result$p.value < 0.001){
        df <- add_row(df, motif_id = i, module = j, enrichment = fisher_result$p.value) ### delete enrichment
      }
    }
    counter <- counter + 1
    print(counter)
  }
  df
}

fimo_simple <- subset(unique(fread('data/fimo_results/fimo_v1.1_all_module_genes_simple.csv')))

circadian_dat <- fread('data/circadian/datasets/circadian_data.txt')

circadian_dat_simplified <- just_modules(circadian_dat)
merged_dat <- merge_fimo_modules(fimo_simple, circadian_dat_simplified)
module_enrichment_per_motif <- unique(calculate_enrichment_motifs(merged_dat, circadian_dat_simplified))
output <- generate_list_of_sig_enriched_motifs(merged_dat, circadian_dat_simplified)[-1,]
fwrite(output, 'data/circadian/produced/all_module_genes/enriched_at_gene_level.csv')
