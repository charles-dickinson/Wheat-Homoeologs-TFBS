### this script cleans up the fimo results

clean_fimo <- function(unique_fimo_results){
  setnames(unique_fimo_results, "sequence_name", "gene")
  unique_fimo_results$gene <- substr(unique_fimo_results$gene, 1, 9) # remove strand directionality (1-9 for At genes)
  fimo_simple <- unique_fimo_results[,.(gene, motif_id)] # remove alt_id
  fimo_simple
}

fimo <- fread("data/fimo_results/fimo_v1.1_all_module_genes.tsv")
#fimo <- fread("data/fimo_results/fimo_v1.1_two_one_module_genes.tsv")

fimo_results_simple <- clean_fimo(fimo)  
fwrite(fimo_results_simple, "data/fimo_results/fimo_v1.1_all_module_genes_simple.csv")
#fwrite(fimo_results_simple, "data/fimo_results/fimo_v1.1_two_one_module_genes_simple.csv")