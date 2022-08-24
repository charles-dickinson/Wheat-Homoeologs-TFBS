### two rhythmic triads analysis
library(data.table)

only_triads <- function(full_dataset){
  triad_genes <- full_dataset[!is.na(triad), .(gene,triad,Ex2mrg_cluster_id)]
  fwrite(triad_genes[,.(gene,triad)], 'Data/Circadian/Used/Produced/triad_mapping.csv')
  triad_genes
}

generate_two_module_data <- function(triad_data){
  two_or_one_modules_triad_list <- triad_data[Ex2mrg_cluster_id != '#N/A',
                                              .N, by=triad][N == 2 | N ==1,]
  true_triads <- triad_gene_data[triad %in% two_or_one_modules_triad_list$triad,
                                            .N, by=triad][N == 3]
  only_true_triads <- triad_gene_data[triad %in% true_triads$triad,]
  only_true_triads$chr <- substring(only_true_triads$gene, 9, 9)
  setnames(only_true_triads, 'Ex2mrg_cluster_id','module')
  fwrite(only_true_triads, 'Data/Circadian/Used/Produced/Two_one_module_genes/two_one_modules_phase_per_gene.csv')
  only_true_triads
}

generate_module_categories <- function(two_one_modules_data, module_genes){
  phase_categs <- two_one_modules_data[module_genes, on=.(gene), nomatch = NULL][,.(gene,triad,module,chr)]
  phases <- phase_categs[,.(phase_list = list(unique(module))), by=triad]
  modules = c('W1','W2','W3','W4','W5','W6','W7','W8','W9')
  phases <- phases[,phase_same := fifelse(phase_list %in% modules, phase_list, list('NA'))][,.(triad,phase_same)]
  phase_df <- phase_categs[,.(chr_phases = sort(paste(chr, collapse = ''))),by=triad][phases, on=.(triad)]
  final <- phase_df[chr_phases == 'A', phase_category := 'A_module']
  final <- phase_df[chr_phases == 'B', phase_category := 'B_module']
  final <- phase_df[chr_phases == 'D', phase_category := 'D_module']
  final <- phase_df[chr_phases == 'BD' | chr_phases == 'DB', phase_category := 'A_no_module']
  final <- phase_df[(chr_phases == 'BD' | chr_phases == 'DB') & phase_same == 'NA', phase_category := 'A_no_two_diff']
  final <- phase_df[chr_phases == 'AD' | chr_phases == 'DA', phase_category := 'B_no_module'] 
  final <- phase_df[(chr_phases == 'AD' | chr_phases == 'DA') & phase_same == 'NA', phase_category := 'B_no_two_diff']
  final <- phase_df[chr_phases == 'AB' | chr_phases == 'BA', phase_category := 'D_no_module']
  final <- phase_df[(chr_phases == 'AB' | chr_phases == 'BA') & phase_same == 'NA', phase_category := 'D_no_two_diff']
  final <- final[,.(triad,phase_category,phase_same)]
  final
}

circadian_dat <- fread('Data/Circadian/Used/Datasets/Circadian_data.txt')
module_gene_list <- fread('Data/Circadian/Used/Produced/All_module_genes/all_module_genes_list.csv')

triad_gene_data <- only_triads(circadian_dat)
two_or_one_modules_dat <- generate_two_module_data(triad_gene_data)
phase_categs <- generate_module_categories(two_or_one_modules_dat, module_gene_list)

fwrite(phase_categs, 'Data/Circadian/Used/Produced/Two_one_module_genes/two_one_modules_phase_categories.csv')