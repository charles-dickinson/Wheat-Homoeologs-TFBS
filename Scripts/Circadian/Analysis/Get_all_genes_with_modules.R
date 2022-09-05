### this script takes the total circadian dataset and generates a file containing
### a list of the genes that show modules, and what module they are. This is for the gene-level 
### module analysis

library(data.table)

circadian_dat <- fread('data/circadian/datasets/circadian_data.txt')

generate_all_module_genes_list <- function(data){
  fwrite(data[Ex2mrg_cluster_id != '#N/A', .(gene,Ex2mrg_cluster_id)], 
         'data/circadian/produced/all_module_genes/all_module_genes_list.csv')
}

generate_all_module_genes_list(circadian_dat)
