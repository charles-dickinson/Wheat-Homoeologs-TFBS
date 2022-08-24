### this script takes the total circadian dataset and generates a file containing
### a list of the genes that show modules, and what module they are. This is for the gene-level 
### module analysis

library(data.table)

circadian_dat <- fread('Data/Circadian/Used/Datasets/Circadian_data.txt')

generate_all_module_genes_list <- function(data){
  fwrite(data[Ex2mrg_cluster_id != '#N/A', .(gene,Ex2mrg_cluster_id)], 
         'Data/Circadian/Used/Produced/All_module_genes/all_module_genes_list.csv')
}

generate_all_module_genes_list(circadian_dat)

