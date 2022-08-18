### Used
### takes Plant PAN family data and maps it to create a 
### file with the motif id, locus and family

library(data.table)

mapping <- fread('Data/TFBS_search/Used/Datasets/PlantPAN/ID_mapping_all_plant.txt')[,.(motif_id, gene)]

families <- fread('Data/TFBS_search/Used/Datasets/PlantPAN/All_species_tfs.txt')[,.(TF_ID, TF_family)][,.(TF_ID = gsub(' ','',TF_ID, fixed = TRUE),TF_family)]

setnames(families, 'TF_ID', 'gene')
setnames(families, 'TF_family', 'family')
mapping_families <- families2[mapping2, on=.(gene), nomatch=0]

fwrite(mapping_families, 'Data/TFBS_search/Used/Produced/PlantPAN/motif_families_all_sp.txt')