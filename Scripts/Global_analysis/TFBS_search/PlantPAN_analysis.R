### this script takes a unique motif hits file and produces
### two dataframes, the overall data and the family data.

library(data.table)

generate_triad_data <- function(fimo, family_mapping, triad_mapping){
  unique_hits_with_families <- motif_hits[unique(family_mapping, by=c('motif_id')), on=.(motif_id), nomatch=0][,n := 1]
  unique_hits_with_families$gene <- substring(unique_hits_with_families$sequence_name, 1, 18)
  by_triad_motifs <- unique_hits_with_families[triad_mapping, on=.(gene),nomatch=0]
  by_triad_motifs <- unique(by_triad_motifs[,.(category = paste(chr,collapse = ''),family),by=.(triad,motif_id)])
  categs_raw <- c('ABD','ADB','ABD','BAD','DBA','DAB','AB','AD','BD','BA','DA','DB','A','B','D')
  categs <- c('all_same','all_same','all_same','all_same','all_same','all_same','D_diff','B_diff','A_diff','D_diff','B_diff','A_diff','A_diff','B_diff','D_diff')
  for (i in seq(1,length(categs_raw))){
    by_triad_motifs <- by_triad_motifs[category == categs_raw[i], simple_category := categs[i]]
  }
  by_triad_motifs
}
generate_motif_data <- function(by_triad_df, triad_category_df){
  by_triad_df <- by_triad_df[simple_category == 'all_same',variation:='No_variation']
  by_triad_df <- by_triad_df[simple_category != 'all_same',variation:='Variation']
  with_balance <- by_triad_df[triad_category_df,on=.(triad)]
  descriptions <- c('A.dominant','B.dominant','D.dominant','A.suppressed','B.suppressed','D.suppressed','Central')
  balances <- c('Unbalanced','Unbalanced','Unbalanced','Unbalanced','Unbalanced','Unbalanced', 'Balanced')
  simple_descriptions <- c('A','B','D','A','B','D','Central')
  for (i in seq(1,length(descriptions))){
    with_balance <- with_balance[description==descriptions[i],balance:=balances[i]][description==descriptions[i],simple_description:=simple_descriptions[i]]
  }
  by_motif <- with_balance[,.(count = .N),by=.(motif_id,balance,variation)]
  fwrite(by_motif,'Data/TFBS_search/Used/Produced/PlantPAN/by_motif_variation_data.csv')
}

generate_by_family_data <- function(by_triad_data,triad_category_data){
  grouped <- by_triad_data[,.(n =.N),by=.(triad,family,simple_category)]
  grouped <- grouped[triad_category_data,.(triad,family,n,description,simple_category),on=.(triad)]
  descriptions <- c('A.dominant','B.dominant','D.dominant','A.suppressed','B.suppressed','D.suppressed','Central')
  balances <- c('Unbalanced','Unbalanced','Unbalanced','Unbalanced','Unbalanced','Unbalanced', 'Balanced')
  simple_descriptions <- c('A','B','D','A','B','D','Central')
  for (i in seq(1,length(descriptions))){
    grouped <- grouped[description==descriptions[i],balance:=balances[i]][description==descriptions[i],simple_description:=simple_descriptions[i]]
  }
  grouped <- grouped[,.(norm = n/sum(n),n,simple_category,balance,description,simple_description),by=.(triad,family)]
  grouped_variation_data <- grouped[simple_category == 'all_same',variation := 'no_variation'][simple_category != 'all_same',variation := 'variation']
  grouped_variation_data <- grouped_variation_data[,.(prop_variation = sum(norm),balance,description,simple_description,simple_category),by=.(triad,variation,family)]
  fwrite(grouped, 'Data/TFBS_search/Used/Produced/PlantPAN/full_family_data.csv')
  fwrite(grouped_variation_data, 'Data/TFBS_search/Used/Produced/PlantPAN/full_family_variation_data.csv')
}

generate_overall_data <- function(by_triad_dat, triad_category_dat){
  overall <- by_triad_dat[,.(n=.N),by=.(triad,simple_category)][triad_category_dat,.(triad,n,description,simple_category),on=.(triad)]
  descriptions <- c('A.dominant','B.dominant','D.dominant','A.suppressed','B.suppressed','D.suppressed','Central')
  balances <- c('Unbalanced','Unbalanced','Unbalanced','Unbalanced','Unbalanced','Unbalanced', 'Balanced')
  simple_descriptions <- c('A','B','D','A','B','D','Central')
  for (i in seq(1,length(descriptions))){
    overall <- overall[description==descriptions[i],balance:=balances[i]][description==descriptions[i],simple_description:=simple_descriptions[i]]
  }
  overall <- overall[,.(norm = n/sum(n),n,simple_category,balance,description,simple_description),by=.(triad)]
  variation_data <- overall[simple_category == 'all_same',variation := 'no_variation'][simple_category != 'all_same',variation := 'variation']
  variation_data <- variation_data[,.(prop_variation = sum(norm),balance,description,simple_description,simple_category),by=.(triad,variation)]
  fwrite(overall, 'Data/TFBS_search/Used/Produced/PlantPAN/full_data.csv')
  fwrite(variation_data, 'Data/TFBS_search/Used/Produced/PlantPAN/full_variation_data.csv')
}


motif_hits <- fread('Data/TFBS_search/Used/Produced/unique_hits_analysis/fimo_unique_v1.0.tsv')
families <- fread('Data/TFBS_search/Used/Produced/PlantPAN/motif_families_all_sp.txt')
setnames(families, 'gene', 'locus')

triads <- fread('Data/TFBS_search/Used/Produced/triads.txt')
triad_categs <- fread('Data/Categorising_triads/Used/Produced/ternary_data_CS.csv') # normal
triad_categs_simple <- triad_categs[,.(triad,group)]
setnames(triad_categs_simple,'group','description')

by_triad <- generate_triad_data(motif_hits, families, triads)
generate_by_family_data(by_triad,triad_categs_simple)
generate_overall_data(by_triad,triad_categs_simple)
generate_motif_data(by_triad,triad_categs_simple)
