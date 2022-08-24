### trying to find the motifs with the biggest difference in proportion of conserved 
### hits in balanced vs unbalanced triads
library(data.table)
full_data <- fread('Data/TFBS_search/Used/Produced/PlantPAN/full_family_data.csv')

by_balance_and_variation <- full_data[,.(count = sum(n)),by=.(family,balance,variation)]

pvalues <- data.frame('',0)
names(pvalues) <- c('Family','P.value')
for(i in unique(by_balance_and_variation$family)){
  filtered <- by_balance_and_variation[family == i,]
  wider <- dcast(filtered,variation~balance,value.var='count')
  wider$Balanced <- as.numeric(wider$Balanced)
  wider$Unbalanced <- as.numeric(wider$Unbalanced)
  matrix <- as.matrix(wider,rownames='variation')
  p <- chisq.test(matrix)$p.value
  row <- c(i, p)
  print(i)
  print(matrix)
  pvalues <- rbind(pvalues, row)
}
pvalues <- pvalues[-1,]
pvalues$P.value <- p.adjust(pvalues$P.value, method='bonferroni')

by_motif_data <- fread('Data/TFBS_search/Used/Produced/PlantPAN/by_motif_variation_data.csv')
pvalues_motifs <- data.frame('',0)
names(pvalues_motifs) <- c('Motif_id','P.value')
for(i in unique(by_motif_data$motif_id)){
  filtered <- by_motif_data[motif_id == i,]
  wider <- dcast(filtered,variation~balance,value.var='count')
  wider$Balanced <- as.numeric(wider$Balanced)
  wider$Unbalanced <- as.numeric(wider$Unbalanced)
  matrix <- as.matrix(wider,rownames='variation')
  p <- chisq.test(matrix)$p.value
  row <- c(i, p)
  pvalues_motifs <- rbind(pvalues_motifs, row)
}
pvalues_motifs <- pvalues_motifs[-1,]
pvalues_motifs$P.value <- p.adjust(pvalues_motifs$P.value, method='bonferroni')
significant_motifs <- pvalues_motifs[pvalues_motifs$P.value < 0.05,]$Motif_id
write.table(significant_motifs, 'Data/TFBS_search/Used/Produced/significant_motifs.csv',quote=F)
