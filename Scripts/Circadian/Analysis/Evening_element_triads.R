### In two-one module triads, check whether EE enriched in genes showing W4/W5/W6 modules

library(data.table)
library(ggplot2)

gene_module_data <- fread('data/circadian/produced/two_one_module_genes/two_one_modules_phase_per_gene.csv')
triad_module_data <- fread('data/circadian/produced/two_one_module_genes/two_one_modules_phase_categories.csv')
just_evening <- gene_module_data[triad_module_data,on=.(triad)][phase_same %in% c('W4','W5')]
motif_mapping <- fread('data/circadian/produced/all_module_genes/enriched_motifs_grouped.txt')
fimo_results <- fread('data/circadian/produced/two_one_module_genes/fimo_v1.1_two_one_module_genes_simple.csv')

fimo_results <- fimo_results[gene %in% just_evening$gene,]

ee_motifs <- unique(motif_mapping[Name=='EE',]$Motif_id)
fimo_results_EE <- fimo_results[motif_id == 'TF_motif_seq_0461',]

hits_with_phase <- unique(fimo_results_EE[just_evening,on=.(gene)],by='gene')

counts_by_module <- hits_with_phase[!is.na(motif_id),count := 1][!is.na(motif_id),present := 'Motif_in'][is.na(motif_id),count := 0][is.na(motif_id),present := 'Motif_not_in']
counts_by_module <- counts_by_module[module == '#N/A',module_present := 'Not_module'][module != '#N/A',module_present := 'Has_module']
grouped <- counts_by_module[,.(number = .N),by=.(module_present,present)]

cont_table <- dcast(grouped,present~module_present,value.var='number')
matrix <- as.matrix(cont_table, rownames='present')

chi_test <- chisq.test(matrix)
chi_test
grouped_with_totals <- grouped[,.(total = sum(number),present,number),by=.(module_present)]
grouped_with_totals <- grouped_with_totals[,prop := number/total]
motif_in <- grouped_with_totals[present == 'Motif_in',]
plot <- ggplot(motif_in,aes(x=module_present,y=prop, fill=module_present)) +
  geom_col(alpha=0.5) +
  ylab('Proportion of genes with EE') +
  geom_segment(aes(x = 1, xend=1, y=0.07, yend = 0.072)) +
  geom_segment(aes(x = 2, xend=2, y=0.07, yend = 0.072)) +
  geom_segment(aes(x = 1, xend=2, y=0.072, yend = 0.072)) +
  annotate('text',x=1.5,y=0.075,label = 'p = 0.00266', size =10) +
  theme_classic() +
  theme(legend.position='none',
        text = element_text(size = 20)) +
  xlab('') + 
  scale_x_discrete(labels = c('Homoeologs in evening modules','Homoeologs without evening modules'))
plot  

ggsave('Figures/Circadian/EE_homoeolog_conservation.png',height=12.5,width=10)
