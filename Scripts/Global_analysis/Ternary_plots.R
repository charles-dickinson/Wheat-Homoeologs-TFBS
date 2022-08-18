### Produce the ternary plots and ternary data for analysis

library(data.table)
library(ggplot2)
library(ggtern)
all_data <- data.table(readRDS('Data/Categorising_triads/Used/Datasets/Triads.rds'))

filter_data <- function(overall_data){
  filtered <- overall_data[dataset=='HC_CS_no_stress' & factor=='all_mean_filter', ]
  filtered
}

make_ternary_data <- function(no_stress_data){
  no_stress_data$description <- as.character(no_stress_data$description)
  triads <- unique(no_stress_data$group_id)
  tern <- data.frame(0,0,0,0,'','','','',0,0)
  colnames(tern) <- c('triad','A','B','D','group','geneA','geneB','geneD','distance','expression_sum')
  for (i in 1:length(triads)){
    print(i)
    by_triad <- no_stress_data[group_id == triads[i], ]
    A_row <- by_triad[chr_group == 'A', ]
    geneA <- A_row$gene
    A <- A_row$normalised_triad
    B_row <- by_triad[chr_group == 'B', ]
    geneB <- B_row$gene
    B <- B_row$normalised_triad
    D_row <- by_triad[chr_group == 'D', ]
    geneD <- D_row$gene
    D <- D_row$normalised_triad
    group <- A_row$description
    distance <- unique(by_triad[[group]])
    expression_sum <- unique(by_triad$triad_sum)
    tern <- rbind(tern, c(triads[i],A,B,D,group,geneA,geneB,geneD,distance,expression_sum))
  }
  tern$A <- as.numeric(tern$A)
  tern$B <- as.numeric(tern$B)
  tern$D <- as.numeric(tern$D)
  tern <- tern[-1,]
}

make_ternary_plots <- function(ternary_data){
  ternary_data <- data.table(ternary_data)
  simple_mapping <- list(c('A.dominant','A.suppressed'),c('B.dominant','B.suppressed'),c('D.dominant','D.suppressed'),c('Central'))
  simple_groups <- c('A_different','B_different','D_different','Central')
  for (i in 1:4){
    ternary_data <- ternary_data[group %in% simple_mapping[[i]],simple_group := simple_groups[i]]
  }
  print(ternary_data)
  unbalanced <- c('A_different','B_different','D_different')
  ternary_data <- ternary_data[simple_group %in% unbalanced,balance:= 'Unbalanced'][!simple_group %in% unbalanced,balance := 'Balanced']
  order <- c('Central','A.dominant', 'A.suppressed','B.dominant', 'B.suppressed','D.dominant', 'D.suppressed')
  simple_order <- c('Central', 'A_different', 'B_different', 'D_different')
  balance_order <- c('Balanced', 'Unbalanced')
  ternary_data$group <- factor(ternary_data$group, levels = order)
  ternary_data$simple_group <- factor(ternary_data$simple_group, levels = simple_order)
  ternary_data$balance <- factor(ternary_data$balance, levels = balance_order)
  simple_group_tern <- ggtern(data = ternary_data, aes(A, B, D, color = simple_group)) +
    geom_point() +
    scale_color_brewer(palette = 'Set1',name = 'Triad Category')
  simple_group_tern
}

CS_no_stress_data <- filter_data(all_data)
tern_data <- make_ternary_data(CS_no_stress_data)
write.csv(tern_data, 'Data/Categorising_triads/Used/Produced/ternary_data_CS2.csv')
overall_tern_plot <- make_ternary_plots(tern_data)

extreme_triads <- tern_data[tern_data$distance < 0.05, ]
extreme_tern_plot <- make_ternary_plots(extreme_triads)
extreme_tern_plot