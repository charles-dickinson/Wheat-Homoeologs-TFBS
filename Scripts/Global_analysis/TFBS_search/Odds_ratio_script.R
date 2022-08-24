### looking to see whether the homoeolog with TFmatrixID_1146
### has higher expression than the other homoeologs lacking TFmatrixID_1146
library(data.table)
library(ggplot2)
dat <- fread('Data/TFBS_search/Used/Produced/unique_hits_analysis/fimo_unique_v1.0.tsv')
motifs <- fread('Data/TFBS_search/Used/Produced/significant_motifs.csv')$x
dat <- unique(dat[motif_id %in% motifs,.(motif_id,sequence_name)][,gene := substring(sequence_name,1,18)][,.(motif_id,gene)],by=c('gene'))

triad_mapping <- fread('Data/Categorising_triads/Used/Produced/CS_no_stress_data.csv')
expression_data <- triad_mapping[,.(group_id,value,gene,chr_group)]
setnames(expression_data,'group_id','triad')

joined <- dat[expression_data,on='gene'][is.na(motif_id),present := 'motif_not_present'][!is.na(motif_id),present := 'motif_present']
joined <- joined[is.na(motif_id),score := 0][!is.na(motif_id),score := 1]

get_group <- function(a_row,b_row,d_row,pair){
  if (pair == 'a_b'){
    if (unique(a_row$present)[1] == 'motif_present' & unique(b_row$present)[1] == 'motif_not_present'){
      group <- 1
    } else {
      group <- 2
    }
    if (unique(a_row$value)[1] > unique(b_row$value)[1]){
      score <- 1
    } else {
      score <- 0
    }
  }
  if (pair == 'b_a'){
    if (unique(b_row$present)[1] == 'motif_present' & unique(a_row$present)[1] == 'motif_not_present'){
      group <- 1
    } else {
      group <- 2
    }
    if (unique(b_row$value)[1] > unique(a_row$value)[1]){
      score <- 1
    } else {
      score <- 0
    }
  }
  if (pair == 'a_d'){
    if (unique(a_row$present)[1] == 'motif_present' & unique(d_row$present)[1] == 'motif_not_present'){
      group <- 1
    } else {
      group <- 2
    }
    if (unique(a_row$value)[1] > unique(d_row$value)[1]){
      score <- 1
    } else {
      score <- 0
    }
  }
  if (pair == 'd_a'){
    if (unique(d_row$present)[1] == 'motif_present' & unique(a_row$present)[1] == 'motif_not_present'){
      group <- 1
    } else {
      group <- 2
    }
    if (unique(d_row$value)[1] > unique(a_row$value)[1]){
      score <- 1
    } else {
      score <- 0
    }
  }
  if (pair == 'b_d'){
    if (unique(b_row$present)[1] == 'motif_present' & unique(d_row$present)[1] == 'motif_not_present'){
      group <- 1
    } else {
      group <- 2
    }
    if (unique(b_row$value)[1] > unique(d_row$value)[1]){
      score <- 1
    } else {
      score <- 0
    }
  }
  if (pair == 'd_b'){
    if (unique(d_row$present)[1] == 'motif_present' & unique(b_row$present)[1] == 'motif_not_present'){
      group <- 1
    } else {
      group <- 2
    }
    if (unique(d_row$value)[1] > unique(b_row$value)[1]){
      score <- 1
    } else {
      score <- 0
    }
  }
  result <- c(group,score)
  result
}


all_triads <- unique(joined$triad)
odds_df <- data.frame('','','','')
colnames(odds_df) <- c('triad','pair','group','score')
counter <- 0
for (i in all_triads){
  filtered <- joined[triad==i,]
  a <- filtered[chr_group == 'A',]
  b <- filtered[chr_group == 'B',]
  d <- filtered[chr_group == 'D',]
  triad_df <- data.frame('','','',0)
  colnames(triad_df) <- c('triad','pair','group','score')
  pairs <- c('a_b','b_a','a_d','d_a','b_d','d_b')
  for (j in pairs){
    group <- get_group(a,b,d,j)[1]
    triad <- i
    pair <- j
    score <- get_group(a,b,d,j)[2]
    row <- c(triad,pair,group,score)
    triad_df <- rbind(triad_df,row)
  }
  odds_df <- rbind(odds_df,triad_df[-1,])
  counter <- counter + 1
  print(counter)
}
odds_df <- data.table(odds_df[-1,])

generate_graphs <- function(odds_data,expression_data){
  pairs <- c('a_b','b_a','a_d','d_a','b_d','d_b')
  pairs_chrs <- list(c('A','B'),c('B','A'),c('A','D'),c('D','A'),c('B','D'),c('D','B'))
  combined_greater <- data.frame(0,0)
  colnames(combined_greater) <- c('group1','group2')
  combined_lesser <- data.frame(0,0)
  colnames(combined_lesser) <- c('group1','group2')
  for (i in seq(1,length(pairs))){
    df <- odds_data[pair==pairs[i],]
    with_greater_without <- df[score == 1,][,.(count=.N),by=group][,expression := 'with_greater']
    with_greater_without <- dcast(with_greater_without,expression~group,value.var='count')
    colnames(with_greater_without) <- c('expression','group1','group2')
    combined_greater <- rbind(combined_greater, with_greater_without[,c('group1','group2')])
    with_lesser_without <- df[score == 0,][,.(count=.N),by=group][,expression := 'with_lesser']
    with_lesser_without <- dcast(with_lesser_without,expression~group,value.var='count')
    colnames(with_lesser_without) <- c('expression','group1','group2')
    combined_lesser <- rbind(combined_lesser, with_lesser_without[,c('group1','group2')])
  }
  group1 <- c(sum(combined_greater$group1),sum(combined_lesser$group1))
  group2 <- c(sum(combined_greater$group2),sum(combined_lesser$group2))
  odds_table <- data.frame(group1, group2)
  print(odds_table)
  odds_ratio_1 <- (odds_table[1,1] * odds_table[2,2]) / (odds_table[2,1] * odds_table[1,2])
  odds_ratio_2 <- (odds_table[2,1] * odds_table[1,2]) / (odds_table[1,1] * odds_table[2,2])
  print(odds_ratio_1)
  print(odds_ratio_2)
  df_for_plot <- data.frame(0,0,0,'',0,0)
  colnames(df_for_plot) <- c('triad','motif_not_present','motif_present','group','present_log10','not_present_log10')
  for (i in seq(1,length(pairs))){
    df <- odds_data[pair==pairs[i],]
    df$triad <- as.integer(df$triad)
    pairs_with_expression <- df[group==1,][expression_data,on=.(triad)]
    pairs_with_expression <- pairs_with_expression[group==1,][chr_group %in% pairs_chrs[[i]],]
    wider <- dcast(pairs_with_expression, triad ~ present,value.var = 'value')
    wider <- wider[motif_present > motif_not_present,group:='Homoeolog with one of \nsignificant motifs has \nhigher expression'][! motif_present > motif_not_present, group:='Homoeolog without any \nsignificant motif has \nhigher or equal expression']
    wider$present_log10 <- log10(wider$motif_present + 1)
    wider$not_present_log10 <- log10(wider$motif_not_present + 1)
    df_for_plot <- rbind(df_for_plot,wider)
  }
  df_for_plot <- df_for_plot[-1,]
  graph <- ggplot(df_for_plot,aes(x = present_log10, y = not_present_log10, colour = group)) +
    geom_point() +
    geom_abline(intercept = c(0,0),slope=1) +
    annotate('text',x=1.25,y=0.5,label=paste('OR =',signif(odds_ratio_1,digits=3)),size=6) +
    annotate('text',x=0.5,y=1.5,label=paste('OR =',signif(odds_ratio_2,digits=3)),size=6) +
    labs(title='Homoeolog 1 has at least one hit for significant motifs \nHomoeolog 2 has no hits for significant motifs') +
    theme_classic() +
    scale_colour_brewer(palette = 'Set1',labels = c('Homoeolog 1 has higher expression',
                                                    'Homoeolog 2 has higher expression')) +
    xlab('Log10(TPM+1) of homoeolog 1') +
    ylab('Log10(TPM+1) of homoeolog 2') +
    theme(text = element_text(size=20),
          legend.title = element_blank(),
          legend.key.size = unit(1.5,'cm'),
          )
    
  ggsave('Figures/TFBS_figure/All_enriched_motifs/odds_graph_overall_fig.png',graph, height=10, width=15)
}
generate_graphs(odds_df, joined)
