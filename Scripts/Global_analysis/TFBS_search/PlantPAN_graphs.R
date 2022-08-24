### produces the graphs used in the report
library(data.table)
library(ggplot2)
library(FSA)
data <- fread('Data/TFBS_search/Used/Produced/PlantPAN/full_data.csv')
family_data <- fread('Data/TFBS_search/Used/Produced/PlantPAN/full_family_data.csv')

variation_data <- fread('Data/TFBS_search/Used/Produced/PlantPAN/full_variation_data.csv')
family_variation_data <- fread('Data/TFBS_search/Used/Produced/PlantPAN/full_family_variation_data.csv')

motif_categories <- c("all_same","A_diff","B_diff","D_diff")
expression_categories <- c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed", "D.dominant","D.suppressed")
balance_categories <- c("Balanced","Unbalanced")
simple_description_categories <- c("Central","A","B","D")
data$simple_category <- factor(data$simple_category, levels = motif_categories)
data$description <- factor(data$description, levels = expression_categories)
data$balance <- factor(data$balance, levels = balance_categories)
data$simple_description <- factor(data$simple_description, levels = simple_description_categories)

variation_data$simple_category <- factor(variation_data$simple_category, levels = motif_categories)
variation_data$description <- factor(variation_data$description, levels = expression_categories)
variation_data$balance <- factor(variation_data$balance, levels = balance_categories)
variation_data$simple_description <- factor(variation_data$simple_description, levels = simple_description_categories)

family_data$simple_category <- factor(family_data$simple_category, levels = motif_categories)
family_data$description <- factor(family_data$description, levels = expression_categories)
family_data$balance <- factor(family_data$balance, levels = balance_categories)
family_data$simple_description <- factor(family_data$simple_description, levels = simple_description_categories)

family_variation_data$simple_category <- factor(family_variation_data$simple_category, levels = motif_categories)
family_variation_data$description <- factor(family_variation_data$description, levels = expression_categories)
family_variation_data$balance <- factor(family_variation_data$balance, levels = balance_categories)
family_variation_data$simple_description <- factor(family_variation_data$simple_description, levels = simple_description_categories)

expression.colors = c('Central'="#9fa1a6", 'A.dominant'="#58a440", 'A.suppressed'="#a7d489",
                      'B.dominant'="#4a2773", 'B.suppressed'="#2bc0f8",'D.dominant'="#f9932e", 'D.suppressed'="#fccf2b")

simple_expression_colours <- c("Central" = "#9fa1a6", 'A'="#58a440", "B" = "#4a2773",
                               "D" = "#f9932e")

balance_colours <- c("Balanced" = "#2bc0f8", "Unbalanced" = "#fccf2b")
for (i in unique(data$simple_category)){
  print(i)
  filtered <- data[simple_category == i,]
  print(kruskal.test(norm~simple_description,filtered))
  print(dunnTest(norm~simple_description,filtered,method='bonferroni'))
}

all_cat_boxplot_norm <- ggplot(data = data, aes(y=norm, x = simple_category, fill = simple_description)) +
  geom_boxplot(alpha = 0.5,linetype='dashed') +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),alpha=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..),alpha=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..),alpha=0.5) +
  theme_classic() +
  ylab('Proportion of Motifs in Triad') +
  xlab('Motif Simple Category') +
  scale_x_discrete(labels = c('All Same', 'A Different', 'B Different', 'D Different')) +
  scale_fill_brewer(palette='Set1',name = "Triad Category", labels = c('Central', 'A dominant/ \nA suppressed', 
                                                     'B dominant/ \nB suppressed', 'D dominant/ \nD suppressed')) +
  theme(legend.position = 'top',
        legend.key.size = unit(0.8, 'cm'),
        text = element_text(size=20)) +
  annotate('text', x = 0.72, y = 1.05, label = 'ABD') +
  annotate('text', x = 0.9, y = 1.05, label = 'C') +
  annotate('text', x = 1.095, y = 1.05, label = 'C') +
  annotate('text', x = 1.28, y = 1.05, label = 'C') +
  annotate('text', x = 1.72, y = 1.05, label = 'A') +
  annotate('text', x = 1.9, y = 1.05, label = 'CBD') +
  annotate('text', x = 2.095, y = 1.05, label = 'A') +
  annotate('text', x = 2.28, y = 1.05, label = 'A') +
  annotate('text', x = 2.72, y = 1.05, label = 'B') +
  annotate('text', x = 2.9, y = 1.05, label = 'B') +
  annotate('text', x = 3.095, y = 1.05, label = 'CAD') +
  annotate('text', x = 3.28, y = 1.05, label = 'B') +
  annotate('text', x = 3.72, y = 1.05, label = 'D') +
  annotate('text', x = 3.9, y = 1.05, label = 'D') +
  annotate('text', x = 4.095, y = 1.05, label = 'D') +
  annotate('text', x = 4.28, y = 1.05, label = 'CAB') +
  annotate('text', x = 1, y = 1.11, label = paste('p =',signif(kruskal.test(norm~simple_description,data[simple_category=='all_same',])$p.value,3)),size=6) +
  annotate('text', x = 2, y = 1.11, label = paste('p =',signif(kruskal.test(norm~simple_description,data[simple_category=='A_diff',])$p.value,3)),size=6) +
  annotate('text', x = 3, y = 1.11, label = paste('p =',signif(kruskal.test(norm~simple_description,data[simple_category=='B_diff',])$p.value,3)),size=6) +
  annotate('text', x = 4, y = 1.11, label = paste('p =',signif(kruskal.test(norm~simple_description,data[simple_category=='D_diff',])$p.value,3)),size=6)
  
all_cat_boxplot_norm
ggsave('Figures/TFBS_figure/overall_specific_cats_boxplot.png',all_cat_boxplot_norm, height = 7.5, width = 10)


no_variation <- variation_data[variation == 'no_variation',]
no_variation_overall_numbers <- no_variation[,.(n = .N),by=.(balance)]
wilcox.test(prop_variation~balance,no_variation)
variation_boxplot_norm <- ggplot(no_variation, aes(x = balance, y = prop_variation, fill=balance)) +
  geom_boxplot(alpha = 0.5,linetype='dashed') +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),alpha=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..),alpha=0.5,width=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..),alpha=0.5,width=0.5) +
  geom_segment(aes(x = 1, xend = 1, y = 0.7, yend = 0.71)) +
  geom_segment(aes(x = 2, xend = 2, y = 0.7, yend = 0.71)) +
  geom_segment(aes(x = 1, xend = 2, y = 0.71, yend = 0.71)) +
  annotate('text', x = 1.5, y = 0.73, label = paste('p =',signif(wilcox.test(prop_variation~balance,no_variation)$p.value,digits=3)), size=7.5) +
  ylab('Proportion of motifs conserved across homeologs') +
  xlab('') +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(size = 20)) +
  scale_fill_brewer(palette='Set1')

ggsave('Figures/TFBS_figure/balance_variation_boxplot.png', variation_boxplot_norm, height = 10, width = 7.5)
kruskal.test(prop_variation~simple_description,no_variation)$p.value
multiple_tests <- dunnTest(prop_variation~simple_description, no_variation, method='bonferroni')$res
variation_subgenome_plot <- ggplot(no_variation, aes(x = simple_description, y=prop_variation, fill = simple_description)) +
  geom_boxplot(alpha = 0.5,linetype='dashed') +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),alpha=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..),alpha=0.5,width=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..),alpha=0.5,width=0.5) +
  geom_segment(aes(x=1,xend=1,y=0.69,yend=0.7)) +
  geom_segment(aes(x=2,xend=2,y=0.69,yend=0.7)) +
  geom_segment(aes(x=1,xend=2,y=0.7,yend=0.7)) +
  geom_segment(aes(x=1,xend=1,y=0.79,yend=0.8)) +
  geom_segment(aes(x=3,xend=3,y=0.79,yend=0.8)) +
  geom_segment(aes(x=1,xend=3,y=0.8,yend=0.8)) +
  geom_segment(aes(x=1,xend=1,y=0.89,yend=0.9)) +
  geom_segment(aes(x=4,xend=4,y=0.89,yend=0.9)) +
  geom_segment(aes(x=1,xend=4,y=0.9,yend=0.9)) +
  annotate('text',x=1.5,y=0.72,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'A - Central',]$P.adj,digits=3)),size=7.5) +
  annotate('text',x=2,y=0.82,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'B - Central',]$P.adj,digits=3)),size=7.5) +
  annotate('text',x=2.5,y=0.92,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'Central - D',]$P.adj,digits=3)),size=7.5) +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(size = 20)) +
  ylab('Proportion of motifs conserved across homoeologs') +
  xlab('Subgenome showing differential expression') +
  scale_fill_brewer(palette='Set1')
ggsave('Figures/TFBS_figure/subgenome_variation_boxplot.png', variation_subgenome_plot, height = 10, width = 7.5)

kruskal.test(prop_variation~simple_description, no_variation)
dunnTest(prop_variation~simple_description, no_variation, method='bonferroni')
no_variation_subgenomes_only <- no_variation[simple_description %in% c('A','B','D')]
kruskal.test(prop_variation~simple_description, no_variation_subgenomes_only)


extreme <- fread('Data/TFBS_search/Used/Produced/extreme_no_stress_triads_under_0.05.txt')
extreme_no_variation <- no_variation[triad %in% extreme$group_id,]

multiple_tests <- dunnTest(prop_variation~simple_description, extreme_no_variation, method='bonferroni')$res
variation_subgenome_plot <- ggplot(extreme_no_variation, aes(x = simple_description, y=prop_variation, fill = simple_description)) +
  geom_boxplot(alpha = 0.5,linetype='dashed') +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),alpha=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..),alpha=0.5,width=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..),alpha=0.5,width=0.5) +
  geom_segment(aes(x=1,xend=1,y=0.69,yend=0.7)) +
  geom_segment(aes(x=2,xend=2,y=0.69,yend=0.7)) +
  geom_segment(aes(x=1,xend=2,y=0.7,yend=0.7)) +
  geom_segment(aes(x=1,xend=1,y=0.79,yend=0.8)) +
  geom_segment(aes(x=3,xend=3,y=0.79,yend=0.8)) +
  geom_segment(aes(x=1,xend=3,y=0.8,yend=0.8)) +
  geom_segment(aes(x=1,xend=1,y=0.89,yend=0.9)) +
  geom_segment(aes(x=4,xend=4,y=0.89,yend=0.9)) +
  geom_segment(aes(x=1,xend=4,y=0.9,yend=0.9)) +
  annotate('text',x=1.5,y=0.72,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'A - Central',]$P.adj,digits=3)),size=7.5) +
  annotate('text',x=2,y=0.82,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'B - Central',]$P.adj,digits=3)),size=7.5) +
  annotate('text',x=2.5,y=0.92,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'Central - D',]$P.adj,digits=3)),size=7.5) +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(size = 20)) +
  ylab('Proportion of motifs conserved across homoeologs') +
  xlab('Subgenome showing differential expression') +
  scale_fill_brewer(palette = 'Set1')
ggsave('Figures/TFBS_figure/extreme_subgenome_variation_boxplot.png', variation_subgenome_plot, height = 10, width = 7.5)

kruskal.test(prop_variation~simple_description, extreme_no_variation)
dunnTest(prop_variation~simple_description, extreme_no_variation, method='bonferroni')

extreme_data <- data[triad %in% extreme$group_id,]
for (i in unique(extreme_data$simple_category)){
  print(i)
  filtered <- extreme_data[simple_category == i,]
  print(kruskal.test(norm~simple_description,filtered))
  print(dunnTest(norm~simple_description,filtered,method='bonferroni'))
}
extreme_all_cat_boxplot_norm <- ggplot(extreme_data, aes(y=norm, x = simple_category, fill = simple_description)) +
  geom_boxplot(alpha = 0.5,linetype='dashed') +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),alpha=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..),alpha=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..),alpha=0.5) +
  theme_classic() +
  ylab('Proportion of Motifs in Triad') +
  xlab('Motif Simple Category') +
  scale_x_discrete(labels = c('All Same', 'A Different', 'B Different', 'D Different')) +
  scale_fill_brewer(palette='Set1',name = "Triad Category", labels = c('Central', 'A dominant/ \nA suppressed', 
                                                                       'B dominant/ \nB suppressed', 'D dominant/ \nD suppressed')) +
  theme(legend.position = 'top',
        legend.key.size = unit(0.8, 'cm'),
        text = element_text(size=20)) +
  annotate('text', x = 0.72, y = 1.05, label = 'ABD') +
  annotate('text', x = 0.9, y = 1.05, label = 'C') +
  annotate('text', x = 1.095, y = 1.05, label = 'C') +
  annotate('text', x = 1.28, y = 1.05, label = 'C') +
  annotate('text', x = 2.72, y = 1.05, label = 'B') +
  annotate('text', x = 2.9, y = 1.05, label = 'B') +
  annotate('text', x = 3.095, y = 1.05, label = 'CAD') +
  annotate('text', x = 3.28, y = 1.05, label = 'B') +
  annotate('text', x = 3.72, y = 1.05, label = 'D') +
  annotate('text', x = 4.28, y = 1.05, label = 'C') +
  annotate('text', x = 1, y = 1.11, label = paste('p =',signif(kruskal.test(norm~simple_description,extreme_data[simple_category=='all_same',])$p.value,3)),size=6) +
  annotate('text', x = 2, y = 1.11, label = paste('p =',signif(kruskal.test(norm~simple_description,extreme_data[simple_category=='A_diff',])$p.value,3)),size=6) +
  annotate('text', x = 3, y = 1.11, label = paste('p =',signif(kruskal.test(norm~simple_description,extreme_data[simple_category=='B_diff',])$p.value,3)),size=6) +
  annotate('text', x = 4, y = 1.11, label = paste('p =',signif(kruskal.test(norm~simple_description,extreme_data[simple_category=='D_diff',])$p.value,3)),size=6)

ggsave('Figures/TFBS_figure/Extreme_overall_boxplot.png',extreme_all_cat_boxplot_norm, height=7.5,width=10)
### family plot
no_variation_family <- family_variation_data[variation == 'no_variation',]
family_plot <- ggplot(no_variation_family, aes(x = simple_description, y = prop_variation)) +
  geom_boxplot() +
  facet_wrap(~family)
family_plot


for(i in c('AP2; ERF', 'NAC; NAM', 'Myb/SANT', 'bZIP')){
  dat <- no_variation_family[family == i]
  multiple_tests <- dunnTest(prop_variation~simple_description, dat, method='bonferroni')$res
  variation_subgenome_plot_family <- ggplot(dat, aes(x = simple_description, y=prop_variation, fill = simple_description)) +
    geom_boxplot(alpha = 0.5,linetype='dashed') +
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),alpha=0.5) +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..),alpha=0.5,width=0.5) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..),alpha=0.5,width=0.5) +
    geom_segment(aes(x=1,xend=1,y=1.09,yend=1.1)) +
    geom_segment(aes(x=2,xend=2,y=1.09,yend=1.1)) +
    geom_segment(aes(x=1,xend=2,y=1.1,yend=1.1)) +
    geom_segment(aes(x=1,xend=1,y=1.19,yend=1.2)) +
    geom_segment(aes(x=3,xend=3,y=1.19,yend=1.2)) +
    geom_segment(aes(x=1,xend=3,y=1.2,yend=1.2)) +
    geom_segment(aes(x=1,xend=1,y=1.29,yend=1.3)) +
    geom_segment(aes(x=4,xend=4,y=1.29,yend=1.3)) +
    geom_segment(aes(x=1,xend=4,y=1.3,yend=1.3)) +
    annotate('text',x=1.5,y=1.12,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'A - Central',]$P.adj,digits=3)),size=6) +
    annotate('text',x=2,y=1.22,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'B - Central',]$P.adj,digits=3)),size=6) +
    annotate('text',x=2.5,y=1.32,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'Central - D',]$P.adj,digits=3)),size=6) +
    theme_classic() +
    theme(legend.position = 'none',
          text = element_text(size = 20)) +
    ylab('Proportion of motifs conserved across homoeologs') +
    xlab('Subgenome showing differential expression') +
    labs(title=i) +
    scale_fill_brewer(palette = 'Set1')
  ggsave(paste('Figures/TFBS_figure/Family/variation_subgenome_plot_family_',gsub('/','_',i),'.png',sep=''), variation_subgenome_plot_family, height = 10, width = 7.5)
  print(i)
  print(kruskal.test(prop_variation~simple_description, dat))
  print(dunnTest(prop_variation~simple_description, dat, method='bonferroni'))
}



### random sequences analysis 

random_data <- fread('Data/TFBS_search/Used/Produced/PlantPAN/random_seq_full_data.csv')
random_variation_data <- fread('Data/TFBS_search/Used/Produced/PlantPAN/random_seq_full_variation_data.csv')

random_data$simple_category <- factor(random_data$simple_category, levels = motif_categories)
random_data$description <- factor(random_data$description, levels = expression_categories)
random_data$balance <- factor(random_data$balance, levels = balance_categories)
random_data$simple_description <- factor(random_data$simple_description, levels = simple_description_categories)

random_variation_data$simple_category <- factor(random_variation_data$simple_category, levels = motif_categories)
random_variation_data$description <- factor(random_variation_data$description, levels = expression_categories)
random_variation_data$balance <- factor(random_variation_data$balance, levels = balance_categories)
random_variation_data$simple_description <- factor(random_variation_data$simple_description, levels = simple_description_categories)


for (i in unique(random_data$simple_category)){
  print(i)
  filtered <- random_data[simple_category == i,]
  print(kruskal.test(norm~simple_description,filtered))
  print(dunnTest(norm~simple_description,filtered,method='bonferroni'))
}

random_no_variation <- random_variation_data[variation == 'no_variation',]
random_no_variation_overall_numbers <- random_no_variation[,.(n = .N),by=.(balance)]
wilcox.test(prop_variation~balance,random_no_variation)

random_no_variation$simple_category <- factor(random_no_variation$simple_category, levels = motif_categories)
random_no_variation$description <- factor(random_no_variation$description, levels = expression_categories)
random_no_variationa$balance <- factor(random_no_variation$balance, levels = balance_categories)
random_no_variation$simple_description <- factor(random_no_variation$simple_description, levels = simple_description_categories)
kruskal.test(prop_variation~simple_description, random_no_variation)
dunnTest(prop_variation~simple_description, random_no_variation, method='bonferroni')

random_all_cat_boxplot_norm <- ggplot(data = random_data, aes(y=norm, x = simple_category, fill = simple_description)) +
  geom_boxplot(alpha = 0.5,linetype='dashed') +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),alpha=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..),alpha=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..),alpha=0.5) +
  theme_classic() +
  ylab('Proportion of Motifs in Triad') +
  xlab('Motif Simple Category') +
  scale_x_discrete(labels = c('All Same', 'A Different', 'B Different', 'D Different')) +
  scale_fill_brewer(palette='Set1',name = "Triad Category", labels = c('Central', 'A dominant/ \nA suppressed', 
                                                                       'B dominant/ \nB suppressed', 'D dominant/ \nD suppressed')) +
  theme(legend.position = 'top',
        legend.key.size = unit(0.8, 'cm'),
        text = element_text(size=20)) +
  annotate('text', x = 0.72, y = 1.05, label = 'ABD') +
  annotate('text', x = 0.9, y = 1.05, label = 'C') +
  annotate('text', x = 1.095, y = 1.05, label = 'C') +
  annotate('text', x = 1.28, y = 1.05, label = 'C') +
  annotate('text', x = 1.72, y = 1.05, label = 'A') +
  annotate('text', x = 1.9, y = 1.05, label = 'CBD') +
  annotate('text', x = 2.095, y = 1.05, label = 'A') +
  annotate('text', x = 2.28, y = 1.05, label = 'A') +
  annotate('text', x = 2.72, y = 1.05, label = 'B') +
  annotate('text', x = 2.9, y = 1.05, label = 'B') +
  annotate('text', x = 3.095, y = 1.05, label = 'CAD') +
  annotate('text', x = 3.28, y = 1.05, label = 'B') +
  annotate('text', x = 3.72, y = 1.05, label = 'D') +
  annotate('text', x = 3.9, y = 1.05, label = 'D') +
  annotate('text', x = 4.095, y = 1.05, label = 'D') +
  annotate('text', x = 4.28, y = 1.05, label = 'CAB') +
  annotate('text', x = 1, y = 1.11, label = paste('p =',signif(kruskal.test(norm~simple_description,random_data[simple_category=='all_same',])$p.value,3)),size=6) +
  annotate('text', x = 2, y = 1.11, label = paste('p =',signif(kruskal.test(norm~simple_description,random_data[simple_category=='A_diff',])$p.value,3)),size=6) +
  annotate('text', x = 3, y = 1.11, label = paste('p =',signif(kruskal.test(norm~simple_description,random_data[simple_category=='B_diff',])$p.value,3)),size=6) +
  annotate('text', x = 4, y = 1.11, label = paste('p =',signif(kruskal.test(norm~simple_description,random_data[simple_category=='D_diff',])$p.value,3)),size=6)

random_all_cat_boxplot_norm


random_variation_boxplot_norm <- ggplot(random_no_variation, aes(x = balance, y = prop_variation, fill=balance)) +
  geom_boxplot(alpha = 0.5,linetype='dashed') +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),alpha=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..),alpha=0.5,width=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..),alpha=0.5,width=0.5) +
  geom_segment(aes(x = 1, xend = 1, y = 0.7, yend = 0.71)) +
  geom_segment(aes(x = 2, xend = 2, y = 0.7, yend = 0.71)) +
  geom_segment(aes(x = 1, xend = 2, y = 0.71, yend = 0.71)) +
  annotate('text', x = 1.5, y = 0.73, label = paste('p =',signif(wilcox.test(prop_variation~balance,random_no_variation)$p.value,digits=3)), size=7.5) +
  ylab('Proportion of motifs conserved across homeologs') +
  xlab('') +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(size = 20)) +
  scale_fill_brewer(palette='Set1')
random_variation_boxplot_norm

multiple_tests <- dunnTest(prop_variation~simple_description, random_no_variation, method='bonferroni')$res
random_variation_subgenome_plot <- ggplot(random_no_variation, aes(x = simple_description, y=prop_variation, fill = simple_description)) +
  geom_boxplot(alpha = 0.5,linetype='dashed') +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),alpha=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..),alpha=0.5,width=0.5) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..),alpha=0.5,width=0.5) +
  geom_segment(aes(x=1,xend=1,y=1.02,yend=1.03)) +
  geom_segment(aes(x=2,xend=2,y=1.02,yend=1.03)) +
  geom_segment(aes(x=1,xend=2,y=1.03,yend=1.03)) +
  geom_segment(aes(x=1,xend=1,y=1.12,yend=1.13)) +
  geom_segment(aes(x=3,xend=3,y=1.12,yend=1.13)) +
  geom_segment(aes(x=1,xend=3,y=1.13,yend=1.13)) +
  geom_segment(aes(x=1,xend=1,y=1.22,yend=1.23)) +
  geom_segment(aes(x=4,xend=4,y=1.22,yend=1.23)) +
  geom_segment(aes(x=1,xend=4,y=1.23,yend=1.23)) +
  annotate('text',x=1.5,y=1.055,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'A - Central',]$P.adj,digits=3)),size=7.5) +
  annotate('text',x=2,y=1.155,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'B - Central',]$P.adj,digits=3)),size=7.5) +
  annotate('text',x=2.5,y=1.255,label=paste('p =',signif(multiple_tests[multiple_tests$Comparison == 'Central - D',]$P.adj,digits=3)),size=7.5) +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(size = 20)) +
  ylab('Proportion of motifs conserved across homoeologs') +
  xlab('Subgenome showing differential expression') +
  scale_fill_brewer(palette='Set1')

random_variation_subgenome_plot

ggsave('Figures/Report/Appendix/random_seq_subgenome.png',random_variation_subgenome_plot,height=7.5,width=10)
