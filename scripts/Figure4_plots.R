library(ggplot2)
library(data.table)
library(magrittr)
library(forcats)
library(tidyr)
library(dplyr)
library(directlabels)
library(patchwork)
library(ggrepel)
library(ggpubr)
#library(ggpmisc)
library(scales)
library(stringr)

df = fread('invivo_metagenomes_FinalPolyTable.tsv', sep='\t')



#### run ep_figure_invivo.py before
lt = fread('tmp/LatePrevalentPolys.csv')
ea = fread('tmp/EarlyPersistentPolys.csv')
idi = fread("tmp/IndividualPersistentPolys.csv")
para = fread("tmp/ParallelPolys.csv")
para_polys <- para %>% pull(PolyID)
lt_polys <- lt %>% pull(PolyID)
ea_polys <- ea %>% pull(PolyID)
idi_polys <- idi %>% pull(PolyID)
#df <- df %>% mutate(Gene_pNpS = case_when(Gene_pNpS == -1 ~ 0, TRUE ~ Gene_pNpS))
ms = fread("tmp/confirmed_invivo.csv")
keep_samples <- df %>% filter(Chrom == 'NC_004663' & Coverage_mean >= 100) %>% distinct(SampleID) %>% pull(SampleID)
df <- df %>% filter(!PolyID %in% ms$PolyID & SampleType == "Feces" & TotalCov >= 100 & SampleID %in% keep_samples)%>% distinct(SampleID,PolyID, .keep_all = T)
#df <- df %>% filter(!PolyID %in% ms$PolyID & SampleType == "Feces" & Coverage_mean >= 100) %>% distinct(SampleID,PolyID, .keep_all = T)
#df <- df %>% mutate(Label=paste(GeneName," (",oldGeneName,")", sep=''))
df <- df %>% mutate(Label = case_when(GeneName != 'Intergenic' & oldGeneName =='Nan' ~ GeneName,
									  GeneName == 'Intergenic'~ as.character(Pos),
									  .default = oldGeneName))

total_samples <- length(df %>% distinct(SampleID) %>% pull(SampleID))
total_samples
#### Look at dataset prevalence and persistence
df <- df %>% group_by(PolyID) %>% mutate(Prevalence = n()) %>% ungroup()
df <- df %>% group_by(PolyID,Mouse) %>% mutate(MousePrevalence = n()) %>% ungroup()
df <- df %>% mutate(PrevalencePercent = Prevalence / total_samples)
df <- df %>% mutate(NotRandomGroup = case_when((MousePrevalence == 4 & PrevalencePercent <= 0.35) | (MousePrevalence == 5 & PrevalencePercent <= 0.55)  ~ 'universal', .default = 'random'))
df <- df %>% mutate(PersistentGroup = ifelse(MousePrevalence >= 4,'yes', 'no'))
df <- df %>% mutate(PrevalentGroup = ifelse(PrevalencePercent > 0.2,'yes', 'no'))
df <- df %>% mutate(TargetGroup = case_when(PolyID %in% lt_polys ~ 'LR',
                                            PolyID %in% ea_polys & !PolyID %in% lt_polys   ~ 'ER',
											.default = 'SP'))


target_genes <- df %>% filter(TargetGroup != 'SP'
								   & Ancestry == 'evolved' & GeneName != 'Intergenic') %>%
	distinct(GeneName) %>% pull(GeneName)

f1 <- ggplot(df %>% filter(oldGeneName == 'BT2172'), aes(x=TargetGroup,y=Freq, fill=TargetGroup)) + geom_boxplot(width=0.5) + theme_bw() +
	stat_compare_means(label='p.signif', size= 4, vjust = 0.5,hjust=-0.5, comparisons = list(c('ER', 'SP')))+
	theme(plot.margin = margin(0.5,0.5,0.1,0.1),legend.position = 'none', aspect.ratio = 2,
		  axis.title = element_text(size=14), axis.text.y = element_text(size=12),
		  axis.text.x = element_text(size=14))+ labs(y='Frequency', x='', title = 'BT2172') +
	scale_fill_manual(values = c('SP' = 'darkgrey', 'ER'='darkred',
								  'LR'='orange', "Parallel"='darkmagenta','Individual persistent'='darkgrey')) +
	coord_cartesian(clip='off')
f1
compare_means(Freq ~ TargetGroup, data=df %>% filter(oldGeneName == 'BT0733'))

f2 <- ggplot(df %>% filter(oldGeneName == 'BT2469'), aes(x=TargetGroup,y=Freq, fill=TargetGroup)) + geom_boxplot(width=0.5) + theme_bw() +
	stat_compare_means(label='p.signif', size= 4, vjust = 0.5,hjust=-0.5, comparisons = list(c('ER', 'SP')))+
	theme(plot.margin = margin(0.5,0.5,0.1,0.1),legend.position = 'none',aspect.ratio = 2,
		  axis.title = element_text(size=14), axis.text.y = element_text(size=12),
		  axis.text.x = element_text(size=14)) + ggtitle('BT2469')+ labs(y='Frequency', x='') +
	scale_fill_manual(values = c('SP' = 'darkgrey', 'ER'='darkred',
								  'LR'='orange', "Parallel"='darkmagenta',
								 'Individual persistent'='darkgrey')) +
	coord_cartesian(clip='off')
f2


f3 <- ggplot(df %>% filter(oldGeneName == 'BT0867'), aes(x=TargetGroup,y=Freq, fill=TargetGroup)) + geom_boxplot(width=0.5) + theme_bw() +
	stat_compare_means(comparisons = list(c('LR', 'SP')),label='p.signif',size= 4, vjust = 0.5,hjust=-0.5) +
	theme(plot.margin = margin(0.5,0.5,0.1,0.1),legend.position = 'none', aspect.ratio = 2,
		  axis.title = element_text(size=14), axis.text.y = element_text(size=12), axis.text.x= element_text(size=14)) + ggtitle('BT0867') + labs(y='Frequency', x='') +
	scale_fill_manual(values = c('SP' = 'darkgrey', 'ER'='darkred',
								  'LR'='orange', "Parallel"='darkmagenta','Individual persistent'='darkgrey')) +
	coord_cartesian(clip = 'off')
f3

compare_means(Freq ~ TargetGroup, data=df %>% filter(oldGeneName == 'BT4029'))


f4 <- ggplot(df %>% filter(oldGeneName == 'BT3621'), aes(x=TargetGroup,y=Freq, fill=TargetGroup)) + geom_boxplot(width=0.5) + theme_bw() +
	stat_compare_means(comparisons = list(c('ER', 'SP')),label='p.signif', size= 4,
					   vjust = 0.5,hjust=-0.5) + theme(plot.margin = margin(0.5,0.5,0.1,0.1),legend.position = 'none', aspect.ratio = 2,
													   axis.title = element_text(size=14), axis.text.y = element_text(size=12),axis.text.x = element_text(size=14))+ ggtitle('BT3621') + labs(y='Frequency', x='') +
	scale_fill_manual(values = c('SP' = 'darkgrey', 'ER'='darkred',
								  'LR'='orange', "Parallel"='darkmagenta','Individual persistent'='darkgrey'))  +
	coord_cartesian(clip = 'off')
f4
f5 <- ggplot(df %>% filter(oldGeneName == 'BT3368'), aes(x=TargetGroup,y=Freq, fill=TargetGroup)) + geom_boxplot(width=0.5) + theme_bw() +
	stat_compare_means(comparisons = list(c('LR', 'SP')),label='p.signif', size= 4,
					   vjust = 0.5,hjust=-0.5) + theme(plot.margin = margin(0.5,0.5,.1,0.1),legend.position = 'none',
													   aspect.ratio = 2, axis.title = element_text(size=14), axis.text.y = element_text(size=12),axis.text.x = element_text(size=14))+ ggtitle('BT3368') + labs(y='Frequency', x='') +
	scale_fill_manual(values = c('SP' = 'darkgrey', 'ER'='darkred',
								  'LR'='orange', "Parallel"='darkmagenta','Individual persistent'='darkgrey')) +
	coord_cartesian(clip = 'off')
f5

freq <- f1 + f2 + f4 + f3 + f5
freq
## get the persistent or very prevalent genes or polys
fig3b <- df %>% filter(GeneName %in% target_genes & Ancestry == 'evolved') %>%
	arrange(PersistentGroup) %>% distinct(PolyID, .keep_all = T) %>%
	ggplot(aes(fct_reorder(Label, PrevalencePercent, .fun = max), y=PrevalencePercent)) + theme_bw() +
	geom_violin() + geom_jitter(size=1,alpha=0.6, aes(color=TargetGroup)) + coord_flip(clip='off')  + guides(colour = guide_legend(override.aes = list(size=3))) +
	theme(plot.margin = margin(0,0,0.5,0), aspect.ratio = 2, legend.position = c(0.75,0.15), legend.box.background = element_rect(colour = "black"), legend.text = element_text(size=12),axis.title = element_text(size=14), axis.text = element_text(size=10)) +
	labs(x='Gene', y='Prevalence') +
	scale_color_manual(name= '',values = c('SP' = 'darkgrey', 'ER'='darkred',
								  'LR'='orange', "Parallel"='darkmagenta','Individual persistent'='darkgrey'))

fig3b
ggsave('Fig3_Prevalent_polys_target_groups.png')
ggsave('Fig3_Prevalent_polys_target_groups.svg', width = 10, height = 10)





gn = fread("tmp/gene_enrichment.csv")
group_color_map_gene = c('Enriched'='#d76650', 'no'='#2065ab')

### to be considered enriched Score > 1 and at least 2 polys present
gn <- gn %>% mutate(Group = case_when(Score > 1 & PassPolyCutoff == 'yes' ~ "Enriched", .default = 'no'))
gn <- gn %>% filter(SampleType == "Feces")
gn <- gn %>% mutate(GroupBinary = case_when(Group == 'Enriched'~1, Group == 'no' ~ 0))
gn <- gn %>% mutate(Label = case_when(oldGeneName =='Nan' ~ GeneName, .default = oldGeneName))
#gn <- gn %>% mutate(Label=paste(GeneName," (",oldGeneName,")", sep=''))
gn <- gn %>% filter(SampleType == "Feces")
highscoregene <- gn %>% filter(Group == 'Enriched') %>% group_by(GeneName) %>% summarize(N = n()) %>%
	filter(N >= 10) %>% arrange(N) %>%  pull(GeneName)
gn <- gn %>% filter(GeneName %in% highscoregene)

unique(gn$Label)

gn <- gn %>% group_by(Label) %>% mutate(EnrichedPrc = sum(GroupBinary/total_samples))

gncp <- ggplot(gn %>% distinct(GeneName,.keep_all = T), aes(fct_reorder(Label,EnrichedPrc), EnrichedPrc)) + theme_bw() +
	geom_bar(stat='identity', na.rm = T, fill='#d76650', alpha=0.6) + coord_flip(clip="off") +
	theme(axis.title = element_text(size=14),panel.grid.minor = element_blank(), aspect.ratio = 5,
		  axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size=10),
		  plot.margin = margin(0,5,0.5,0)) +
	scale_x_discrete(na.translate = F) +
	labs(x='',y='Prevalence')
gncp
gep <- gn %>% ggplot(aes(x=fct_reorder(Label,EnrichedPrc), y=log10(Score))) +
  geom_boxplot(outlier.size = .1, outlier.shape = NA) + geom_jitter(aes(color=Group),size=1, width=.1, alpha=0.4)+
  scale_color_manual(values =group_color_map_gene) + labs(y=expression('log'[10]*' Enrichment score' ), x= 'Gene') + theme_bw()+
  #scale_color_gradient2(midpoint = 1, high=muted("red"), low=muted("blue"), mid = 'lightgrey') +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=9),panel.grid.minor = element_blank(),
		aspect.ratio = 3, legend.position = 'none', plot.margin = margin(0,0,0.5,0)) + coord_flip(clip='off')


### This plot shows genes that have been enriched in polymorphisms in at least 10 samples in the dataset.
### Enrichement score is estimated and a gene is considered enriched if at least 5 (or 10) polymorphisms are detected in a given sample.
### see ep_cog.py for more details
fig3e <- gep + gncp
fig3e
ggsave('/Fig4_Enriched_genes.png')
ggsave('Fig4_Enriched_genes.svg', height=10, width=10)

design <- "
AABCC
DEFGH
"


(fig3e + fig3b ) / (f1 + f2 + f4 + f3 + f5 + plot_layout(nrow=1)) + plot_layout(heights = c(3,1)) + plot_annotation(tag_levels = 'A')# + theme(plot.margin = margin(.5,.5,.5,.5))

ggsave("Fig4.svg", width = 790, height = 1120*0.6, units = 'px', dpi = 'screen')

gep | gncp | fig3b #+ plot_layout( nrow =1, ncol = 3, heights = c(1,1,1), widths = c(1,1,1))
ggsave("Fig4_p1.png")
ggsave("Fig4_p1.svg")
f1 + f2 + f3 + f4 + f5 + plot_layout(nrow=1)
ggsave("Fig4_p2.svg")


colnames(df)



df <- df %>% group_by(SampleID, GeneName) %>% mutate(GenePolyCount = n()) %>% ungroup()



target_genes
print(df %>% filter(GeneName %in% target_genes) %>% distinct(Label), n=30)


for (i in target_genes){
	print(i)
}
as.character(target_genes)


#### show genes that are associated with persistent polys and have mean gene_pNpS > 0
mean_positive_gene_pNpS <- df %>% group_by(GeneName) %>% mutate(Mean = mean(Gene_pNpS)) %>%
	filter(Mean > 1) %>% distinct(GeneName) %>% pull(GeneName)

mean_positive_gene_pNpS


hot_genes <- df %>% filter(Ancestry == 'evolved' & GeneName != 'Intergenic') %>% group_by(GeneName) %>% summarise(median = median(GenePolyCount, na.rm = TRUE)) %>% arrange(desc(median)) %>%
	head(20) %>% distinct(GeneName) %>% pull(GeneName)



df %>% filter(Ancestry == 'evolved' & GeneName %in% hot_genes) %>%
	ggplot(aes(x=fct_reorder(Label, GenePolyCount, .fun = median), y=GenePolyCount)) +
	geom_violin() + geom_jitter(alpha=0.2, aes(color = TargetGroup)) + coord_flip() +
	scale_color_manual(values = c('Sporadic' = 'darkgrey', "Early persistent"='darkred','Late global'='darkmagenta'))


df %>% filter(Ancestry == 'evolved' & str_detect(GeneName, 'BT_')) %>%
	distinct(SampleID, GeneName, .keep_all = T) %>% arrange(HighpnpsGroup) %>%
	ggplot(aes(y=GenePolyCount, x=Gene_pNpS, color=HighpnpsGroup)) + geom_point(alpha=0.2)+
	scale_color_manual(values = c('no' = 'darkgrey', 'yes'='darkmagenta'))

df %>% filter(Ancestry == 'evolved' & str_detect(GeneName, 'BT_')) %>%
	distinct(SampleID, GeneName, .keep_all = T) %>% arrange(desc(TargetGroup)) %>%
	ggplot(aes(y=as.factor(TrajectoriesDetected), x=Gene_pNpS, color=TargetGroup)) + geom_point(alpha=0.2)+

	scale_color_manual(values = c('Sporadic' = 'darkgrey', "Early persistent"='darkred','Late global'='darkmagenta'))



aa



## enrichment
cog =fread('tmp/cog_cat_enrichment.csv')
cog <- cog %>% mutate(Group = case_when(Score > 1 ~ "greater", Score < 1 ~ 'less'))
cog <- cog %>% filter(SampleType == "Feces")
group_color_map = c('greater'='#d76650', 'equal'='lightgrey','less'='#2065ab')


## count how many times (samples) a COG category is enriched in the dataset as a percentage of all samples.
total_samples <- length(cog %>% distinct(SampleID) %>% pull(SampleID))
cog <- cog %>% mutate(GroupBinary = ifelse(Group =='greater', 1, 0))
enc <- cog %>% group_by(CogCategoryDesc) %>% summarize(EnrichedPrc = sum(GroupBinary)/total_samples) %>% arrange(EnrichedPrc)
encp <- ggplot(enc, aes(fct_reorder(CogCategoryDesc,EnrichedPrc), EnrichedPrc)) + theme_bw() +
	geom_bar(stat='identity', na.rm = T, fill='#d76650', alpha=0.6) + coord_flip() +
	theme(panel.grid.minor = element_blank(), aspect.ratio = 5, axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
	scale_x_discrete(na.translate = F) +
	labs(x='',y='Enrichment prevalence')
encp


## relevel cog object based on the order of enriched COGs from previous plot to place them side by side..
cog <- cog %>% mutate(CogCategoryDesc = fct_relevel(CogCategoryDesc, as.character(enc$CogCategoryDesc)))
#fct_reorder(CogCategory,Score, .fun = 'median', .desc = F)
cop <- ggplot(cog, aes(x=CogCategoryDesc, y=Score)) +
  geom_boxplot(width=0.5 , outlier.size = .1, outlier.shape = NA) + geom_jitter(aes(color=Group),alpha=0.4,size=1, width=.1)+
  scale_color_manual(values =group_color_map) + labs(y='Enrichment score', x= 'COG Category') + theme_bw() +
  #scale_color_gradient2(midpoint = 1, high=muted("red"), low=muted("blue"), mid = 'lightgrey') +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_text(angle=0, size=8), aspect.ratio = 2, legend.position = 'none') + coord_flip()
cop + encp


cig = fread('tmp/cog_id_enrichment.csv')
cig <- cig %>% mutate(Group = case_when(Score > 1  & Count >= 2 ~ "greater",  Score > 1 & Count < 2 ~ 'greater_uni',Score < 1 ~ 'less', .default = 'less'))

highscoreCogID <- cig %>% filter(Score > 1 & Count >= 2) %>% group_by(CogID) %>% summarize(N = n()) %>% filter(N >= 10) %>% arrange(N) %>%  pull(CogID)
total_samples <- length(cig %>% distinct(SampleID) %>% pull(SampleID))
cig <- cig %>% mutate(GroupBinary = case_when(Group =='greater' ~ 1, Group == 'greater_uni' ~ 0, Group == 'less' ~ 0))
cig <- cig %>% filter(CogID %in% highscoreCogID & SampleType == 'Feces')

enci <- cig %>% group_by(CogID) %>% summarize(EnrichedPrc = sum(GroupBinary)/total_samples) %>% arrange(EnrichedPrc)
encip <- ggplot(enci, aes(fct_reorder(CogID,EnrichedPrc), EnrichedPrc)) +
	geom_bar(stat='identity', na.rm = T, fill='#d76650', alpha=.6) + coord_flip() + theme_bw()+
	theme(panel.grid.minor = element_blank(), aspect.ratio = 5, axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
	scale_x_discrete(na.translate = F) +
	labs(x='',y='Enrichment prevalence')
encip
cig <- cig %>% mutate(CogID = fct_relevel(CogID, as.character(enci$CogID)))

cip <- ggplot(cig, aes(x=CogID, y=Score)) +
  geom_boxplot(outlier.size = .1, outlier.shape = NA) + geom_jitter(aes(color=Group),size=1, width=.1, alpha=0.4)+
  scale_color_manual(values =group_color_map) + labs(y='Enrichment score', x= 'COG ID') +
	theme_bw()+
  #scale_color_gradient2(midpoint = 1, high=muted("red"), low=muted("blue"), mid = 'lightgrey') +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_text(angle=0, size=10), aspect.ratio = 2, legend.position = 'none') + coord_flip()

cip + encip



design = "AABCCC
          AABCCC
          DEFGH#
          "


(fig3e + fig3b) / (freq + plot_layout(nrow=1)) + plot_layout(heights = c(2,1))

ggsave("Fig4.svg", width=14, height=12)


cog_plot = cop + encp
cog_plot + cip + encip + gep + gncp + plot_layout(design = 'ABCDEF')
