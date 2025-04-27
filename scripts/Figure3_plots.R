library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(forcats)
library(ggpubr)
library(vegan)
library(ggExtra)
library(ggpointdensity)
library(viridis)
library(MASS)


#setwd('')

df = fread('invivo_metagenomes_FinalPolyTable.tsv', sep='\t')
mask = fread('tmp/confirmed_invivo.csv')
mask <- mask %>% distinct(PolyID)
mask_polys <- mask$PolyID

## filter data
df <- df %>% filter(!PolyID %in% mask_polys) #### mask polymorphisms within shuffling regions
keep_samples <- df %>% filter(Chrom == 'NC_004663' & Coverage_mean >= 100) %>% distinct(SampleID) %>% pull(SampleID)
df <- df %>% filter(TotalCov >= 100 & SampleID %in% keep_samples)

mf <- df %>% filter(SampleProfile == 'metagenome' & SampleType == "Feces") %>% distinct(SampleID, PolyID, .keep_all = T)

### estimate fraction of the ancestral polys in the experiment
mf[Dataset=='bern_metagenomes' & Day == 3] %>% group_by(Ancestry) %>% distinct(PolyID, .keep_all = T) %>% summarise(Count = n()) %>% mutate(fraction = Count / sum(Count))

###### FIGURE 1A - Poly counts by Type over time ##########w
pc <- mf %>% filter(Type %in% c("SNP","DEL",'INS')) %>%
  count(Type,SampleID, Mouse, Day, Cage, Dataset)
pc <- pc %>% mutate(Type =fct_relevel(Type, "SNP", "DEL", "INS"))
pc$Day = as.factor(pc$Day)
pc$Mouse = as.factor(pc$Mouse)
pc$Cage = as.factor(pc$Cage)
pc$Dataset = as.factor(pc$Dataset)

pc2 <- mf %>% filter(Type %in% c("SNP","DEL","INS")) %>% count(Type,SampleID,Mouse,Day,Cage)



ta1 <- ggplot(pc %>% filter(Type == 'SNP'), aes(as.factor(Day), log10(n))) +
	geom_boxplot(outlier.shape = NA) + geom_point(alpha=0.5, size=1) + geom_smooth(aes(group=NA)) +
	labs(x="Time (Days)", y=expression("log"[10]*" Count"), title = 'SNP') + theme_bw() + theme(axis.title = element_text(size=14), axis.text = element_text(size=12), legend.text = element_text(size=14),
         legend.title = element_text(size=14), aspect.ratio = 1, plot.title = element_text(size=14,hjust = 0.5)) + ylim(0,4)
ta1
ta2 <- ggplot(pc %>% filter(Type == 'DEL'), aes(as.factor(Day), log10(n))) +
	geom_boxplot(outlier.shape = NA) + geom_point(alpha=0.5, size=1) + geom_smooth(aes(group=NA)) +
	labs(x="Time (Days)", y=expression("log"[10]*" Count"), title='DEL') +
	theme_bw() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title = element_text(size=14), axis.text = element_text(size=12), legend.text = element_text(size=14),
         legend.title = element_text(size=14), strip.text = element_text(size=16), aspect.ratio = 1, plot.title = element_text(size=14,hjust = 0.5)) + ylim(0,4)
ta3 <- ggplot(pc %>% filter(Type == 'INS'), aes(as.factor(Day), log10(n))) +
	geom_boxplot(outlier.shape = NA) + geom_point(alpha=0.5, size=1) + geom_smooth(aes(group=NA)) +
	labs(x="Time (Days)", y=expression("log"[10]*" Count"), title='INS') +
	theme_bw() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title = element_text(size=14), axis.text = element_text(size=12), legend.text = element_text(size=14),
         legend.title = element_text(size=14), strip.text = element_text(size=16), aspect.ratio = 1, plot.title = element_text(size=14,hjust = 0.5)) + ylim(0,4)


f1a = ta1 + ta2 + ta3
f1a



##### pi diversity
pi = fread("tmp/genome_pi.csv")
## filter for samples according to previous plot
pi <- pi %>% filter(SampleID %in% keep_samples & SampleType == 'Feces')

pip <- ggplot(pi, aes(as.factor(Day), log10(pi))) + geom_boxplot(outlier.shape = NA, width=0.5) + geom_point(alpha=0.5, size=1) +
  labs(y=expression("log"[10]*"  "*pi[genes]), x='Time (Days)') + theme(aspect.ratio = 1) + geom_smooth(aes(group=NA)) +
  theme_bw() + theme(axis.title = element_text(size=14), axis.text = element_text(size=12), legend.text = element_text(size=14),
        legend.title = element_text(size=14), strip.text = element_text(size=16), aspect.ratio = 1)
pip
pi %>% group_by(Day) %>% summarise(Median=median(pi))
pi %>% filter(Day %in% c(7,14))
compare_means(pi ~ Day, data= pi %>% filter(Day %in% c(7,3)) %>% arrange(Mouse), paired = F)


##### stats
# ns
pcStats37 <- as.data.table(compare_means(n ~ Day, data = pc %>% filter(Type == "SNP" & Day %in% c(3,7)) %>% arrange(Mouse), method='wilcox.test', paired = T))
pcStats314 <- as.data.table(compare_means(n ~ Day, data = pc %>% filter(Type == "SNP" & Day %in% c(3,14))  %>%arrange(Mouse), method='wilcox.test', paired = T))
pcStats321 <- as.data.table(compare_means(n ~ Day, data = pc %>% filter(Type == "SNP" & Day %in% c(3,21)) %>%arrange(Mouse), method='wilcox.test', paired = T))
pcStats328 <- as.data.table(compare_means(n ~ Day, data = pc %>% filter(Type == "SNP" & Day %in% c(3,28)) %>%arrange(Mouse), method='wilcox.test', paired = T))
pcStats328
# ns
pcStats714 <- as.data.table(compare_means(n ~ Day, data = pc %>% filter(Type == "SNP" & Day %in% c(7,14)) %>%arrange(Mouse), method='wilcox.test', paired = T))
pcStats721 <- as.data.table(compare_means(n ~ Day, data = pc %>% filter(Type == "SNP" & Day %in% c(7,21)) %>%arrange(Mouse), method='wilcox.test', paired = T))
pcStats728 <- as.data.table(compare_means(n ~ Day, data = pc %>% filter(Type == "SNP" & Day %in% c(7,28)) %>%arrange(Mouse), method='wilcox.test', paired = T))
# ns
pcStats1421 <- as.data.table(compare_means(n ~ Day, data = pc %>% filter(Type == "SNP" & Day %in% c(21,14))%>%arrange(Mouse), method='wilcox.test', paired = T))
pcStats1428 <- as.data.table(compare_means(n ~ Day, data = pc %>% filter(Type == "SNP" & Day %in% c(28,14)) %>%arrange(Mouse), method='wilcox.test', paired = T))
pcStats1428
# ns
pcStats2128 <- as.data.table(compare_means(n ~ Day, data = pc %>% filter(Type == "SNP" & Day %in% c(21,28)) %>%arrange(Mouse), method='wilcox.test', paired = T))

# ### Significant differces between 3 - 14 and 3 - 21 pairs
#
# pcStats3 <- as.data.table(compare_means(n ~ Dataset, group.by = c('Day','Type'), data = pc %>% select(Mouse,Day,n,Dataset,Type) %>% arrange(Mouse), method='wilcox.test', paired = F))
# print(pcStats3[p.adj <= 0.05])

### PERMANOVA (distance matrix from ep_figure_pcoa.py)
dm <- read.table("tmp/DistanceMatrix.csv", header=T, sep=',', row.name="SampleID",check.names = F) ## run ep_figure_1.py before
me <- read.table("tmp/DM_metadata.csv", header=T, sep=',', row.name='SampleID',check.names = F)
me$Day <- factor(me$Day)
### PERMANOVA ( permanova clearly indicates a Dataset signal - if day 3 is removed then there is no variation explained by Day
### however Dataset seems to explain 7% of the variability in polymorphisms
me_sub <- subset(me, me$Day != 3)
dm_sub <- dm[rownames(me_sub), rownames(me_sub)]
adonis2(dm_sub ~ Day + Dataset * Cage * Mouse, data=me_sub, stata=Mouse)
adonis2(dm ~ Day + Dataset ,data=me, stata=Mouse)
adonis2(dm ~ me$Day + me$Dataset/me$Cage/me$Mouse)



tr = fread('tmp/Turnover.csv')
#tr <- tr %>% filter(Dataset == 'bern_metagenomes')
tr <- tr %>% filter(!Mouse %in% c("p1",'p4','p5','p6'))


ggplot(tr, aes(x=as.factor(Day), y=, group=Mouse )) +
	geom_smooth(aes(y=log10(Total),group=NA), color='darkred') +
	geom_smooth(aes(y=log10(Denovo), group=NA),color='darkblue')+
	geom_smooth(aes(y=log10(FirstEmerge), group=NA), color='darkgreen') + ylim(0,4)

tr <- tr %>% mutate(StandingFreq = (Total - Denovo) /Total)
tr <- tr %>% mutate(DenovoFreq = (Denovo - FirstEmerge) / Total)
tr <- tr %>% mutate(DenovoPrc = Denovo/Total)
tr <- tr %>% mutate(FirstEmergeFreq = FirstEmerge / Total)


compare_means(DenovoPrc ~ Day, data = tr)


tr %>% filter(Day > 3) %>% arrange(DenovoPrc) %>% pull(DenovoPrc)
tr %>% filter(Day > 3) %>% pull(DenovoPrc) %>% min()

ggplot(tr %>% pivot_longer(cols=c('StandingFreq','DenovoFreq','FirstEmergeFreq')), aes(x=as.factor(Day), y=value,fill=name)) + geom_bar(position = 'stack', stat='identity') + facet_wrap(~Mouse)

total <- ggplot(tr, aes(x=as.factor(Day), group=Mouse)) +
	geom_line(aes(y=log10(Total)),alpha=0.2) +
	geom_smooth(aes(y=log10(Total),group=NA), alpha=0.5) +
	ylim(0,4) + labs(title='Total',x='Time (Days)',y=expression('log'[10]*' Count')) +
	theme_bw() + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),aspect.ratio = 1,
		  axis.title = element_text(size=14),axis.text = element_text(size=12) , plot.title = element_text(size=14,hjust = 0.5))
total
b <- ggplot(tr, aes(x=as.factor(Day), group=Mouse)) +
	geom_line(aes(y=DenovoPrc),alpha=0.2, color='darkblue') +
	geom_smooth(aes(y=DenovoPrc,group=NA), color='darkblue',fill='blue', alpha=0.5) +
	ylim(0,1) + labs(title=expression("\u0394N"),x='Time (Days)',y='Total') +
	theme_bw() + theme(aspect.ratio = 1,axis.title = element_text(size=14),
		  axis.text = element_text(size=12) , plot.title = element_text(size=14,hjust = 0.5)) +
	scale_y_continuous(labels = scales::percent, limits = c(0,1), breaks = c(0,0.5,1))
c <- ggplot(tr, aes(x=as.factor(Day), group=Mouse)) +
	geom_line(aes(y=FirstEmergeFreq),alpha=0.2, color='darkgreen') +
	geom_smooth(aes(y=FirstEmergeFreq,group=NA), color='darkgreen',fill='green', alpha=0.5) +
	ylim(0,1)  + labs(title='Mouse de novo',x='Time (Days)',y='Total') +
	theme_bw() + theme(aspect.ratio = 1,axis.title = element_text(size=14),axis.text = element_text(size=12) ,
		  axis.ticks.y = element_blank(),plot.title = element_text(size=14,hjust = 0.5)) + scale_y_continuous(labels = scales::percent, limits = c(0,1), breaks=c(0,0.5,1))
b | c
c


d + e + f
cor.test(method = 'spearman', tr$Day, tr$Total )
cor.test(method = 'spearman', tr$Day, tr$Denovo )
cor.test(method = 'spearman', tr$Day, tr$FirstEmerge )
compare_means(data = tr, formula = Total ~ Day, paired = T, method = 'kruskal.test')
compare_means(data = tr, formula = Denovo ~ Day, paired = T, method = 'kruskal.test')
compare_means(data = tr, formula = FirstEmerge ~ Day, paired = T, method = 'kruskal.test')




## Synonymy by day
st <- fread('tmp/SynonymyCount.csv')
st <- st[SampleID %in% (df$SampleID)]
st <- st %>% mutate(Dataset=case_when(Dataset == 'bern_metagenomes' ~ 'Facility A', Dataset == 'ploen' ~ 'Facility B'))
st <- st %>% mutate(Synonymy=case_when(Synonymy == 'Nonsynonymous' ~ 'Non-synonymous', Synonymy == 'Synonymous'~'Synonymous',
									   Synonymy == 'undefined'~'Intergenic'))
st$Synonymy <- factor(st$Synonymy, levels = c("Non-synonymous",'Synonymous','Intergenic'))


stp <- ggplot(st, aes(as.factor(Day), log10(Count), color=Synonymy, fill=Synonymy)) +
  geom_boxplot(alpha=0.4, width=0.6, outlier.shape = NA) + geom_point(alpha=0.4,position = position_jitterdodge(jitter.width = 0.15)) +#+ facet_grid(Dataset~., scales = "fixed") +
  theme_bw() + theme(legend.background = element_blank(),legend.position=c(.7,.95), aspect.ratio = 1,axis.title = element_text(size=14), axis.text = element_text(size=12),
        legend.text = element_text(size=14)) + guides(color=guide_legend(ncol=3)) +
  labs(x='Time (Days)',y=expression('log'[10]*' Count'), color='', fill="") + guides(fill='none') +
	scale_color_manual(labels =c("N", "S",'I'), values = c("Non-synonymous"='#c42929','Synonymous'='#007d6b','Intergenic'='#6a6a6a')) +
	scale_fill_manual(values = c("Non-synonymous"='#c42929','Synonymous'='#007d6b','Intergenic'='#6a6a6a')) +
	geom_smooth(aes(group=Synonymy)) + ylim(0,4)
stp
### stats
compare_means(Count ~ Synonymy , data = st[Synonymy != 'Intergenic'] %>% arrange(Day), group.by = "Day", method = 'wilcox.test', paired=F)

syn_ratio <- st %>% filter(Synonymy != 'Intergenic') %>% pivot_wider(names_from = Synonymy, values_from=Count) %>% mutate(Ratio = `Non-synonymous` / Synonymous)
synStats3 <- as.data.table(compare_means(Count ~ Synonymy, data = st %>% filter(Synonymy != "Intergenic" & Day %in% c(3)) %>% arrange(Mouse), method='wilcox.test', paired = T))
synStats7 <- as.data.table(compare_means(Count ~ Synonymy, data = st %>% filter(Synonymy != "Intergenic" & Day %in% c(7)) %>% arrange(Mouse), method='wilcox.test', paired = T))
synStats14 <- as.data.table(compare_means(Count ~ Synonymy, data = st %>% filter(Synonymy != "Intergenic" & Day %in% c(14)) %>% arrange(Mouse), method='wilcox.test', paired = T))
synStats21 <- as.data.table(compare_means(Count ~ Synonymy, data = st %>% filter(Synonymy != "Intergenic" & Day %in% c(21)) %>% arrange(Mouse), method='wilcox.test', paired = T))
synStats28 <- as.data.table(compare_means(Count ~ Synonymy, data = st %>% filter(!SampleID %in% c(35241) & Synonymy != "Intergenic" & Day %in% c(28)) %>% arrange(Mouse), method='wilcox.test', paired = T))

synStats28

syn_ratio %>% group_by(Day) %>% summarise(Mean = mean(Ratio, na.rm = T), Median = median(Ratio, na.rm = T))
wilcox.test(syn_ratio %>% filter(Day == 3) %>% pull(Ratio), mu=2, alternative = 'greater')
wilcox.test(syn_ratio %>% filter(Day == 7) %>% pull(Ratio), mu=2, alternative = 'greater')
wilcox.test(syn_ratio %>% filter(Day == 14) %>% pull(Ratio), mu=2, alternative = 'greater')
wilcox.test(syn_ratio %>% filter(Day == 21) %>% pull(Ratio), mu=2, alternative = 'greater')
wilcox.test(syn_ratio %>% filter(Day == 28) %>% pull(Ratio), mu=2, alternative = 'greater')

model <- lm(Ratio ~ 1, syn_ratio)
confint(model,level=0.95)


### Final Figure

traj = b + c
traj
fig1 = ta1 + ta2 + ta3 + total + stp + pip + traj + plot_layout(design = "ABCD\nEEFF\nEEFF\nGGGG\nGGGG", nrow = 3) +
	plot_annotation(tag_levels = 'A')
fig1
ggsave('manuscript/Figure3.svg', width=790, height = 1120*0.7, unit = 'px', dpi = 'screen')







