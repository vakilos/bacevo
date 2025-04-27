library(ggplot2)
library(data.table)
library(forcats)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
library(patchwork)
library(ggExtra)
library(ggpointdensity)


df = fread('tmp/invivo_metagenomes_FinalPolyTable_comp28.tsv')
keep_samples <- df %>% filter(Chrom == 'NC_004663' & Coverage_mean >= 100) %>% distinct(SampleID) %>% pull(SampleID)
remove_samples <- df %>% filter(Chrom == 'NC_004663' & Coverage_mean < 100) %>% distinct(SampleID) %>% pull(SampleID)
remove_samples
mask <- fread("tmp/confirmed_invivo.csv")
masked_polys <- unique(mask$PolyID)
df <- df %>% filter(!PolyID %in% masked_polys & TotalCov >= 100 & SampleID %in% keep_samples &
  Ancestry == 'evolved' & Day == 28 & Dataset == 'ploen') %>% distinct(SampleID,PolyID, .keep_all = T)
df <- df %>% mutate(SampleTypeGroup = ifelse(SampleType %in% c("SIP","SID"), "SI", "LI"))
sampletype_color = c('Cecum' = '#f1c50d', 'Feces' = '#59594b', 'Mucus'='#7dc7b0',
                     'SIP' = '#2065ab','SID'='#f57089')
sampletypegroup_color = c('LI'="#d2b53e", "SI"="#2065ab")


pc <- fread('tmp/ploen_binary_pcoa.csv')
pc <- pc %>% mutate(SampleTypeGroup = case_when(SampleType %in% c('SIP','SID') ~ 'SI', .default = 'LI'))
#pc <- left_join(pc, df %>% filter(Chrom == "NC_004663") %>% distinct(SampleID,Coverage_mean), by="SampleID")
com1 <- ggplot(pc %>% filter(SampleID %in% keep_samples), aes(x=PC1, y=PC2, color=SampleTypeGroup, shape=SampleType, label=Day)) +
  geom_point(alpha=0.7, size=5) + scale_color_manual(values=sampletypegroup_color) + labs(x='PC1 (17.56 %)',y='PC2 (3.94 %)') +
  theme(axis.title = element_text(size=14), legend.position = 'bottom') + coord_fixed(ratio = 1)
com1
df$SampleID <- as.character(df$SampleID)



### Choose one of the _coverage_mean.csv , _chromreads.csv , _readcount.csv
sign <- fread("tmp/compartments_stat_sign_polys_readcount.csv")
# sign <- fread("~/bacevo/invivo/tmp/compartments_stat_sign_polys_freq.csv") ### 10 more significant Polys based on Frequency data compared to normalized AltCov
dim(sign %>% filter(padj <= 0.05)) ### 133503
dim(sign %>% filter(qval <= 0.05)) ### 133503


rich <- df %>% group_by(SampleID,SampleType) %>% summarise(PolyCount = n())
rich <- rich %>% mutate(SampleTypeGroup = ifelse(SampleType %in% c("SIP","SID"), "SI", "LI"))
com2 <- ggplot(rich, aes(x=SampleType, y=log10(PolyCount), fill=SampleTypeGroup)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2, alpha=0.7) +
  theme(aspect.ratio = 1.5, axis.title = element_text(size=14), axis.text.x = element_text(size=12, angle = -30), legend.position = 'none') +
  scale_fill_manual(values=sampletypegroup_color) + labs(x='', y= expression("log"[10]*" Count"))
com2

prev2 = fread('tmp/SampleTypeGroupPrev.csv')

prev2 <- merge(prev2, sign %>% select(padj,qval,meanSI,meanLI, PolyID), by='PolyID', all.x = T)
prev2 <- prev2 %>% mutate(GeneLabel = case_when(oldGeneName =='Nan' ~ GeneName, .default = oldGeneName))
setkey(prev2,NULL)
prev2 <- prev2 %>% mutate(ColorGroup = case_when(padj <= 0.05 & is.na(qval) ~ 'Mann-Whitney',padj <= 0.05 & qval > 0.05 ~ 'Mann-Whitney',
                                            is.na(padj) & qval <= 0.05 ~ "Maaslin2", padj > 0.05 & qval <= 0.05 ~ 'Maaslin2',
                                            padj <= 0.05 & qval <= 0.05 ~ 'Consensus', .default='Other'))
prev2 <- prev2 %>% mutate(Appearance = case_when(preSI == 0 & preLI > 0 ~ 'LI',
                                                 preSI > 0 & preLI == 0 ~ 'SI', .default = 'common'))
prev2 <- prev2 %>% mutate(Significant = case_when((ColorGroup == 'Maaslin2') | (ColorGroup == 'Consensus') ~'Significant association', .default = 'Other'))
prev2 <- prev2 %>% mutate(Diff = meanSI - meanLI)
prev2 <- prev2 %>% mutate(Diffreq = freqSI - freqLI)

#### Normalized read support for alternative allele correlates with frequency in the population.
rn <- prev2 %>% filter(Significant == 'Significant association' & meanSI > 0 & meanLI > 0) %>% pivot_longer(cols = c('meanSI', 'meanLI'), names_to = 'ReadsNorm',values_to ='ReadsNormValue') %>% select(PolyID,ReadsNorm,ReadsNormValue)
fr <- prev2 %>% filter(Significant == 'Significant association' & meanSI > 0 & meanLI > 0) %>% pivot_longer(cols = c('freqSI', 'freqLI'), names_to = 'Freq',values_to ='FreqValue') %>% select(PolyID,Freq, FreqValue)
ou <- left_join(rn, fr, by='PolyID', relationship = 'many-to-many')
correlation <- cor(prev2$Diff, prev2$Diffreq, method = "pearson")

ggplot(prev2, aes(x=log10(Diff), y=Diffreq)) + geom_point(alpha=0.3) + geom_smooth(method='lm') +
annotate("text", x = -5, y = 1.1,     # Add the correlation coefficient text annotation
           label = paste("Pearson r =", round(correlation, 2)),
           color = "red", size = 5)

p1 <- ggplot(prev2, aes(x=meanLI, y=meanSI, color=ColorGroup)) + geom_point(alpha=.4) +
  scale_color_manual(values=c('Maaslin2'='red', 'Consensus'='green', "Mann-Whitney"='blue', "Other"='darkgrey'))
p1
p2 <- ggplot(prev2, aes(x=freqLI, y=freqSI, color=ColorGroup)) + geom_point(alpha=.4) +
  scale_color_manual(values=c('Maaslin2'='red', 'Consensus'='green', "Mann-Whitney"='blue', "Other"='darkgrey'))
p2
p3 <- ggplot(prev2, aes(x=preLI, y=preSI, color=ColorGroup)) + geom_point(alpha=.4) +
  scale_color_manual(values=c('Maaslin2'='red', 'Consensus'='green', "Mann-Whitney"='blue', "Other"='darkgrey'))
p3
p1 + p2 + p3 + plot_layout(guides='collect')
significant_hits <- of %>% filter(ColorGroup == 'Consensus') %>% pull(PolyID)


### make boxplots with poly count of unique to the LI, SI and shared polys.
prev2 %>% group_by(Appearance, Significant) %>% summarize(Count = n())
com3 <- ggplot(prev2 %>% group_by(Appearance, Significant) %>% summarize(Count = n()), aes(x=Appearance, y=Count, fill=Significant)) +
  geom_bar(width=0.5, stat = 'identity') +
  scale_fill_manual(name='',values=c('Significant association'='#d76650', 'Other'='darkgrey')) + scale_y_continuous(labels=scales::label_number(scale=1e-3, suffix="K")) +
  theme(legend.position = 'bottom',aspect.ratio = 1.5, axis.text.y = element_text(size=14),  axis.text.x = element_text(size=10),
        axis.title = element_text(size=14)) + labs(x='', y='Polymorphisms')
com3

### The portion of mutations that is LI unique is not in the statistically significant differentially abundant set of Maaslin2,
### because they does not satisfy the 0.2 minimum prevalence threshold.
#ggplot(prev2 %>% filter(Appearance == 'LI unique'), aes(x=preLI)) + geom_density()


### select present - absent polys between sampletypgroups - there are no cases of present in LI and absent from SI which might be also indicative
## of some directionality - i.e genetic diversity emerged in upper parts of the gut can be selected against in lower parts and not the other way round
selection <- prev2 %>% filter(PolyID %in% significant_hits & (preSI == 0 & preLI > 0) | (preSI > 0 & preLI == 0) & preSI > 0.5 & freqSI >= 0.1 & GeneName != "Intergenic") %>%
  arrange(Diff) %>% mutate(PolyID=factor(PolyID, PolyID))

prev2 <- prev2 %>% mutate(GeneLabel = case_when(oldGeneName =='Nan' ~ GeneName, .default = oldGeneName))
prev2 <- prev2 %>% mutate(Label = paste(paste(PolyID,oldGeneName, sep=' ('),')',sep=''))
prev2 <- prev2 %>% mutate(Log2 = log2( (meanSI + (10^-10)) / meanLI + (10^-10)))
#ggplot(prev2, aes(x=Log2)) + geom_density()


prev2 <- prev2 %>% mutate(Genegroup = case_when(GeneName == 'Intergenic' ~ 'no', .default = 'yes'))

library(ggforce)
prev.sum <- prev2 %>% filter(Appearance == 'SI') %>% group_by(Significant) %>% summarize(SD = sd(meanSI), MEAN = mean(meanSI))

genes_of_significant <- prev2 %>% filter(Significant == 'Significant association') %>% distinct(GeneName) %>% pull(GeneName)
length(genes_of_significant)



ggplot(prev2 %>% filter(str_detect(PolyID,'^NC_004663')), aes(x=meanSI)) + geom_density() + xlim(0,1e-07)

prev2 %>% filter( ColorGroup == 'Consensus' & Diff >= 2.5e-06 & GeneName != 'Intergenic') %>% distinct(oldGeneName)
a <- ggplot(prev2 %>% filter( ColorGroup == 'Consensus' & Diff >= 2.5e-06 & GeneName != 'Intergenic'))  +
  geom_segment( aes(x=fct_reorder(Label,Diff), xend=fct_reorder(Label,Diff), y=meanLI, yend=meanSI), color="grey") +
  geom_point( aes(x=fct_reorder(Label,Diff), y=meanSI), color="#2065ab", alpha=0.8,size=3 ) +
  geom_point( aes(x=Label, y=meanLI), color="#d2b53e", alpha =.8 , size=3) +
  coord_flip() + theme(legend.position = "bottom") + xlab("") + ylab("Normalized read depth")+
  theme(axis.text.y = element_text(size=8), aspect.ratio = 2, axis.title = element_text(size=14))
a
a2 <- ggplot(prev2 %>% filter( ColorGroup == 'Consensus' & Diff >= 2.5e-06 & GeneName != 'Intergenic'))  +
  geom_bar( aes(x=fct_reorder(Label,Diff), y=preSI),stat='identity', fill='#2065ab') +
  coord_flip() + theme(legend.position = "none") + xlab("") + ylab("Prevalence")+
  scale_y_continuous(labels = scales::label_percent(), limits=c(0,1))  +
  theme(axis.title = element_text(size=14), axis.text.x = element_text(size=8),aspect.ratio = 4, axis.text.y = element_blank(), axis.ticks.y = element_blank())
a2
a + a2


ap <- ggplot(prev2 %>% filter(str_detect(PolyID,'^NC_004703') & ColorGroup == 'Consensus' & Diff >= 2.1e-06 & GeneName != 'Intergenic'))  +
  geom_segment( aes(x=fct_reorder(Label,Diff), xend=fct_reorder(Label,Diff), y=meanLI, yend=meanSI), color="grey") +
  geom_point( aes(x=fct_reorder(Label,Diff), y=meanSI), color="#2065ab", alpha=0.8,size=3 ) +
  geom_point( aes(x=Label, y=meanLI), color="#d2b53e", alpha =.8 , size=3) +
  coord_flip() + theme(legend.position = "bottom") + xlab("") + ylab("Normalized read depth")+
  theme(axis.text.y = element_text(size=8), aspect.ratio = 2, axis.title = element_text(size=14))

ap2 <- ggplot(prev2 %>% filter(str_detect(PolyID,'^NC_004703') & ColorGroup == 'Consensus' & Diff >= 2.1e-06 & GeneName != 'Intergenic'))  +
  geom_bar( aes(x=fct_reorder(Label,Diff), y=preSI),stat='identity', fill='#2065ab') +
  coord_flip() + theme(legend.position = "none") + xlab("") + ylab("Prevalence")+
  scale_y_continuous(labels = scales::label_percent(), limits=c(0,1))  +
  theme(axis.title = element_text(size=14), axis.text.x = element_text(size=8),aspect.ratio = 4, axis.text.y = element_blank(), axis.ticks.y = element_blank())

a + a2 + ap + ap2

prev2 %>% filter(ColorGroup =='Consensus' & meanSI > 0 & meanLI > 0 & GeneName != 'Intergenic') %>% distinct(oldGeneName)
b <- ggplot(prev2 %>% filter(ColorGroup =='Consensus' & meanSI > 0 & meanLI > 0 & GeneName != 'Intergenic'))  +
  geom_segment( aes(x=fct_reorder(Label,Diff), xend=fct_reorder(Label,Diff), y=meanLI, yend=meanSI), color="grey") +
  geom_point( aes(x=fct_reorder(Label,Diff), y=meanSI), color="#2065ab", alpha=0.8, size=4) +
  geom_point( aes(x=fct_reorder(Label,Diff), y=meanLI), color="#d2b53e", alpha=0.8, size=4) +
  coord_flip()  + xlab("") + ylab("Normalized read depth")+
   theme(axis.title = element_text(size=14),legend.position = 'bottom',axis.text.y = element_text(size=10), aspect.ratio = 2)
b

d <- ggplot(prev2 %>% filter(ColorGroup == 'Consensus' & meanSI > 0 & meanLI > 0 & GeneName != 'Intergenic') %>%
pivot_longer(cols=starts_with('pre'))) +
  geom_bar( aes(x=fct_reorder(Label,Diff), y=value,fill=name), stat='identity', position = 'dodge', width=0.5) +
  coord_flip()  + xlab("") + ylab("Prevalence")+
  scale_fill_manual(values=c("preSI"="#2065ab", 'preLI'="#d2b53e")) +
  scale_y_continuous(labels = scales::label_percent(), limits=c(0,1)) +
  theme(axis.title = element_text(size=14),axis.text.x = element_text(size=8),legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), aspect.ratio = 4)

d
b + d


(com1 + com2 + com3 + plot_layout(guides='auto', heights = c(20,1,1), widths=c(8,4,4))) / (a + a2 + b + d + plot_layout(heights = c(200,1,0),widths = c(1,0.5,1,0.5)))
ggsave("~Figure5_compartments.svg", height=10, width=14)

(com1 + com2 + com3 + plot_layout(guides='auto', heights = c(20,1,1), widths=c(8,4,4)))
ggsave('Figure5_compartments_p1.svg',height = 10,width=14)
ggsave('~Figure5_compartments_p1.png')
(a + a2 + b + d + plot_layout(heights = c(200,100,0),widths = c(1,0.5,1,0.5)))
ggsave('Figure5_compartments_p2.svg',height = 10,width=14)
ggsave('Figure5_compartments_p2.png')

