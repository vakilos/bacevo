library(data.table)
library(ggplot2)
library(ggrepel)
library(forcats)
library(ggpubr)
library(dplyr)
library(patchwork)
library(glue)
library(ggvenn)
setwd("~/bacevo/invitro")


df = fread("invitro_masked_FinalPolyTable.tsv", sep='\t')

anc = fread("AncestralMutations.csv")
ggplot(anc, aes(x=Freq)) + geom_histogram()

ggplot(df %>% filter(TotalCov >= 100), aes(x=as.factor(Replicate),y=Freq, color=Ancestry)) +
  geom_boxplot() + facet_wrap(~Mix, scales='free')

df %>% filter(Freq > 0.5 & Ancestry == 'evolved') %>% arrange(Chrom,Pos) %>% distinct(Chrom,Pos, GeneName,oldGeneName)



uq = fread("tmp/Unique_SNPs_per_library.csv")
uq = fread("tmp/InvitroFinalPolyTable_unique_per_library.tsv", sep='\t')
mu = fread('tmp/MixUnique.tsv', sep='\t')
mu <- mu %>% mutate(MixUnique = case_when(MixUnique == "0" ~ "2 Media",
                                          MixUnique == "1" ~ "Medium 1 only",
                                          MixUnique == "2" ~ "Medium 2 only",
                                          MixUnique == "3" ~ "Medium 3 only",
                                          MixUnique == "all" ~ "All Media"))
mu <- mu %>% mutate(Mix = case_when(Mix == 1 ~ "Medium 1",
                                    Mix == 2 ~ "Medium 2",
                                    Mix == 3 ~ "Medium 3"))
mu

f1a <- ggplot(mu[SampleProfile == 'metagenome'], aes(x=MixPrevalence,y=MixMeanFreq, color=MixUnique)) +
  geom_jitter(size=1, alpha=0.4) + facet_wrap(~Mix) +
  scale_color_manual(values=c("All Media"="#67b27e",
                              "Medium 1 only"="#e99427",
                              "Medium 2 only"="#37aced",
                              "Medium 3 only"="#eba1cf",
                              "2 Media" = '#b6d595')) +
  theme(axis.title = element_text(size=24), axis.text = element_text(size=16), strip.text = element_text(size=18),
        legend.title = element_text(size=16), legend.text = element_text(size=14), aspect.ratio = 1) +
  labs(x='Medium Prevalence', y='Mean Frequency in Medium', color='Medium Specificity')

f1a
ggsave("Figures/MixMeanFreq_MixPrevalence.svg")
ggsave("Figures/MixMeanFreq_MixPrevalence.pdf")
f1b <- ggplot(mu[SampleProfile == 'metagenome'], aes(x=MixPrevalence, fill=MixUnique)) +
  geom_bar( color='black', alpha=0.8)  + facet_wrap(~Mix) +
  scale_fill_manual(values=c("All Media"="#67b27e",
                              "Medium 1 only"="#e99427",
                              "Medium 2 only"="#37aced",
                              "Medium 3 only"="#eba1cf",
                              "2 Media" = '#b6d595')) +
  theme(axis.title = element_text(size=24), axis.text = element_text(size=16), strip.text = element_text(size=18),
        legend.title = element_text(size=16), legend.text = element_text(size=14), aspect.ratio = 1) +
  labs(x='Medium Prevalence', y='Polymorphisms', fill='Medium Specificity')
f1b
ggsave("Figures/MixPrevalenceCount.svg")
ggsave("Figures/MixPrevalenceCount.pdf")
### Venn diagram
ven = list(`Medium 1`=unique(mu[Mix == "Medium 1", PolyID ]),
           `Medium 2`= unique(mu[Mix == "Medium 2", PolyID]),
           `Medium 3`= unique(mu[Mix == "Medium 3", PolyID]))

f1c <- ggvenn(ven, fill_color = c("#e99427", "#37aced", "#eba1cf"),
stroke_size = 0.1,, fill_alpha = 0.6, set_name_size = 12, text_size= 10 )
f1c
ggsave("Figures/VennDiagramm_polys_per_Mix.svg")
ggsave("Figures/VennDiagramm_polys_per_Mix.pdf")

