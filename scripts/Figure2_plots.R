library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gggenes)
library(ggfittext)
library(stringr)
library(forcats)
library(patchwork)
library(ggrepel)
library(ggpattern)

#setwd('gene_order')
cl = fread('tmp/my_clusters_stats.csv')
tr = fread('tmp/confirmed_invivo.csv')
target_genes <- tr %>% distinct(GeneName, .keep_all = T) %>% pull(GeneName)
cl <- cl[Cluster != 'no cluster' & Signature == TRUE & GeneName %in% target_genes] %>% distinct(SampleID,PolyID, .keep_all = T)
cl %>% filter(oldGeneName == 'BT0436' & Cluster != 'no cluster') %>%
  ggplot(aes(x=as.factor(Day),y=Freq,color=oldGeneName, group=PolyID)) +
  geom_line(alpha=1) + facet_wrap(~Mouse)

### From the reported shufflons BT4520 signature polys appear in very low frequency only in one specimen..similartly BT0436 appears in one specimen.
colorBT2268 = c('BT2259'='#3aaca3','BT2260'='#f1c50d','BT2261'='grey','BT2262'='grey','BT2263'='#3aaca3','BT2264'='#f1c50d',
              'BT2265'='grey','BT2266'='grey','BT2267'='#cf6b60','BT2268'='#f1c50d','BT2269'='#3aaca3',
                'BT2270'='grey', na.value='grey')
colorBT1042 = c('BT1040' = '#f1c50d', 'BT1041'='#cf6b60', "BT1042"='#f1c50d',"BT1043"='#3aaca3',
                "BT1044"='grey','BT1045'='grey','BT1046'='#f1c50d', na.value='grey')
colorBT3239 = c("BT3235" = 'grey', "BT3236" = 'grey', 'BT3237'='grey', 'BT3238'='#3aaca3','BT3239'='#f1c50d','BT3240'='#f1c50d',
                "BT3241"='#3aaca3', "BT3242"='grey','BT3243'='grey','BT3244'='grey',
                'BT3245'='grey','BT3246'='grey')

colorBT0436 = c("BT0433"='grey', "BT0434"='grey','BT0435'='grey','BT0436'='#67b27e', 'BT0437'='#b2258f','BT0438'='grey',
                'BT0439'='#f1c50d', 'BT0440'='#3aaca3','BT0441'='grey','BT0442'='grey',"BT0443"='grey','BT0444'='grey',
                'BT0445'='grey', 'BT0446'='grey','BT0447'='grey','BT0448'='grey','BT0449'='grey','BT0450'='grey',
                'BT0451'='#3aaca3','BT0452'='#f1c50d','BT0453'='#b2258f','BT0454'='grey','BT0455'='grey', na.value='grey')

colorBT3477 = c("BT3477"='grey','BT3478'='#cf6b60','BT3479'='#a24057','BT3480'='#a24057','BT3481'='grey','BT3482'='grey','BT3483'='#f1c50d',
                'BT3484'='#3aaca3','BT3485'='grey', "BT3486"='grey', 'BT3487'='grey','BT3488'='grey','BT3489'='grey','BT3490'='grey',
                'BT3491'='grey', "BT3492"='grey',"BT3494"='#f1c50d', "BT3495"='#3aaca3',"BT3496"='grey',"BT3497"='grey', "BT3498"='grey',
                "BT3499"='grey', "BT3500"='grey', "BT3501"='purple',na.value='grey')

colorBT4540 = c("BT4540"='#67b27e',"BT4541"='#67b27e',"BT4542"='#67b27e', "BT4543"='#67b27e',"BT4544"='#cf6b60', na.value='grey')

colorBT4520 = c('BT4519'='grey', 'BT4520'='#67b27e', 'BT4521'='#cf6b60', 'BT4522'='#67b27e', "BT4523"='#67b27e', na.value='grey')

colorDiftype = c("INS" = 'green',"DEL" ='red', na.value='black')



colorShufflon = c("BT0436"='pink',
                  "BT1040"='green',"BT1042"='green',
                  "BT2260"='blue',"BT2264"='blue',"BT2268"='blue',
                  "BT3239"='purple','BT3240'='purple',
                  "BT4520"='brown',
                  'BT3477'='orange',
                  "BT4540"='yellow',"BT4541"='yellow',"BT4543"='yellow',"Nan"='grey')
ggplot(cl %>% filter(oldGeneName != "Nan" & oldGeneName %in% c("BT4522")) , aes(x=as.factor(Day),y=Freq,color=oldGeneName, group=PolyID))  +
geom_line(alpha=0.5) + facet_wrap(~Mouse)+
  theme(legend.position = 'bottom') + scale_color_manual(values=colorShufflon) +
  labs(title='Poly trajectories associated with a shufflon signature',
       color='In color genes mutated only invivo', x="Day post inoculation", y='Frequency in population')

shufflons =  c('BT4520','BT1042',"BT2268",'BT0436',
           "BT3477",'BT4540','BT3239')
shufflon = shufflons[4]



vgm = fread(paste0(shufflon,"_variants_frequency_per_mouse.csv"))
cb= fread(paste0(shufflon,'_relative_coords.csv'))
cb <- cb %>% mutate(NodeCount = max(nodeN), .by=SampleID)
cb <- cb %>% mutate(SampleID = str_replace(SampleID,'S35254','Ancestor'))

### do not show variants with freq < 0.05 in the population
cb <- cb %>% filter(VariantFreq >= 0.01 & str_detect(VariantName,'V') | str_detect(VariantName,'V',negate = T))

Ngenes <- cb %>% filter(SampleID == 'RefSeq') %>% arrange(generelstart) %>% distinct(oldName) ## genes within shufflon


##### Plots the frequency of each detected variant per mouse
ggplot(vgm[Mouse != 0], aes(x=fct_reorder(VariantName,Freq, .fun = mean), y=Freq)) +
  geom_bar(stat='identity') + facet_wrap(~Mouse, scales = 'free_y') + coord_flip() +
  labs(title=paste0(shufflon," gene cluster"))

### Variant frequency in the population (isolates) - only for one node (complete) variants -> will not sum up to 1
p1 <- cb %>% filter(SampleID != "RefSeq" & NodeCount == 1) %>% distinct(string, .keep_all = T) %>%
  filter(VariantFreq >= 0.0 & str_detect(VariantName,'V') | str_detect(VariantName,'V',negate = T)) %>%
  ggplot(aes(x=fct_reorder(VariantName,VariantFreq),y=VariantFreq, fill=VariantName, label=string)) +theme_classic()+
  geom_bar(stat='identity') + labs(y='Variant frequency', x='Variant',
                                   title=paste0(shufflon," gene cluster (",first(Ngenes)," - ",last(Ngenes),")")) +
  theme(legend.position = 'none', aspect.ratio = 1)  +
  coord_flip()
p2 <- cb %>% filter(SampleID != "RefSeq" & NodeCount == 1) %>% distinct(VariantName, .keep_all = T) %>%
  filter(VariantFreq >= 0.0 & str_detect(VariantName,'V') | str_detect(VariantName,'V',negate = T)) %>%
  ggplot(aes(x=fct_reorder(VariantName,VariantFreq),y=NodeCount, fill=VariantName, group=VariantName)) +
  geom_bar(stat='identity') + theme_classic() +
  theme(legend.position = 'none', aspect.ratio = 1, axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(x='') +
  coord_flip() + geom_hline(yintercept = 1, linetype='dashed', color='black')
p1 + p2 + plot_layout(widths = c(1,1), guides = 'collect')




#####
p1_v1 <- cb %>% distinct(string, .keep_all = T) %>%
  filter(NodeCount == 1 & str_detect(VariantName,'V') | str_detect(VariantName,'V',negate = T)) %>%
  ggplot(aes(x=fct_reorder(VariantName,VariantFreq),y=VariantFreq, label=string)) +theme_classic()+
  geom_bar(stat='identity',fill='grey', width=0.7) + labs(y='Variant frequency', x='') +
  theme(legend.position = 'none', aspect.ratio = 3, axis.ticks.y = element_blank(),axis.text.y = element_blank())  +
  coord_flip()

p1_v2 <- cb %>% distinct(string, .keep_all = T) %>%
  filter(NodeCount == 1 ) %>%
  ggplot(aes(x=fct_reorder(VariantName,VariantFreq),y=VariantFreq, label=string)) +theme_classic()+
  geom_bar(stat='identity',fill='grey', width=0.7) + labs(y='Variant frequency', x='') +
  theme(legend.position = 'none', aspect.ratio = 3, axis.ticks.y = element_blank(),axis.text.y = element_blank())  +
  coord_flip()
p1_v1
#cb %>% distinct(string, .keep_all = T) %>% select(VariantName, string)


p3 <- cb %>% mutate(genestrand = ifelse(genestrand == "-", 0,1)) %>%
  ggplot(aes(xmin=generelstart, xmax=generelend, y=fct_reorder(VariantName,VariantFreq), fill=as.factor(nodeN),label=RefBestHit,
                   forward = genestrand)) +
  geom_gene_arrow(arrowhead_width = unit(0.1,'cm'), arrow_body_height = unit(0.5,'cm')) + geom_fit_text(size=8) +
  scale_x_continuous(breaks=seq(0,max(cb$generelend), 500))  +
  theme(axis.text.x = element_text(angle = 50,hjust=1),
        legend.position = 'bottom') +
  labs(title = paste0(shufflon," gene cluster (",first(Ngenes)," - ",last(Ngenes),")"), y='Variant', fill='Sample contig / node')

p4 <- cb %>% mutate(genestrand = ifelse(genestrand == "-", 0,1)) %>%
  ggplot(aes(xmin=generelstart, xmax=generelend, y=fct_reorder(VariantName,VariantFreq), fill=as.factor(nodeN),label=oldName,
                   forward = genestrand)) +
  geom_gene_arrow(arrowhead_width = unit(0.1,'cm'), arrow_body_height = unit(0.5,'cm')) + geom_fit_text(size=8) +
  scale_x_continuous(breaks=seq(0,max(cb$generelend), 500)) +
  theme(axis.text.x = element_text(angle = 50,hjust=1), legend.position = 'bottom') +
  labs(title = paste0(shufflon," gene cluster (",first(Ngenes)," - ",last(Ngenes),")"),y = "Variant",fill='Sample contig / node')

layout = "AAAAAB"
p3 + p1_v1 + plot_layout(guides = 'keep', design = layout)
p4 + p1_v1 + plot_layout(guides = 'keep', design = layout)



p6 <- cb %>% filter(NodeCount == 1) %>% mutate(genestrand = ifelse(genestrand == "-", 0,1)) %>%
  ggplot(aes(xmin=generelstart, xmax=generelend, y=fct_reorder(VariantName, VariantFreq), fill=info,label=oldName,
                   forward = genestrand)) +
  geom_gene_arrow(arrowhead_width = unit(0.1,'cm'), arrow_body_height = unit(0.5,'cm')) + geom_fit_text(size=8) +
  scale_x_continuous(breaks=seq(0,max(cb$generelend), 500)) +
  theme(legend.position = 'bottom',axis.text.x = element_text(angle = 50,hjust=1)) +  guides(fill = guide_legend(nrow = 2)) +
  labs(title = paste0(shufflon," gene cluster"))
p6 + p1_v2 + plot_layout(guides = 'keep', design = layout)


### BT1042
shufflon = shufflons[2]
fra = fread(paste0(shufflon,'_fragments_motifs_rel_coords.csv'))
fra <- fra %>% mutate(NodeCount = max(nodeN), .by=SampleID)
fra <- fra %>% mutate(genestrand = ifelse(genestrand == "-", 0,1))

Ngenes <- fra %>% filter(SampleID == 'RefSeq' & oldName != '') %>% arrange(generelstart) %>% distinct(oldName)
oneNodeSamples = fra %>% group_by(SampleID) %>% distinct(nodeID) %>% mutate(Nodes = n()) %>%
             filter(Nodes == 1) %>% distinct(SampleID) %>% select(SampleID)
oneNodeSamples = unlist(oneNodeSamples)

library(extrafont)
library(svglite)
ggplot(fra %>% filter(SampleID %in% oneNodeSamples) %>% distinct(VariantName, RefBestHit, .keep_all = T),aes(xmin=generelstart, xmax=generelend, y=VariantName,label=oldName,
                   forward = genestrand,fill=oldName)) +
  scale_fill_manual(values=colorBT1042, breaks=c("BT1040",'BT1043','BT1041'),labels=c('susC','susD','site-specific integrase')) +
  geom_gene_arrow(arrowhead_width = unit(0.1,'cm'), arrow_body_height = unit(0.5,'cm')) +
  #geom_rect_pattern(data= fra %>% filter(str_detect(motif, 'BT1042_rs')), aes(xmin=motifrelstart,xmax=motifrelstart+1, ymin=VariantName, ymax=VariantName, pattern_fill=motif)) +
  geom_gene_arrow(data=fra[fragName != ''],arrowhead_width = unit(0.0,'cm'),arrow_body_height = unit(0.5,'cm'),
                  aes(xmin=fragrelstart, xmax=fragrelend, fill=fragorigin, color=diftype),alpha=0.7, size=1) +
  scale_color_manual(values=colorDiftype, breaks=c("INS",'DEL'), labels = c('insertion','deletion')) +
  geom_fit_text(size=10, fontface = 'bold.italic') +
  #geom_feature(data= fra %>% filter(str_detect(motif, 'BT1042_rs')) ,aes(x = motifrelstart, y =VariantName , forward = motifori),feature_height = unit(4,'mm'), feature_width = unit(4,'mm')) +
  #geom_text_repel(data=fra %>% distinct(SampleID, fragName, .keep_all = T),aes(x=fragrelstart, y=VariantName, label=fragName), angle=00,size=2.5, nudge_y =-0.2,segment.colour = NA) +theme_genes() +
  #geom_text(data = fra %>% filter(str_detect(motif, 'BT1042_rs')), aes(x=motifrelstart,y=VariantName, label=motif), size=3, nudge_y=0.3, angle=20) +
  scale_x_continuous(breaks=seq(0,max(fra[!is.na(generelend)]$generelend), 500)) +
  theme_genes() + theme(axis.text.x = element_text(angle = 50,hjust=1, size=8),
        legend.position = 'right', aspect.ratio = 0.5) +
  labs(title = paste0(shufflon," gene cluster (",first(Ngenes)," - ",last(Ngenes),")"), y='Variant', fill='', x= "Genomic distance (bp)")
ggsave(paste0("shufflons/",shufflon,".svg"))

## BT2268
shufflon = shufflons[3]
fra = fread(paste0(shufflon,'_fragments_motifs_rel_coords.csv'))
fra <- fra %>% mutate(NodeCount = max(nodeN), .by=SampleID)
fra <- fra %>% mutate(genestrand = ifelse(genestrand == "-", 0,1))
Ngenes <- fra %>% filter(SampleID == 'RefSeq' & oldName != '') %>% arrange(generelstart) %>% distinct(oldName)
oneNodeSamples = fra %>% group_by(SampleID) %>% distinct(nodeID) %>% mutate(Nodes = n()) %>%
             filter(Nodes == 1) %>% distinct(SampleID) %>% select(SampleID)
oneNodeSamples = unlist(oneNodeSamples)
ggplot(fra %>% filter(SampleID %in% oneNodeSamples) %>% distinct(VariantName, RefBestHit, generelstart, .keep_all = T),aes(xmin=generelstart, xmax=generelend, y=VariantName,label=oldName,
                   forward = genestrand,fill=oldName)) +
  scale_fill_manual(values=colorBT2268, breaks=c("BT2260",'BT2263','BT2267'),labels=c('susC','susD','site-specific integrase')) +
  geom_gene_arrow(arrowhead_width = unit(0.1,'cm'), arrow_body_height = unit(0.5,'cm')) +
  geom_gene_arrow(data=fra[fragName != ''],arrowhead_width = unit(0.0,'cm'),arrow_body_height = unit(0.5,'cm'),
                  aes(xmin=fragrelstart, xmax=fragrelend, fill=fragorigin, color=diftype),alpha=0.7, size=1) +
  scale_color_manual(values=colorDiftype, breaks=c("INS",'DEL'), labels = c('insertion','deletion')) +
  geom_fit_text(size=10, fontface = 'bold.italic') +
  geom_feature(data= fra %>% filter(str_detect(motif, 'olga')) ,aes(x = motifrelstart, y =VariantName , forward = motifori),feature_height = unit(4,'mm'), feature_width = unit(4,'mm')) +
  geom_text_repel(data=fra %>% distinct(SampleID, fragName, .keep_all = T),aes(x=fragrelstart, y=VariantName, label=fragName), angle=00,size=2.5, nudge_y =-0.2,segment.colour = NA) +theme_genes() +
  geom_text_repel(data = fra %>% filter(str_detect(motif, 'olga')), aes(x=motifrelstart,y=VariantName, label=motif), size=3, nudge_y=0.3, angle=20, segment.colour =NA) +
  scale_x_continuous(breaks=seq(0,max(fra[!is.na(generelend)]$generelend), 500)) +
  theme_genes() + theme(axis.text.x = element_text(angle = 50,hjust=1, size=8),
        legend.position = 'right', aspect.ratio = 0.5) +
  labs(title = paste0(shufflon," gene cluster (",first(Ngenes)," - ",last(Ngenes),")"), y='Variant', fill='', x= "Relative position")
ggsave(paste0("/shufflons/",shufflon,".svg"), height=7.161, width=13.357)

## BT3239
shufflon = shufflons[7]
fra = fread(paste0(shufflon,'_fragments_motifs_rel_coords.csv'))
fra <- fra %>% mutate(NodeCount = max(nodeN), .by=SampleID)
fra <- fra %>% mutate(genestrand = ifelse(genestrand == "-", 0,1))
Ngenes <- fra %>% filter(SampleID == 'RefSeq' & oldName != '') %>% arrange(generelstart) %>% distinct(oldName)
oneNodeSamples = fra %>% group_by(SampleID) %>% distinct(nodeID) %>% mutate(Nodes = n()) %>%
             filter(Nodes == 1) %>% distinct(SampleID) %>% select(SampleID)
oneNodeSamples = unlist(oneNodeSamples)
ggplot(fra %>% filter((SampleID %in% oneNodeSamples & VariantFreq >= 0.01) | (VariantName == "RefSeq variant")) %>% distinct(VariantName, RefBestHit, generelstart, .keep_all = T),aes(xmin=generelstart, xmax=generelend, y=VariantName,label=oldName,
                   forward = genestrand,fill=oldName)) +
  scale_fill_manual(values=colorBT3239, breaks=c("BT3239",'BT3241'),labels=c('susC','susD')) +
  geom_gene_arrow(arrowhead_width = unit(0.1,'cm'), arrow_body_height = unit(0.5,'cm')) +
  #geom_gene_arrow(data=fra[fragName != ''],arrowhead_width = unit(0.0,'cm'),arrow_body_height = unit(0.5,'cm'),
  #                aes(xmin=fragrelstart, xmax=fragrelend, fill=fragorigin, color=diftype),alpha=0.7, size=1) +
  scale_color_manual(values=colorDiftype, breaks=c("INS",'DEL'), labels = c('insertion','deletion')) +
  geom_fit_text(size=10, fontface = 'bold.italic') +
  #geom_feature(data= fra %>% filter(str_detect(motif, 'olga')) ,aes(x = motifrelstart, y =VariantName , forward = motifori),feature_height = unit(4,'mm'), feature_width = unit(4,'mm')) +
  #geom_text_repel(data=fra %>% distinct(SampleID, fragName, .keep_all = T),aes(x=fragrelstart, y=VariantName, label=fragName), angle=00,size=2.5, nudge_y =-0.2,segment.colour = NA) +theme_genes() +
  #geom_text_repel(data = fra %>% filter(str_detect(motif, 'olga')), aes(x=motifrelstart,y=VariantName, label=motif), size=3, nudge_y=0.3, angle=20, segment.colour =NA) +
  scale_x_continuous(breaks=seq(0,max(fra[!is.na(generelend)]$generelend), 500)) +
  theme_genes() + theme(axis.text.x = element_text(angle = 50,hjust=1, size=8),
        legend.position = 'right', aspect.ratio = 0.5) +
  labs(title = paste0(shufflon," gene cluster (",first(Ngenes)," - ",last(Ngenes),")"), y='Variant', fill='', x= "Relative position")
ggsave(paste0("shufflons/",shufflon,".svg"), height=7.161, width=13.357)


#### BT4540
shufflon = shufflons[6]
fra = fread(paste0(shufflon,'_fragments_motifs_rel_coords.csv'))
fra <- fra %>% mutate(NodeCount = max(nodeN), .by=SampleID)
fra <- fra %>% mutate(genestrand = ifelse(genestrand == "-", 0,1))
Ngenes <- fra %>% filter(SampleID == 'RefSeq' & oldName != '') %>% arrange(generelstart) %>% distinct(oldName)
oneNodeSamples = fra %>% group_by(SampleID) %>% distinct(nodeID) %>% mutate(Nodes = n()) %>%
             filter(Nodes == 1) %>% distinct(SampleID) %>% select(SampleID)
oneNodeSamples = unlist(oneNodeSamples)
ggplot(fra %>% filter((SampleID %in% oneNodeSamples & VariantFreq >= 0.01) | (VariantName == "RefSeq variant")) %>% distinct(VariantName, RefBestHit, generelstart, .keep_all = T),aes(xmin=generelstart, xmax=generelend, y=VariantName,label=oldName,
                   forward = genestrand,fill=oldName)) +
  scale_fill_manual(values=colorBT4540, breaks=c("BT4540",'BT4544'),labels=c('restriction endonuclease','site-specific integrase')) +
  geom_gene_arrow(arrowhead_width = unit(0.1,'cm'), arrow_body_height = unit(0.5,'cm')) +
  #geom_gene_arrow(data=fra[fragName != ''],arrowhead_width = unit(0.0,'cm'),arrow_body_height = unit(0.5,'cm'),
  #                aes(xmin=fragrelstart, xmax=fragrelend, fill=fragorigin, color=diftype),alpha=0.7, size=1) +
  scale_color_manual(values=colorDiftype, breaks=c("INS",'DEL'), labels = c('insertion','deletion')) +
  geom_fit_text(size=10, fontface = 'bold.italic') +
  #geom_feature(data= fra %>% filter(str_detect(motif, 'olga')) ,aes(x = motifrelstart, y =VariantName , forward = motifori),feature_height = unit(4,'mm'), feature_width = unit(4,'mm')) +
  #geom_text_repel(data=fra %>% distinct(SampleID, fragName, .keep_all = T),aes(x=fragrelstart, y=VariantName, label=fragName), angle=00,size=2.5, nudge_y =-0.2,segment.colour = NA) +theme_genes() +
  #geom_text_repel(data = fra %>% filter(str_detect(motif, 'olga')), aes(x=motifrelstart,y=VariantName, label=motif), size=3, nudge_y=0.3, angle=20, segment.colour =NA) +
  scale_x_continuous(breaks=seq(0,max(fra[!is.na(generelend)]$generelend), 500)) +
  theme_genes() + theme(axis.text.x = element_text(angle = 50,hjust=1, size=8),
        legend.position = 'right', aspect.ratio = 0.5) +
  labs(title = paste0(shufflon," gene cluster (",first(Ngenes)," - ",last(Ngenes),")"), y='Variant', fill='', x= "Relative position")
ggsave(paste0("shufflons/",shufflon,".svg"), height=7.161, width=13.357)

##### BT0436
shufflon = shufflons[4]
fra = fread(paste0(shufflon,'_fragments_motifs_rel_coords.csv'))
fra <- fra %>% mutate(NodeCount = max(nodeN), .by=SampleID)
fra <- fra %>% mutate(genestrand = ifelse(genestrand == "-", 0,1))
Ngenes <- fra %>% filter(SampleID == 'RefSeq' & oldName != '') %>% arrange(generelstart) %>% distinct(oldName)
oneNodeSamples = fra %>% group_by(SampleID) %>% distinct(nodeID) %>% mutate(Nodes = n()) %>%
             filter(Nodes == 1) %>% distinct(SampleID) %>% select(SampleID)
oneNodeSamples = unlist(oneNodeSamples)
ggplot(fra %>% filter((SampleID %in% oneNodeSamples & VariantFreq >= 0.01) | (VariantName == "RefSeq variant")) %>% distinct(VariantName, RefBestHit, generelstart, .keep_all = T),aes(xmin=generelstart, xmax=generelend, y=VariantName,label=oldName,
                   forward = genestrand,fill=oldName)) +
  scale_fill_manual(values=colorBT0436,breaks=c("BT0436",'BT0437','BT0439', 'BT0440'),labels=c('MFS transporter','AGE epimerase','susC','susD') ) +
  geom_gene_arrow(arrowhead_width = unit(0.1,'cm'), arrow_body_height = unit(0.5,'cm')) +
  #geom_gene_arrow(data=fra[fragName != ''],arrowhead_width = unit(0.0,'cm'),arrow_body_height = unit(0.5,'cm'),
  #                aes(xmin=fragrelstart, xmax=fragrelend, fill=fragorigin, color=diftype),alpha=0.7, size=1) +
  scale_color_manual(values=colorDiftype, breaks=c("INS",'DEL'), labels = c('insertion','deletion')) +
  geom_fit_text(size=10, fontface = 'bold.italic') +
  #geom_feature(data= fra %>% filter(str_detect(motif, 'olga')) ,aes(x = motifrelstart, y =VariantName , forward = motifori),feature_height = unit(4,'mm'), feature_width = unit(4,'mm')) +
  #geom_text_repel(data=fra %>% distinct(SampleID, fragName, .keep_all = T),aes(x=fragrelstart, y=VariantName, label=fragName), angle=00,size=2.5, nudge_y =-0.2,segment.colour = NA) +theme_genes() +
  #geom_text_repel(data = fra %>% filter(str_detect(motif, 'olga')), aes(x=motifrelstart,y=VariantName, label=motif), size=3, nudge_y=0.3, angle=20, segment.colour =NA) +
  scale_x_continuous(breaks=seq(0,max(fra[!is.na(generelend)]$generelend), 500)) +
  theme_genes() + theme(axis.text.x = element_text(angle = 50,hjust=1, size=8),
        legend.position = 'right', aspect.ratio = 0.5) +
  labs(title = paste0(shufflon," gene cluster (",first(Ngenes)," - ",last(Ngenes),")"), y='Variant', fill='', x= "Relative position")
ggsave(paste0("shufflons/",shufflon,".svg"), height=7.161, width=13.357)


##### BT3477
shufflon = shufflons[5]
fra = fread(paste0(shufflon,'_fragments_motifs_rel_coords.csv'))
fra <- fra %>% mutate(NodeCount = max(nodeN), .by=SampleID)
fra <- fra %>% mutate(genestrand = ifelse(genestrand == "-", 0,1))
Ngenes <- fra %>% filter(SampleID == 'RefSeq' & oldName != '') %>% arrange(generelstart) %>% distinct(oldName)
oneNodeSamples = fra %>% group_by(SampleID) %>% distinct(nodeID) %>% mutate(Nodes = n()) %>%
             filter(Nodes == 1) %>% distinct(SampleID) %>% select(SampleID)
oneNodeSamples = unlist(oneNodeSamples)
# twoNodeSamples = fra %>% group_by(SampleID) %>% distinct(nodeID) %>% mutate(Nodes = n()) %>%
#              filter(Nodes <= 2) %>% distinct(SampleID) %>% select(SampleID)
# twoNodeSamples = unlist(twoNodeSamples)

ggplot(fra %>% filter((SampleID %in% oneNodeSamples & VariantFreq >= 0.01) | (VariantName == "RefSeq variant")) %>% distinct(VariantName, RefBestHit, generelstart, .keep_all = T),aes(xmin=generelstart, xmax=generelend, y=VariantName,label=oldName,
                   forward = genestrand,fill=oldName)) +
  scale_fill_manual(values=colorBT3477,breaks=c("BT3483",'BT3484',"BT3478", "BT3479"),labels=c('susC','susD', 'site-specific integrase','transposase') ) +
  geom_gene_arrow(arrowhead_width = unit(0.1,'cm'), arrow_body_height = unit(0.5,'cm')) +
  #geom_gene_arrow(data=fra[fragName != ''],arrowhead_width = unit(0.0,'cm'),arrow_body_height = unit(0.5,'cm'),
  #                aes(xmin=fragrelstart, xmax=fragrelend, fill=fragorigin, color=diftype),alpha=0.7, size=1) +
  scale_color_manual(values=colorDiftype, breaks=c("INS",'DEL'), labels = c('insertion','deletion')) +
  geom_fit_text(size=10, fontface = 'bold.italic') +
  #geom_feature(data= fra %>% filter(str_detect(motif, 'olga')) ,aes(x = motifrelstart, y =VariantName , forward = motifori),feature_height = unit(4,'mm'), feature_width = unit(4,'mm')) +
  #geom_text_repel(data=fra %>% distinct(SampleID, fragName, .keep_all = T),aes(x=fragrelstart, y=VariantName, label=fragName), angle=00,size=2.5, nudge_y =-0.2,segment.colour = NA) +theme_genes() +
  #geom_text_repel(data = fra %>% filter(str_detect(motif, 'olga')), aes(x=motifrelstart,y=VariantName, label=motif), size=3, nudge_y=0.3, angle=20, segment.colour =NA) +
  scale_x_continuous(breaks=seq(0,max(fra[!is.na(generelend)]$generelend), 500)) +
  theme_genes() + theme(axis.text.x = element_text(angle = 50,hjust=1, size=8),
        legend.position = 'right', aspect.ratio = 0.5) +
  labs(title = paste0(shufflon," gene cluster (",first(Ngenes)," - ",last(Ngenes),")"), y='Variant', fill='', x= "Relative position")
ggsave(paste0("shufflons/",shufflon,".svg"), height=7.161, width=13.357)

##### BT4520
shufflon = shufflons[1]
fra = fread(paste0(shufflon,'_fragments_motifs_rel_coords.csv'))
fra <- fra %>% mutate(NodeCount = max(nodeN), .by=SampleID)
fra <- fra %>% mutate(genestrand = ifelse(genestrand == "-", 0,1))
Ngenes <- fra %>% filter(SampleID == 'RefSeq' & oldName != '') %>% arrange(generelstart) %>% distinct(oldName)
oneNodeSamples = fra %>% group_by(SampleID) %>% distinct(nodeID) %>% mutate(Nodes = n()) %>%
             filter(Nodes == 1) %>% distinct(SampleID) %>% select(SampleID)
oneNodeSamples = unlist(oneNodeSamples)
# twoNodeSamples = fra %>% group_by(SampleID) %>% distinct(nodeID) %>% mutate(Nodes = n()) %>%
#              filter(Nodes <= 2) %>% distinct(SampleID) %>% select(SampleID)
# twoNodeSamples = unlist(twoNodeSamples)

ggplot(fra %>% filter((SampleID %in% oneNodeSamples & VariantFreq >= 0.01) | (VariantName == "RefSeq variant")) %>% distinct(VariantName, RefBestHit, generelstart, .keep_all = T),aes(xmin=generelstart, xmax=generelend, y=VariantName,label=oldName,
                   forward = genestrand,fill=oldName)) +
  scale_fill_manual(values=colorBT4520,breaks=c("BT4520", "BT4521"),labels=c('restriction endonuclease','site-specific integrase') ) +
  geom_gene_arrow(arrowhead_width = unit(0.1,'cm'), arrow_body_height = unit(0.5,'cm')) +
  #geom_gene_arrow(data=fra[fragName != ''],arrowhead_width = unit(0.0,'cm'),arrow_body_height = unit(0.5,'cm'),
  #                aes(xmin=fragrelstart, xmax=fragrelend, fill=fragorigin, color=diftype),alpha=0.7, size=1) +
  scale_color_manual(values=colorDiftype, breaks=c("INS",'DEL'), labels = c('insertion','deletion')) +
  geom_fit_text(size=10, fontface = 'bold.italic') +
  #geom_feature(data= fra %>% filter(str_detect(motif, 'olga')) ,aes(x = motifrelstart, y =VariantName , forward = motifori),feature_height = unit(4,'mm'), feature_width = unit(4,'mm')) +
  #geom_text_repel(data=fra %>% distinct(SampleID, fragName, .keep_all = T),aes(x=fragrelstart, y=VariantName, label=fragName), angle=00,size=2.5, nudge_y =-0.2,segment.colour = NA) +theme_genes() +
  #geom_text_repel(data = fra %>% filter(str_detect(motif, 'olga')), aes(x=motifrelstart,y=VariantName, label=motif), size=3, nudge_y=0.3, angle=20, segment.colour =NA) +
  scale_x_continuous(breaks=seq(0,max(fra[!is.na(generelend)]$generelend), 500)) +
  theme_genes() + theme(axis.text.x = element_text(angle = 50,hjust=1, size=8),
        legend.position = 'right', aspect.ratio = 0.5) +
  labs(title = paste0(shufflon," gene cluster (",first(Ngenes)," - ",last(Ngenes),")"), y='Variant', fill='', x= "Relative position")
ggsave(paste0("shufflons/",shufflon,".svg"), height=7.161, width=13.357)