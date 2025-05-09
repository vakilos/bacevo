# bacevo
Scripts and workflow for the analysis of experimental evolution datasets.

This code has been developed for the data analysis in the following publications:
<br>&nbsp;&nbsp; i. _Nutrient landscape shapes adaptive diversification of the prominent human
gut commensal Bacteroides thetaiotaomicron_ 
<br>&nbsp;&nbsp; ii. _Rapid genetic diversification of Bacteroides thetaiotaomicron in mono-associated mice 
revealed through deep population-level sequencing_

<br>

## Polymorphism detection
#### Preparation
Download RefSeq genome and annotation files, for example here for B. theta VPI:
[RefSeq link](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Bacteroides_thetaiotaomicron/all_assembly_versions/GCF_000011065.1_ASM1106v1/) <br>

Prepare a _csv_ file that contains samples meta-data and add it in the project directory. It is expected that this files has the columns below:
- &nbsp;&nbsp;**SampleID** (sequencing library ID)
- &nbsp;&nbsp;**MetagenomeID** (same with SampleID for metagenomes - ID of the respective metagenomic sample for isolates)
- &nbsp;&nbsp;**Isolate** (ID for isolate if any)
- &nbsp;&nbsp;**Mouse** or **Replicate** (Mouse for _in vivo_ experiments, Replicate for _in vitro_)
- &nbsp;&nbsp;**Day** (days post inoculation)
- &nbsp;&nbsp;**Cage** or **Mix** (Cage for _in vivo_ experiments, Mix for _in vitro_)
- &nbsp;&nbsp;**SampleType** (label the type of sample)

#### Execution
1. Run [`ep_annotable.py`](./scripts/ep_annotable.py) to create a table that stores RefSeq's annotation data.
2. Run [`ep_start.py`](./scripts/ep_start.py) in a SLURM configured system to submit an array job comprising 
multiple instances of the main pipeline which preprocesses sequencing reads and calls the breseq pipeline for a given number of samples.  
3. Run [`ep_align_stats.py`](./scripts/ep_align_stats.py) to create reports of alignment statistics (parsed from _**Breseq**_'s output).
4. Run [`ep_parse.py`](./scripts/ep_parse.py) to combine **_Breseq_** output with meta-data. (exports  _<name\>\_AnnotatedPolyTable.txt_)
5. Apply additional filters on polymorphisms and calculate metrics/stats 
by running [`ep_lyse_invivo.py`](./scripts/ep_lyse_invivo.py) for an _in vivo_ 
or [`ep_lyse_invitro.py`](./scripts/ep_lyse_invitro.py) for 
an _in vitro_ context respectively.  (exports  _<name\>\_FinalPolyTable.tsv_)
6. Merge multiple _<name\>\_FinalPolyTable.tsv_ tables
and estimate prevalence stats in the concatenated table, by running [`ep_merge_tables.py`](./scripts/ep_merge_tables.py)

## Phase variation detection
This pipeline evaluates the presence of phase variants in genome assemblies. Genome assemblies and annotation should be organised in directories per sample. 

Pipeline Execution order:  

1. [`ai_start.py`](./scripts/ai_start.py)
2. [`ai_reverse_nodes.py`](./scripts/reverse_nodes.py) 
3. [`ai_gene_order.py`](./scripts/ai_gene_order.py)
4. [`ai_gene_alignment.py`](./scripts/ai_gene_alignment.py)

## _dNdS_ analysis
This analysis relies on alignment files after mapping reads to a reference genome (BAM files).
Files should be organised in directories per sample. The first step of the pipeline will export files in a created (if not present) _dnds_ directory. Files in this directory will need to be removed if need to be replaced or run analysis again, as they are NOT overwritten. Then, `dnds_concat.py` will concatenate all sample _.dnds_ files and a concatenated table will be exported by default at the parent directory.

Pipeline Execution order:  

1. [`dnds_from_table_invivo.py`](./scripts/dnds_from_table_invivo.py) or [`dnds_from_table_invitro.py`](./scripts/dnds_from_table_invitro.py)
2. [`dnds_concat.py`](./scripts/dnds_concat.py)
3. [`dnds_figure_invivo.py`](./scripts/dnds_figure_invivo.py) or [`dnds_figure_invitro.py`](./scripts/dnds_figure_invitro.py)



## Scripts for analysis specific to the manuscripts:
### (i) _in vitro_:
[`ep_figure_invitro.py`](./scripts/ep_figure_invitro.py) - exports various tables for plotting in R.
<br>[`plot_invitro_figures.R`](./scripts/plot_invitro_figures.R) - creates plots 


### (ii) _in vivo_:
[`ep_figures_invivo.py`] calls `ep_enrichment.py` (Enrichment analysis), `ep_turnover.py` (Polymorphism turnover and dynamics) and `ep_figure_pcoa.py` (PCoA on Sampletypes)

[`ep_enrichment.py`](./scripts/ep_enrichment.py)
<br>[`ep_turnover.py`](./scripts/ep_turnvover.py)
<br>[`ep_figure_invivo.py`](./scripts/ep_figure_invivo.py)
<br>[`ep_figure_pcoa.py`](./scripts/ep_figure_pcoa.py)

#### **Intestinal compartments**

Pipeline Execution order:  

1. [`ep_comp_significant.py`](./scripts/ep_comp_significant.py)
2. [`Maaslin2.R`](./scripts/Maaslin2.R)
3. [`ep_comp_select.py`](./scripts/ep_comp_select.py)

#### **Poly clusters and shufflon signatures**
We selected polymorphisms present in 80% of sampled timepoints
for each mouse, which had read support of at least 100 reads. We then calculated Pearson Correlation
of frequencies overtime for each pair of polymorphisms for each mouse (missing values were represented as zeros).
For highly correlated pairs r > 0.8, we created clusters in a step process which maximized mean correlation within cluster at each step,
recursively, until no polymorphisms remained unclustered. Then, shufflon signatures were defined as 
polymorphisms under the same cluster within 1000 bp distance on the genome. This was a primer for manual inspection of isolate genome assemblies. All polymorphisms within shuffling regions that were also validated from isolate genomes (Phase variation pipeline)
were masked from the rest of the analysis - as it was hard to distinguish _de novo_ mutations and false positive calls that occured due to inversions/ recombinations in these regions.

Pipeline Execution order:  

1. [`pc_correlations.py`](./scripts/pc_correlations.py)
2. [`pc_signatures.py`](./scripts/pc_signatures.py)
3. [`pc_shufflons.py`](./scripts/pc_shufflons.py)

#### Plotting 
[`Figure2_plots.R`](./scripts/Figure2_plots.R)
[`Figure3_plots.R`](./scripts/Figure3_plots.R)
[`Figure4_plots.R`](./scripts/Figure4_plots.R)
[`Figure5_plots.R`](./scripts/Figure5_plots.R)


