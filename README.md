# bacevo
Scripts and workflow for the analysis of experimental evolution datasets.

This code has been developed for the data analysis in the following publications:
<br>&nbsp;&nbsp; i. _Nutrient landscape shapes adaptive diversification of the prominent human
gut commensal Bacteroides thetaiotaomicron_ 
<br>&nbsp;&nbsp; ii. _Rapid genetic diversification of Bacteroides thetaiotaomicron in mono-associated mice revealed through deep population-level sequencing_

<br>

## Polymorphism detection
Download RefSeq annotation files, for example here for B. theta VPI:
[RefSeq link](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Bacteroides_thetaiotaomicron/all_assembly_versions/GCF_000011065.1_ASM1106v1/) <br>

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
This pipeline identifies phase-variable regions in genome assemblies.

<span style="color:red">Pipeline Execution order</span>  

1. [`ai_start.py`](./scripts/ai_start.py)
2. [`ai_reverse_nodes.py`](./scripts/reverse_nodes.py) 
3. [`ai_gene_order.py`](./scripts/ai_gene_order.py)
4. [`ai_gene_alignment.py`](./scripts/ai_gene_alignment.py)

## _dNdS_ analysis

<span style="color:red">Pipeline Execution order</span>  

1. [`dnds_from_table_invitro.py`](./scripts/dnds_from_table_invitro.py)
2. [`dnds_concat.py`](./scripts/dnds_concat.py)
3. [`dnds_figure_invivo.py`](./scripts/dnds_figure_invivo.py) or [`dnds_figure_invitro.py`](./scripts/dnds_figure_invitro.py)



## Scripts for analysis specific to the manuscripts:
### (i) _in vitro_:
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - [`ep_figure_invitro.py`](./scripts/ep_figure_invitro.py) - exports various tables for plotting in R.
<br> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - [`plot_invitro_figures.R`](./scripts/plot_invitro_figures.R) - creates plots 


### (ii) _in vivo_:
[`ep_figures_invivo.py`] calls `ep_enrichment.py` (Enrichment analysis), `ep_turnover.py` (Polymorphism turnover and dynamics) and `ep_figure_pcoa.py` (PCoA on Sampletypes)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - [`ep_enrichment.py`](./scripts/ep_enrichment.py)
<br> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - [`ep_turnover.py`](./scripts/ep_turnvover.py)
<br> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - [`ep_figure_invivo.py`](./scripts/ep_figure_invivo.py)
<br> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - [`ep_figure_pcoa.py`](./scripts/ep_figure_pcoa.py)

#### **Intestinal compartments**

<span style="color:red">Pipeline Execution order</span>  

1. [`ep_comp_significant.py`](./scripts/ep_comp_significant.py)
2. [`Maaslin2.R`](./scripts/Maaslin2.R)
3. [`ep_comp_select.py`](./scripts/ep_comp_select.py)

#### **Poly clusters and shufflon signatures**
We selected polymorphisms present in 80% of sampled timepoints
for each mouse, which had read support of at least 100 reads. We then calculated Pearson Correlation
of frequencies overtime for each pair of polymorphisms for each mouse (missing values were represented as zeros).
For highly correlated pairs r > 0.8, we created clusters in a step process which maximized mean correlation within cluster at each step,
recursively, until no polymorphisms remained unclustered. Then, shufflon signatures were defined as 
polymorphisms under the same cluster within 1000 bp distance on the genome. All polymorphisms within shuffling regions
were masked for the rest of the analysis - as it is hard to distinguish _de novo_ mutations from false positive calls that appear due to
inversions/ recombinations in these regions.

<span style="color:red">Pipeline Execution order</span>  

1. [`pc_correlations.py`](./scripts/pc_correlations.py)
2. [`pc_signatures.py`](./scripts/pc_signatures.py)
3. [`pc_shufflons.py`](./scripts/pc_shufflons.py)



