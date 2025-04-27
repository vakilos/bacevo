import os
import sys

import pandas as pd

from bacevo_fun import  *

parser = argparse.ArgumentParser(description='Creates a mask file with polymorphisms within shufflons and also determines the activity of shufflons in metagenomes..')
parser.add_argument('-m', dest='mode', choices=['invivo','invitro'], default='invivo', help='Default=invivo')
parser.add_argument('-d', dest='basedir', help='Project directory. Default = current working directory', default=os.getcwd())
parser.add_argument('-a', dest='anno', help='Path to AnnotationTable.txt (see ep_annotable.py)', required=True)

args = parser.parse_args()
mode = args.mode
os.chdir(args.basedir)
console = Console()

df = pd.read_csv("tmp/my_clusters_stats.csv")

import polars as pl
q = (
pl.scan_csv("invivo_metagenomes_FinalPolyTable.tsv", separator='\t').filter((pl.col("Coverage_mean") >= 100)  & (pl.col("Dataset") == 'ploen') & (pl.col("Day") == 28)))
com = q.collect()
com = com.with_columns(com.select(pl.all(), pl.when(pl.col("SampleType").is_in(['SIP','SID'])).then(pl.lit('SI')).otherwise(pl.lit('LI')).alias('SampleTypeGroup')))
com = com.to_pandas()

ano = pd.read_csv(args.anno, sep='\t')

### list of genes from metagenomic evidence for shufflons
me_genes = df.loc[(df.Cluster != 'no cluster') & (df.Signature == True)].drop_duplicates("GeneName")[["GeneName", 'oldGeneName']].sort_values('GeneName')

shufflons = df.copy()
target = pd.read_csv("gene_order/target_genes.csv")


reported = ['BT1042','BT2268'] ## those reported previously
confirmed = ["BT1042",'BT2268',"BT3239", "BT4540", "BT4520", "BT0436", 'BT3477'] ## support from both assemblies and metagenomes
putative = ['BT0072','BT0318',"BT1922",'BT3746','BT4778'] ## supported only from metagenomes (possible detection limit on low frequency variants in isolates, due to extreme downsampling.
notshuffling = ['BT0681', 'BT0981','BT1553', 'BT1636','BT2172','BT2240','BT2469', "BT3439", 'BT4029', 'BT4440']
### BT4029 gene is reversed in the ancestral (not a shufflon though), BT4440: partially duplicated, some variants have short versions, some have BT4441 instead
target['novelty'] = ["novel" for x in target.index] ## if it is described for the first time
target['evidence'] = ['no' for x in target.index] ## if it is fully supported



target.loc[target.GeneGroup.isin(reported), 'novelty'] = 'published'
target.loc[target.GeneGroup.isin(confirmed),'evidence'] = 'confirmed'
#print(target.set_index("GeneName").to_dict()['GeneGroup'])
confirmed_genes = list(target.loc[target.evidence == 'confirmed', 'GeneName'].unique())
print(f'\nNumber of confirmed genes: {len(confirmed_genes)}')
shufflons['GeneGroup'] = shufflons.GeneName.map(target.set_index("GeneName").to_dict()['GeneGroup']) ## GeneGroup contains all GeneNames for a shuffling regions
shufflons['ShufflonGroup'] = ["no shufflon" for x in shufflons.index] ## shufflonGroup contains all polys for a shufflon
shufflons_length = {}
shufflons_locations = {}
shufflons_start = {}
shufflons_end = {}
for name, group in target.groupby("GeneGroup"):
    Min = ano.loc[ano['locus_tag'].isin(group.GeneName.values),'start'].astype(int).min()
    Max = ano.loc[ano['locus_tag'].isin(group.GeneName.values),'end'].astype(int).max()
    shufflons_length[name] = Max - Min
    shufflons_locations[name] = round(np.mean([Max,Min]))
    #shufflons_start[name] = round()

polys_within_shufflons = [] ## for each GeneGroup (shufflon name) get min and max coords and include all polys within its range

for name, group in shufflons.loc[shufflons['GeneGroup'].isin(confirmed)].groupby(['GeneGroup','Dataset']): ## group each PolyID in shufflons and update the shufflons dataframe
    genegroup, dataset = name
    Min = group.Pos.min()
    Max = group.Pos.max()
    this_polys = shufflons.loc[(shufflons.Pos.between(Min,Max)) & (shufflons.Dataset == dataset), "PolyID"].tolist()
    this_polys_sign = shufflons.loc[(shufflons.Pos.between(Min,Max)) & (shufflons.Dataset == dataset) & (shufflons.Signature == True),"PolyID"].tolist()

    shufflons.loc[(shufflons.Pos.between(Min,Max)) & (shufflons.Dataset == dataset),'ShufflonGroup'] = genegroup

    polys_within_shufflons += this_polys

    console.print(f"{dataset}: Shufflon {genegroup} || range {Min} - {Max}|: {len(set(this_polys))} "
              f"polys (catalogue)[{len(set(this_polys_sign))/len(set(this_polys)):.2%} with signature ({len(set(this_polys_sign))} polys)] - {len(this_polys)} mutations [{len(this_polys_sign)/len(this_polys):.2%} with signature].")

console.print(f'\nTotal polys (catalogue) within shufflon regions (will be masked): '
      f'{len(set(polys_within_shufflons))}/{shufflons.drop_duplicates("PolyID").shape[0]} - total mutations {len(polys_within_shufflons)}/{shufflons.shape[0]} (* fecal samples)\n')


shufflons = shufflons.loc[shufflons.PolyID.isin(polys_within_shufflons)] ### subset for only polys within shufflon regions


shufflons['ShufflonLength'] = shufflons.ShufflonGroup.map(shufflons_length)
shufflons['ShufflonLocation'] = shufflons.ShufflonGroup.map(shufflons_locations)


shufflons.to_csv('tmp/confirmed_invivo.csv', index=False) ## this to be used for masking


print('\n\nDone!!')
