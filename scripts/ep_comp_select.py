import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
from PIL.ImagePalette import sepia

from bacevo_fun import *

"""
Adds freqSI, freqLI, preSI, preLI (frequency and prevalence values for each poly in the two intestinal compartments
"""

os.chdir("/home/christos/bacevo/invivo")
console = Console(color_system='windows')

#import polars as pl
#mask = pl.read_csv("tmp/confirmed_invivo.csv").select(pl.col("PolyID")).unique() ### this is polymorphisms that have been masked from shuffling regions..
mask = pd.read_csv("tmp/confirmed_invivo.csv")['PolyID'].unique()
#sign = pl.read_csv("tmp/compartments_stat_sign_polys.csv").filter((pl.col('qval') <= 0.05) & (pl.col('padj') <= 0.05)).select(pl.col("PolyID")).unique()

## choose between _readcount.csv and _chromreads.csv

# sign = sign.loc[(sign.qval <= 0.05) & (sign.padj <= 0.05)]['PolyID'].unique()
#df= pl.read_csv("invivo_metagenomes_FinalPolyTable.tsv", separator='\t')
# q = (
# pl.scan_csv("tmp/invivo_metagenomes_FinalPolyTable_comp28.tsv", separator='\t'))
# df = q.collect()
df = pd.read_csv('tmp/invivo_metagenomes_FinalPolyTable_comp28.tsv', sep='\t').drop_duplicates(['SampleID','PolyID'])
#sign = pd.read_csv("tmp/compartments_stat_sign_polys_chromreads.csv")
sign = pd.read_csv("tmp/compartments_stat_sign_polys_readcount.csv")

#df = df.loc[df.PolyID.isin(sign.PolyID.unique())]
# sign = pd.read_csv("tmp/compartments_stat_sign_polys_readcount.csv")

# df = pd.read_csv("invivo_metagenomes_FinalPolyTable.tsv", sep='\t')
# keep_samples = df.loc[(df.Chrom == 'NC_004663') & (df['Coverage_mean'] >= 100), 'SampleID'].unique()
# df = df.loc[(~df.PolyID.isin(mask)) & (df.Ancestry == 'evolved') & (df.SampleID.isin(keep_samples))].drop_duplicates(['SampleID','PolyID'])
# df['SampleTypeGroup'] = df.SampleType.apply(lambda x:'SI' if x in ["SIP",'SID'] else 'LI')
#df = df.with_columns(df.select(pl.all(), pl.when(pl.col("SampleType").is_in(['SIP','SID'])).then(pl.lit('SI')).otherwise(pl.lit('LI')).alias('SampleTypeGroup')))

#sample_type_size = df.unique(subset='SampleID').group_by("SampleType").agg(pl.col('PolyID').len()).rename({'PolyID':'size'})
#sample_type_group_size = df.unique(subset='SampleID').groupby("SampleTypeGroup").agg({'PolyID':'size'}).rename({'PolyID':'size'})
sampletypegroup_sizes = df.drop_duplicates('SampleID').groupby("SampleTypeGroup").agg({'PolyID':'size'})['PolyID'].to_dict()
sampletype_sizes = df.drop_duplicates("SampleID").groupby('SampleType').agg({'PolyID':'size'})['PolyID'].to_dict()
poly_catalog = df.PolyID.unique()
tested_polys = sign.PolyID.unique()
samples = df.SampleID.unique()
pre = pd.DataFrame(index = poly_catalog, columns=["freqSI",'freqLI','preSI','preLI'], data = np.zeros(shape=[len(poly_catalog), 4]))

console.print('Parsing frequency and prevalence data...')
for name,group in df.groupby('PolyID'):
   meanSI = group.loc[group.SampleTypeGroup == 'SI']['Freq'].mean()
   meanLI = group.loc[group.SampleTypeGroup == 'LI']['Freq'].mean()
   preSI = group.loc[group.SampleTypeGroup == 'SI'].shape[0] / sampletypegroup_sizes['SI']
   preLI = group.loc[group.SampleTypeGroup == 'LI'].shape[0] / sampletypegroup_sizes['LI']
   pre.loc[name, ['freqSI','freqLI','preSI','preLI']] = [meanSI, meanLI, preSI, preLI]

pre = pre.reset_index().rename(columns={'index':'PolyID'}).fillna(0)
pre = df.drop_duplicates("PolyID")[['PolyID','GeneName','oldGeneName', 'Chrom', 'Synonymy']].merge(pre, on='PolyID', how='right')
print(pre.isnull().values.any())
### This merge adds all not significant (from Maaslin2) polys
out = sign.merge(pre, on=['PolyID','GeneName'], how='right')


console.print('Adding columns...')
out[['padj','qval']] = out[['padj','qval']].fillna(1)
print(out[['freqSI','freqLI']].isnull().values.any())
print(out[out.freqSI != out.freqLI].shape)
out['Diff'] = out.meanSI - out.meanLI
out['Diffreq'] = out.freqSI - out.freqLI

out['Change'] = out.PolyID.apply(lambda x:"→".join(x.split(":")[2:]))
out['Chromosome'] = out.PolyID.apply(lambda x:'c'+x.split(":")[1] if x.split(":")[0] == 'NC_004663' else 'p')
out['PolyLabel'] = out.PolyID.apply(lambda x:'c'+x.split(":")[1]+"  "+"→".join(x.split(":")[2:]) if x.split(":")[0] == 'NC_004663' else 'p'+x.split(":")[1]+":"+"→".join(x.split(":")[2:]))

out['ColorGroup'] = out.apply(lambda x: 'Mann-Whitney' if x.padj <= 0.05 < x.qval else "Maaslin2" if x.qval <= 0.05 < x.padj else "Consensus" if x.qval <= 0.05 and x.padj <= 0.05 else 'Other', axis=1)
print(out.ColorGroup.unique())
out['Appearance'] = out.apply(lambda x: "LI" if x.preSI == 0 and x.preLI > 0 else "SI" if x.preSI > 0 and x.preLI == 0 else 'LI+SI', axis=1)
print(out[out.freqSI != out.freqLI].shape)

console.print('Estimating poly distances...')
"""
Estimate distance between polys (for the significant ones) 
"""
distance = pd.DataFrame()
for name, group in out.loc[(out.ColorGroup == 'Consensus')].groupby('GeneName'):
   pairwise_differences = [abs(i - j) for i, j in itertools.combinations(group.Pos, 2)]
   distance = pd.concat([distance, pd.DataFrame({'GeneName':[name for x in pairwise_differences],
                                                             "BPDistance":pairwise_differences})])

distance.to_csv('tmp/Comp_significant_polys_BPdistance_by_gene.csv',index=False)
out.to_csv('tmp/Compartment_polys_stats.csv', index=False)
console.print("THIS IS THE END!! ")
sys.exit()

























# for st in sample_type_group_size.iter_rows(named=True):
#    #print(st['SampleType'], st['size'])
#    for _,poly in zip(track(range(len(all_polys)), description=f'Working on {st["SampleTypeGroup"]}'),all_polys):
#       count = df.unique(subset=['SampleID',"PolyID"]).filter((pl.col("PolyID") == poly) & (pl.col('SampleTypeGroup') == st['SampleTypeGroup'])).shape[0]
#       prev = count / st['size']
#       #print(count, prev, st['size'])
#       st_mean = df.unique(['SampleID',"PolyID"]).filter((pl.col('PolyID') == poly) & (pl.col('SampleTypeGroup') == st['SampleTypeGroup'])).select(pl.mean('Freq').fill_null(0)).to_numpy()[0]
#       #print(st_mean)
#       #print(pl.DataFrame({"PolyID":poly, "SampleType":[st['SampleType']],'SampleTypePrev':[prev], "SampleTypeMean":st_mean}))
#       pre = pl.concat([pre, pl.DataFrame({"PolyID":poly, "SampleTypeGroup":[st['SampleTypeGroup']],'SampleTypeGroupPrev':[prev], "SampleTypeGroupMean":st_mean})])
# print(pre)

df = pd.read_csv('~/bacevo/invivo/tmp/SampleTypeGroupPrev.csv')
prev = pre.pivot(index='PolyID', columns='SampleTypeGroup', values='SampleTypeGroupPrev').rename(columns={"SI":"preSI",'LI':'preLI'}).reset_index()
fre = pre.pivot(index='PolyID', columns='SampleTypeGroup', values='SampleTypeGroupMean').rename(columns={"SI":"freqSI",'LI':'freqLI'}).reset_index()
prev = prev.merge(fre, on='PolyID', how='left')
prev = df.drop_duplicates("PolyID")[['PolyID','GeneName','oldGeneName','Pos']].merge(prev, on='PolyID', how='right')
prev.write_csv('tmp/SampleTypeGroupPrev.csv')
sys.exit()







for st in sample_type_size.iter_rows(named=True):
   #print(st['SampleType'], st['size'])
   for _,poly in zip(track(range(len(all_polys)), description=f'Working on {st["SampleType"]}'),all_polys):
      count = df.unique(subset=['SampleID',"PolyID"]).filter((pl.col("PolyID") == poly) & (pl.col('SampleType') == st['SampleType'])).shape[0]
      prev = count / st['size']
      #print(count, prev, st['size'])
      st_mean = df.unique(['SampleID',"PolyID"]).filter((pl.col('PolyID') == poly) & (pl.col('SampleType') == st['SampleType'])).select(pl.mean('Freq').fill_null(0)).to_numpy()[0]
      #print(st_mean)
      #print(pl.DataFrame({"PolyID":poly, "SampleType":[st['SampleType']],'SampleTypePrev':[prev], "SampleTypeMean":st_mean}))
      pre = pl.concat([pre, pl.DataFrame({"PolyID":poly, "SampleType":[st['SampleType']],'SampleTypePrev':[prev], "SampleTypeMean":st_mean})])
print(pre)
pre.write_csv('tmp/SampleTypePrev.csv')

sys.exit()

# df = pd.read_csv("invivo_metagenomes_FinalPolyTable.tsv", sep='\t')
# #df = df[~df.SampleID.isin([50915, 50917, 50919, 50918, 50913])]
# df = df.loc[df['Coverage_mean'] >= 100]
#
# df = df[~df.PolyID.isin(mask.PolyID.unique())]
# df = df.loc[(df.Dataset == 'ploen') & (df.Day == 28)]
#
# sample_types = df.SampleType.unique()
# print(df.drop_duplicates("SampleID")['SampleType'].value_counts())
### Find polys enriched in certain compartments.
stdf = pd.DataFrame()

#table = pd.pivot(data=df, index='SampleID', columns='PolyID', values="Freq").fillna(0)
#table[].to_csv("tmp/compartments_table1.csv", index=False)
# for _,(name,group) in zip(track(range(len(df.groupby("PolyID").groups))),df.groupby('PolyID')):
#     for st in sample_types:
#         if not group.loc[group.SampleType == st].empty:
#             st_mean_freq = group.loc[group.SampleType == st, "Freq"].mean()
#             stdf = pd.concat([stdf, pd.DataFrame({"PolyID": [name], 'SampleType': [st], 'MeanFreq': [st_mean_freq]})])
#         else:
#             stdf = pd.concat([stdf, pd.DataFrame({"PolyID": [name], 'SampleType': [st], 'MeanFreq': [0]})])
# print(stdf)

# for name,group in df[['SampleID','PolyID','Day','Mouse']].groupby('SampleID'):
#     for st in sample_types:
#         if group.loc[group.SampleType == st].empty:
#             df = pd.concat([df, pd.DataFrame({})])


pcoa_columns = ["SampleID",'Mouse', 'Day','Cage','MetagenomeID','Dataset', "SampleType"]


#make_pcoa(df, binary=True, plot=False, color="SampleType", style='Day', keep_columns= pcoa_columns, NameTag="invivo", metric='braycurtis')

#importance = make_random_forest_model(df, variable='SampleType', binary=True, oneHotEncode=True)
stg = {"SI":["SIP"]}
df['SampleTypeGroup'] = df['SampleType'].apply(lambda x:"SI" if x in ["SIP","SID"] else "LI")

df['SampleID'] = df['SampleID'].astype(str).apply(lambda x: "S"+x)


""" Prevalence of Polys per SampleType and SampleTypeGroup"""
pre = pd.DataFrame()

sample_types_size = df.drop_duplicates("SampleID").groupby('SampleType').agg({'SampleID':'size'}).to_dict()['SampleID']
sample_type_groups_size = df.drop_duplicates("SampleID").groupby('SampleTypeGroup').agg({'SampleID':'size'}).to_dict()['SampleID']

all_polys = df.PolyID.unique()
print(f"Number of total polys to be analyzed: {len(all_polys)}")

for st in sample_types_size:
   for _,poly in zip(track(range(len(all_polys)), description=f'Working on {st}'),all_polys):
      count = df.drop_duplicates(['SampleID',"PolyID"]).loc[(df.PolyID == poly) & (df.SampleType == st)].shape[0]
      prev = count / sample_types_size[st]
      st_mean = df.drop_duplicates(['SampleID',"PolyID"]).loc[(df.PolyID == poly) & (df.SampleType == st),'Freq'].mean
      pre = pd.concat([pre, pd.DataFrame({"PolyID":[poly], "SampleType":[st],'SampleTypePrev':[prev], "SampleTypeMean":[st_mean]})])
#    for poly, subgroup in group.groupby("PolyID"):
#       prev = round(subgroup.shape[0] /size, 4)
#       st_mean = np.mean(subgroup.Freq.tolist() + [0]*(size - subgroup.shape[0]))
#
#       ### add zeros based on 1 - prev (fraction of absences)
#       pre = pd.concat([pre, pd.DataFrame({"PolyID":[poly], "SampleType":[st],'SampleTypePrev':[prev], "SampleTypeMean":[st_mean]})])
pre.to_csv("tmp/SampleTypePrev.csv", index=False)

pre = pd.DataFrame()
for st in sample_type_groups_size:
   for poly in all_polys:
      count = df.drop_duplicates(['SampleID',"PolyID"]).loc[(df.PolyID == poly) & (df.SampleTypeGroup == st)].shape[0]
      prev = count / sample_type_groups_size[st]
      st_mean = df.drop_duplicates(['SampleID',"PolyID"]).loc[(df.PolyID == poly) & (df.SampleTypeGroup == st),'Freq'].mean
      pre = pd.concat([pre, pd.DataFrame({"PolyID":[poly], "SampleTypeGroup":[st],'SampleTypeGroupPrev':[prev], "SampleTypeGroupMean":[st_mean]})])
#    for poly, subgroup in group.groupby("PolyID"):
#       prev = round(subgroup.shape[0] /size, 4)
#       st_mean = np.mean(subgroup.Freq.tolist() + [0]*(size - subgroup.shape[0]))
#
#       ### add zeros based on 1 - prev (fraction of absences)
#       pre = pd.concat([pre, pd.DataFrame({"PolyID":[poly], "SampleType":[st],'SampleTypePrev':[prev], "SampleTypeMean":[st_mean]})])
pre.to_csv("tmp/SampleTypePrev.csv", index=False)

sns.scatterplot(pre, x='SampleTypeGroupPrev', y='SampleTypeGroupMean')
plt.show()

sys.exit()

for st, group in df.drop_duplicates(['SampleID','PolyID']).groupby('SampleTypeGroup'):
   size = group.shape[0]
   for poly, subgroup in group.groupby("PolyID"):
      prev = round(subgroup.shape[0] / size, 4)
      st_mean = np.mean(subgroup.Freq.tolist() + [0] * (size - subgroup.shape[0]))

      ### add zeros based on 1 - prev (fraction of absences)

      pre = pd.concat([pre, pd.DataFrame({"PolyID": [poly], 'SampleTypeGroup':[st], 'SampleTypeGroupPrev': [prev], "SampleTypeGroupMean": [st_mean]})])
pre.to_csv("tmp/SampleTypeGroupPrev.csv", index=False)
sys.exit()





