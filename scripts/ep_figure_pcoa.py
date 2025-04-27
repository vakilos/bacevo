import sys

import pandas as pd


#### exported files are fed to plot_figures2.R


from bacevo_fun import *

parser = argparse.ArgumentParser(description="Estimate samples distances and runs a PCoA")
parser.add_argument('-i', dest='polytable', help='Path to invivo_FinalPolyTable.tsv', required=True)
parser.add_argument('-m', dest='mask', help='File for masking polymorphisms.', default=None)
args = parser.parse_args()


df = pd.read_csv(args.polytable,sep='\t') ### this is the concatenation of bern_metagenomes and ploen FinalPolyTables with the ep_merge_tables.py script
if args.mask:
    mask = pd.read_csv(args.mask)
    df = df[~df.PolyID.isin(mask.PolyID.unique())]

keep_samples = df.loc[(df.Chrom == 'NC_004663') & (df['Coverage_mean'] >= 100), 'SampleID'].unique()
df = df.loc[(df.TotalCov.ge(100)) &  (df.Ancestry == 'evolved')  & (df.SampleID.isin(keep_samples))].drop_duplicates(['SampleID','PolyID'])

meta = df.drop_duplicates("SampleID")[["SampleID",'Day','Mouse', 'Cage','Dataset', "SampleType",'PolyCount','SampleProfile']]

# add the ancestral libraries (not in the final table)
meta = pd.concat([meta,pd.DataFrame({"Day":[0,0],
                                 'Mouse':[0,0],
                                 "Cage":[0,0],
                                 "Dataset":['bern_metagenomes','ploen'],
                                 "SampleID":[35254,50995], 'SampleProfile':['metagenome','metagenome'],
                                     "PolyCount":[0,0],
                                 "SampleType":["Feces","Feces"]})])

meta.set_index("SampleID",inplace=True)

"""Estimate distance between samples based on Poly presence/frequency."""
from skbio.stats.distance import DistanceMatrix
metric = 'jensenshannon'
metric = 'braycurtis'
#metric = 'euclidean'
make_binary = True


# this dataset USES masked polys.. but need to set ancestral frequencies to zero.
dmd = df.loc[df['SampleType'] == 'Feces']
table = pd.pivot(data=dmd, index='SampleID', columns='PolyID', values="Freq").fillna(0)
## add ancestral samples as zeros..
anc = pd.DataFrame({35254:[0 for x in range(table.shape[1])],50995: [0 for x in range(table.shape[1])]}).T.rename(columns={x:y for x,y in zip(range(table.shape[1]),table.columns.values)})
#anc = pd.DataFrame({35254:[0 for x in range(table.shape[1])]}).T.rename(columns={x:y for x,y in zip(range(table.shape[1]),table.columns.values)})

table = pd.concat([table,anc])
binary_table = table.copy()
if make_binary:
    binary_table[binary_table > 0] = 1

di = DistanceMatrix(scipy.spatial.distance.pdist(binary_table.values, metric=metric)).data
dm = pd.DataFrame(index=binary_table.index.values, columns=binary_table.index.values, data=di)
dm = dm.stack().reset_index().rename(columns={'level_0': 's1', 'level_1': 's2', 0: 'Distance'})

#di = di.merge(df.drop_duplicates("SampleID")[["SampleID","Day","Mouse",'Dataset','SampleType']], on='SampleID', how='left')
time_dict = {0:0,3:1, 7:2, 14:3, 21:4, 28:5}
dm["DayA"] = [ meta.loc[x,"Day"] for x in dm['s1'].astype(int)]
dm["DayB"] = [ meta.loc[x,"Day"] for x in dm['s2'].astype(int)]
dm['TimepointA'] = dm.DayA.map(time_dict)
dm['TimepointB'] = dm.DayB.map(time_dict)
dm['TimepointLatestSample'] = dm.apply(lambda x: x.s1 if x.TimepointA > x.TimepointB else x.s2, axis=1)
dm['TimepointLatest'] = dm.apply(lambda x: x.TimepointA if x.TimepointA > x.TimepointB else x.TimepointB, axis=1)
dm['MouseA'] = [ meta.loc[x,"Mouse"] for x in dm['s1']]
dm['MouseB'] = [ meta.loc[x,"Mouse"] for x in dm['s2']]
dm['CageA'] = [meta.loc[x, "Cage"] for x in dm['s1']]
dm['CageB'] = [meta.loc[x, "Cage"] for x in dm['s2']]
dm['SampleTypeA'] = [meta.loc[x,"SampleType"] for x in dm['s1']]
dm['SampleTypeB'] = [meta.loc[x,"SampleType"] for x in dm['s2']]
dm["DatasetA"] = [meta.loc[x, "Dataset"] for x in dm["s1"]]
dm["DatasetB"] = [meta.loc[x, "Dataset"] for x in dm["s2"]]
dm['DatasetGroup'] = dm.filter(regex='Dataset', axis='columns').apply(lambda x: 1 if x.DatasetA == x.DatasetB else 0, axis=1)
dm['SampleTypeGroup'] = dm.filter(regex="SampleType",axis='columns').apply(lambda x: 1 if x.SampleTypeA == x.SampleTypeB else 0, axis=1)
dm['DayGroup'] = dm.filter(regex="Day",axis='columns').apply(lambda x: 1 if x.DayA == x.DayB else 0, axis=1)
dm['DayDistance'] = dm.filter(regex='Day', axis='columns').apply(lambda x: abs(x.DayA - x.DayB), axis=1)
dm['TimepointDistance'] = (dm.TimepointA - dm.TimepointB).abs()
dm['TimepointInterval'] = dm.apply(lambda x:min(x.TimepointA, x.TimepointB),axis=1)
dm['MouseGroup'] = dm.filter(regex="Mouse",axis='columns').apply(lambda x: 1 if x.MouseA == x.MouseB else 0 , axis=1)
dm['CageGroup'] = dm.filter(regex='Cage', axis='columns').apply(lambda x: 1 if x.CageA == x.CageB else 0, axis=1)
dm.sort_values(['s1', 's2'], inplace=True)

dm0 = dm.loc[(dm.TimepointDistance == 1) & ((dm.DayA == 0) | (dm.DayB == 0)) & (dm.DatasetGroup == 1) &
             (dm.SampleTypeA == 'Feces') & (dm.SampleTypeGroup == 1)]
#dm0 = dm.loc[(dm.TimepointDistance == 1) & ((dm.DayA == 0) | (dm.DayB == 0))]

dm_timepoint3 = dm0.rename(columns={"TimepointLatestSample":"SampleID"})[["Distance", "SampleID"]]
dm_timepoint3 = dm_timepoint3.merge(meta, on='SampleID',how='left')
dm_timepoints = dm.loc[(dm.TimepointDistance == 1) & (dm.MouseGroup == 1) & (dm.SampleTypeGroup == 1) &
            (dm.TimepointLatest > 1) & (dm.SampleTypeA == "Feces")].rename(columns={'TimepointLatestSample':"SampleID"})[['Distance','SampleID']]
#dm_timepoints = dm.loc[(dm.TimepointDistance == 1) & (dm.MouseGroup == 1)  &
#            (dm.TimepointLatest > 1) ].rename(columns={'TimepointLatestSample':"SampleID"})[['Distance','SampleID']]

dm_timepoints = dm_timepoints.merge(meta, on='SampleID',how='left')
dm_timepoints = pd.concat([dm_timepoint3,dm_timepoints], ignore_index=True)
dm_timepoints.to_csv('tmp/Distance_timepoints.csv', index=False)

dm1 = dm.loc[(dm.DatasetGroup.isin([0,1])) & (dm.DayGroup == 1) & (dm.SampleTypeGroup == 1) & (dm.SampleTypeA == 'Feces')].rename(columns={'TimepointLatestSample':"SampleID"})[['Distance', 'DayGroup', 'DatasetGroup',"SampleID", "s1","s2"]]
dm1 =dm1.merge(df.drop_duplicates('SampleID')[['SampleID','Cage', "Mouse",'Day',"Dataset", "SampleType"]], on='SampleID', how='left')
dm1.loc[dm1.DatasetGroup == 0 , "Dataset"] = 'between Facilities'

dm1.loc[dm1.s1 != dm1.s2].to_csv('tmp/Distance_day_dataset_comparisons.csv', index=False)
dm.to_csv("tmp/DistanceMatrix_all.csv", index=False)

#### compare within Mouse distance to between
#
dm2 = dm.loc[(dm.MouseGroup.isin([0,1]) & (dm.DayGroup.isin([0,1])) & (dm.DatasetGroup.isin([0,1])) & (dm.SampleTypeGroup == 1) & (dm.SampleTypeA == 'Feces'))].rename(columns={'TimepointLatestSample':"SampleID"})[['Distance', 'DayGroup', 'MouseGroup',"SampleID", "s1","s2"]]
#dm2 = dm.loc[(dm.MouseGroup.isin([0,1]) & (dm.DayGroup.isin([0,1])))].rename(columns={'TimepointLatestSample':"SampleID"})[['Distance', 'DayGroup', 'MouseGroup',"SampleID", "s1","s2"]]
dm2 =dm2.merge(df.drop_duplicates('SampleID')[['SampleID','Cage', "Mouse",'Day',"Dataset", "SampleType"]], on='SampleID', how='left')
dm2.loc[dm2.s1 != dm2.s2].to_csv('tmp/Distance_day_mouse_comparisons.csv',index=False)


""" PCoAs """
b = pd.read_csv("bern_metagenomes_Evopipe_master.txt", sep='\t')
b["Dataset"] = ['bern_metagenomes' for x in b.index]
b = b[~b.PolyID.isin(mask.PolyID.unique())]
#b = b.loc[(b['Coverage_mean'] >= 100) & (b.TotalCov >= 100)]
i = pd.read_csv("bern_isolates_Evopipe_master.txt", sep='\t')
i["Dataset"] = ['bern_isolates' for x in i.index]
i = i[~i.PolyID.isin(mask.PolyID.unique())]
p = pd.read_csv('ploen_Evopipe_master.txt', sep='\t')
p['Dataset'] = ['ploen' for x in p.index]
p = p[~p.PolyID.isin(mask.PolyID.unique())]

a = pd.concat([b,i,p], ignore_index=True)

a = a[~a.PolyID.isin(mask.PolyID.unique())]

pcoa_columns = ["SampleID",'Mouse', 'Day','Cage','MetagenomeID','Dataset', "SampleType"]
make_pcoa(df.loc[(df.Day == 28) & (df.Dataset == 'ploen')], binary=True, plot=False,color="", keep_columns= pcoa_columns, NameTag="ploen", metric='jensenshannon')


print('\n\nDone!')

## --> plot.pcoa.R
