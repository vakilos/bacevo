import argparse
import statistics
import sys

import pandas as pd
import scipy.stats

from bacevo_fun import *
from scipy.spatial.distance import squareform, pdist

parser = argparse.ArgumentParser()
parser.add_argument("-d", dest='proj_dir', help='A project directory.')
currentDir = args.proj_dir
import scipy.stats as st

os.chdir(currentDir)
if not os.path.isdir(currentDir+'/tmp'): os.mkdir('tmp')

df = pd.read_csv('invivo_metagenomes_FinalPolyTable.tsv',sep='\t') ### this is the concatenation of bern_metagenomes and ploen FinalPolyTables with the ep_merge_tables.py script
df = df.loc[df.SampleType == 'Feces']

print(df.drop_duplicates("SampleID")['Coverage_mean'].median())
print(df.loc[df.Day == 3].drop_duplicates("SampleID")['PolyCount'].mean())
print(df.loc[df.Day == 3].drop_duplicates("SampleID")['PolyCount'].median())
print(st.t.interval(0.95, len(df.loc[df.Day == 3].drop_duplicates("SampleID")['PolyCount'])-1, loc=np.mean(df.loc[df.Day == 3].drop_duplicates("SampleID")['PolyCount']), scale=st.sem(df.loc[df.Day == 3].drop_duplicates("SampleID")['PolyCount'])))
print(df.loc[df.Day >= 7].drop_duplicates("SampleID")['PolyCount'].mean())
print(df.loc[df.Day >= 7].drop_duplicates("SampleID")['PolyCount'].median())
for name,group in df.drop_duplicates("SampleID").groupby("Day"):
    print(name)
    print(group['PolyCount'].agg(['mean','median', 'max']))
    print(st.t.interval(0.95, len(group.PolyCount)-1, loc=np.mean(group.PolyCount), scale=st.sem(group.PolyCount)))
for name,group in df.drop_duplicates("SampleID").groupby("Mouse"):
    print(name)
    print(group['PolyCount'].agg(['mean','median']))


mask = pd.read_csv("tmp/confirmed_invivo.csv") ### this is polymorphisms that have been masked from shuffling regions..
df = df[~df.PolyID.isin(mask.PolyID.unique())]

keep_samples = df.loc[(df.Chrom == 'NC_004663') & (df['Coverage_mean'] >= 100), 'SampleID'].unique()

df = df.loc[(df.TotalCov.ge(100)) & (df.SampleID.isin(keep_samples))].drop_duplicates(['SampleID','PolyID'])


console = Console(color_system='windows')

#### Find late prevalent polys
late = find_late_prevalent_polys(df)
early = find_early_persistent_polys(df)




para = find_parallel_polys(df)
late.to_csv("tmp/LatePrevalentPolys.csv", index=False)
para.to_csv("tmp/ParallelPolys.csv", index=False)
indi = find_individual_peristent_polys(df)
early.to_csv('tmp/EarlyPersistentPolys.csv', index=False)
indi.to_csv('tmp/IndividualPersistentPolys.csv',index=False)


"""                   FIGURE 1                                 """

###########################################            Distance matrix for PERMANOVA (Fig 1A)
matrix = pd.pivot(df.loc[(df.SampleType == "Feces") & (df.SampleProfile == 'metagenome')], columns="PolyID", index='SampleID', values="Freq").fillna(0)
matrix[matrix > 0] = 1 ## binarize matrix
dm = pd.DataFrame(data=squareform(pdist(matrix, metric='euclidean')), index=matrix.index, columns=matrix.index).reset_index().\
    merge(df.drop_duplicates('SampleID')[['SampleID','Day',"Mouse",'Cage','Dataset']], on='SampleID',
          how='left')
dm.loc[:, ~dm.columns.isin(['Day','Cage','Mouse','Dataset'])].to_csv('tmp/DistanceMatrix.csv', index=False)
dm[["SampleID",'Mouse','Day','Cage','Dataset']].to_csv("tmp/DM_metadata.csv", index=False)

console.print(f'Exported tmp/DistanceMatrix.csv and tmp/DM_metadata.csv')

############################################            per library non-syn, syn and other SNP counts (Fig 1B)
synonymy_count = pd.DataFrame()
for name, group in df.loc[(df.SampleType == "Feces") & (df.SampleProfile == 'metagenome') & (df.Type == 'SNP')].groupby(
    ["SampleID", "Day", 'Mouse', 'Dataset', 'Synonymy']):
    synonymy_count = pd.concat([synonymy_count,
                                pd.DataFrame({"SampleID":[group.SampleID.values[0]],
                                              "Day":[group.Day.values[0]],
                                              "Mouse":[group.Mouse.values[0]],
                                              "Dataset":[group.Dataset.values[0]],
                                              "Synonymy":[group.Synonymy.values[0]],
                                              "Count":[group.shape[0]]})])



synonymy_count.to_csv('tmp/SynonymyCount.csv', index=False)
console.print(f'Exported tmp/SynonymyCount.csv')


############################################   per library MutationType (Transversions and Transitions) count (Fig 1C)
df['MutationType2'] = ["" for x in df.index]
df.loc[((df.Ref == "A") & (df.Alt == "G")) | ((df.Ref == "G") & (df.Alt == 'A')), "MutationType2"] = 'A \u2194 G'
df.loc[((df.Ref == "C") & (df.Alt == "T")) | ((df.Ref == 'T') & (df.Alt == 'C')), "MutationType2"] = 'C \u2194 T'
df.loc[((df.Ref == "A") & (df.Alt == 'C')) | ((df.Ref == "C") & (df.Alt == 'A')), "MutationType2"] = 'A \u2194 C'
df.loc[((df.Ref == "A") & (df.Alt == 'T')) | ((df.Ref == "T") & (df.Alt == 'A')), "MutationType2"] = 'A \u2194 T'
df.loc[((df.Ref == "C") & (df.Alt == 'G')) | ((df.Ref == "G") & (df.Alt == 'C')), "MutationType2"] = 'C \u2194 G'
df.loc[((df.Ref == "G") & (df.Alt == 'T')) | ((df.Ref == "T") & (df.Alt == 'G')), "MutationType2"] = 'G \u2194 T'

mutationtype_count = df.loc[(df.SampleType == "Feces") & (df.SampleProfile == 'metagenome') & (df.Type == "SNP")].groupby(["SampleID","Day","Cage","Mouse","Dataset","MutationType2",'Ancestry'])\
    .agg({"PolyID":'size'}).reset_index().rename(columns={"PolyID":'Count'})

print(df.loc[(df.SampleType == "Feces") & (df.SampleProfile == 'metagenome') & (df.Type == "SNP") & (df.Dataset == 'ploen') &(df.Ancestry == 'ancestral'), 'MutationType2'].value_counts(normalize=True))
print(df.loc[(df.SampleType == "Feces") & (df.SampleProfile == 'metagenome') & (df.Type == "SNP") & (df.Dataset == 'ploen') &(df.Ancestry == 'ancestral')].drop_duplicates("PolyID").shape[0])
print(df.loc[(df.SampleType == "Feces") & (df.SampleProfile == 'metagenome') & (df.Type == "SNP") & (df.Dataset == 'bern_metagenomes') &(df.Ancestry == 'ancestral'), 'MutationType2'].value_counts(normalize=True))
print(df.loc[(df.SampleType == "Feces") & (df.SampleProfile == 'metagenome') & (df.Type == "SNP") & (df.Dataset == 'bern_metagenomes') &(df.Ancestry == 'ancestral')].drop_duplicates("PolyID").shape[0])
mutationtype_count.to_csv('tmp/MutationTypesCount.csv', index=False)
console.print(f'Exported tmp/MutationTypesCount.csv')



###################### Enrichment_analysis

os.system("python ../scripts/ep_enrichment.py")
## Will export tmp/Cog20-CategoryCounts.csv
console.print(f'Exported tmp/Cog20-CategoryCounts.csv')
console.print(f'\n\nUse the exported files with plot_figure1.R')


print('THIS IS THE END!!')
