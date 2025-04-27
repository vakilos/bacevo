import re, sys, os, argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance
import seaborn as sns
import pickle
from bacevo_fun import *

parser = argparse.ArgumentParser()
parser.add_argument('-d', dest='proj_dir', help='Project directory. Default= CWD', default=os.getcwd())
args = parser.parse_args()

currentDir = args.proj_dir
os.chdir(currentDir)
if os.path.isdir(currentDir+'/tmp') is False:
    os.mkdir('tmp')




df = pd.read_csv('invitro_FinalPolyTable.tsv',sep='\t', dtype={'RefCov':int,'AltCov':int})
t = pd.read_csv("~/bacevo/invitro/gene_order/target_genes.csv")

t = t.loc[~t.GeneName.isin(["BT_RS02150", "BT_RS02195",'BT_RS02220',
                            "BT_RS11455", "BT_RS16445", "BT_RS17590",
                            "BT_RS17635", "BT_RS17695", "BT_RS17700", "BT_RS22910", 'BT_RS07695'])]
invitro_mask_genes = t.loc[t.GeneGroup.isin(['BT1042',"BT2268", 'BT0436',"BT3477", 'BT4540', 'BT3239'])]

the_mask = df.loc[df.GeneName.isin(invitro_mask_genes.GeneName.values)]
the_mask.to_csv("~/bacevo/invitro/tmp/confirmed_invitro.csv",index=False)
sample_order = df.drop_duplicates("SampleID").sort_values(['Mix'])['SampleID'].tolist()
add_cols = ["PolyID","Chrom", "Pos", "Ref", "Alt", 'oldGeneName','GeneName', 'Synonymy','Ancestry','Type','MutationType',"product", 'cddID','cddName','cddDescription',  "CogID", "CogAnnotation","CogFamily" ]
piv = pd.pivot(data=the_mask, index="PolyID",columns='SampleID',values='Freq').fillna(0)
piv = piv.reset_index().merge(the_mask.drop_duplicates(['PolyID',"Chrom", "Pos", "Ref", "Alt"])[add_cols], on='PolyID',how='left')
piv.to_csv("invitro_the_mask_frequencies.tsv",sep='\t',index=False)

df = df.loc[~df.GeneName.isin(invitro_mask_genes.GeneName.values)]
#df.to_csv("~/bacevo/invitro/tmp/confirmed_invitro.csv",index=False)

df.to_csv('invitro_masked_FinalPolyTable.tsv', sep='\t', index=False)


### portion of ancestal polys

print(df.drop_duplicates("PolyID").Ancestry.value_counts(normalize=True))
print(df.Ancestry.value_counts(normalize=True))



#### Created Frequency-table
sample_order = df.drop_duplicates("SampleID").sort_values(['Mix'])['SampleID'].tolist()
add_cols = ["PolyID","Chrom", "Pos", "Ref", "Alt", 'oldGeneName','GeneName', 'Synonymy','Ancestry','Type','MutationType',"product", 'cddID','cddName','cddDescription',  "CogID", "CogAnnotation","CogFamily" ]
piv = pd.pivot(data=df, index="PolyID",columns='SampleID',values='Freq').fillna(0)
piv = piv.reset_index().merge(df.drop_duplicates(['PolyID',"Chrom", "Pos", "Ref", "Alt"])[add_cols], on='PolyID',how='left')
piv.to_csv("invitro_FinalPolyTable_frequencies.tsv",sep='\t',index=False)

meta = df.drop_duplicates("SampleID")[["SampleID",'Replicate','Mix']]



"""                                Unique variants per Mix                    """
mdf = pd.DataFrame() ## store unique variants per mix
df["MixUnique"] = [0 for x in df.index]
for name, group in df.groupby("Mix"):
    mix_SNPs = group["PolyID"].values
    unique_mix_SNPs = [x for x in mix_SNPs if x not in df.loc[df.Mix != name, "PolyID"].values] ## unique SNPs to the mix - not present in any other mix

    #mdf = mdf.append(df.loc[(df['PolyID'].isin(unique_mix_SNPs)) & (df["Mix"] == name)])
    df.loc[df["PolyID"].isin(unique_mix_SNPs), "MixUnique"] = name

for n, g in df.groupby("PolyID"):
    if len(g.Mix.unique()) == 3:### present in all three mixes
        df.loc[df.PolyID == n, "MixUnique"] = 'all'
df.drop_duplicates(["PolyID",'Mix']).to_csv("tmp/MixUnique.tsv",sep='\t',index=False)

print('\n\nDone!')








