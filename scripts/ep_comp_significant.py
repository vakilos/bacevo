import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats

from bacevo_fun import *
parser = argparse.ArgumentParser()
parser.add_argument('-d', dest='proj_dir', help='Project directory. Defauly = CWD.', default=os.getcwd())
args = parser.parse_args()

os.chdir(args.proj_dir)
console = Console(color_system='windows')
import polars as pl

mask = pd.read_csv('tmp/confirmed_invivo.csv')['PolyID'].unique()

df = pd.read_csv("invivo_metagenomes_FinalPolyTable.tsv", sep='\t',
                 usecols=['PolyID','SampleID','Ancestry','Coverage_mean', 'Freq',
                          'AltCov','TotalCov','SampleType','oldGeneName','GeneName','Mouse','Cage','Pos',
                          'Dataset', "Day",'ChromReads','ReadCount','Chrom', 'Synonymy'])

keep_samples = df.loc[(df.Chrom == 'NC_004663') & (df['Coverage_mean'] >= 100), 'SampleID'].unique()
df = df.loc[(df.TotalCov >=100) & (~df.PolyID.isin(mask)) & (df.Dataset =='ploen') & (df.Day == 28) & (df.Ancestry == 'evolved') &
            (df.SampleID.isin(keep_samples))].drop_duplicates(['SampleID','PolyID'])

df['SampleTypeGroup'] = df.SampleType.apply(lambda x:'SI' if x in ["SIP",'SID'] else 'LI')


sg_dict = df[['SampleID','SampleTypeGroup']].set_index('SampleID')['SampleTypeGroup'].to_dict()

####
df.to_csv('tmp/invivo_metagenomes_FinalPolyTable_comp28.tsv', sep='\t', index=False)
#### Mapped = ReadCount is the total mapped reads for the library
### ChromReads is reads mapped to each chrom
df['AltSupportNorm'] = df['AltCov'] / df['ReadCount']

prepare_maaslin_files = True

if prepare_maaslin_files:

    #gf = df.loc[df.Chrom == 'NC_004663']
    #pf = df.loc[df.Chrom == 'NC_004703']
    prevalence_cutoff = 0.2
    data = pd.pivot(data=df, index='SampleID', columns='PolyID', values="AltSupportNorm").fillna(0)
    print(data.shape)
    ### filter columns (polys) on ther prevalence
    data = data[data.columns[data.astype(bool).sum(axis=0)  >= len(data) * prevalence_cutoff]]
    print(data.shape)
    cor_thres = 0.9
    data_cor = data.corr()
    upper_triangle = data_cor.where(np.triu(np.ones(data_cor.shape), k=1).astype(bool))
    # Find features with correlation higher than the threshold
    to_drop = [column for column in upper_triangle.columns if any(upper_triangle[column] > cor_thres)]
    console.print(f'Dropping {len(to_drop)} highly correlated features..')
    # Drop the identified columns
    data = data.drop(columns=to_drop)
    data.to_csv('maaslin/maaslin_data_genome.csv')

    del data
    console.print('Files for Maaslin2 have been exported!!')
    sys.exit()

df.drop_duplicates('SampleID')[['SampleID','Mouse',"SampleTypeGroup",'SampleType','Cage']].to_csv('maaslin/maasling_meta.csv', index=False)
table = pd.pivot(data=df, index= 'PolyID',columns =['SampleID','SampleTypeGroup'], values='AltSupportNorm').fillna(0)

imp1 = pd.read_csv("maaslin/genome/significant_results.tsv", sep='\t')
imp1['PolyID'] = imp1.feature.apply(lambda x: x.replace("." , ":"))


imp1.drop(columns=['feature','metadata'],inplace=True)

imp1 = imp1.rename(columns={'value':'SampleTypeGroup'})

""" Run Mann-Whitney test between LI and SI for each PolyID on frequency data. """



""" Run Mann-Whitney test between LI and SI for each PolyID on normalized AltCov. """
sign = pd.DataFrame()
from scipy.stats import mannwhitneyu, f_oneway, false_discovery_control

### remove anything with prevalence < 0.2
prevalence_cutoff = 0.2




for _,poly in zip(track(range(table.shape[0])), table.index.values):
   if table.loc[poly,:].astype(bool).sum() / table.shape[1] > prevalence_cutoff:

       stat, pvalue = mannwhitneyu(table.loc[poly, (slice(None),'SI')].values, table.loc[poly, (slice(None),'LI')].values)
       meanSI, meanLI = np.mean(table.loc[poly, (slice(None),'SI')].values), np.mean(table.loc[poly, (slice(None),'LI')].values)

       sign = pd.concat([sign, pd.DataFrame({"pvalue": [pvalue], 'PolyID': [poly], 'Pos':[df.loc[df.PolyID == poly, 'Pos'].values[0]],
                                             'GeneName':[df.loc[df.PolyID == poly, 'GeneName'].values[0]], 'meanSI':[meanSI], 'meanLI':[meanLI]})])

sign['padj'] = false_discovery_control(sign.pvalue)
print(sign)

sign = sign.merge(imp1, on='PolyID', how='outer')

sign.to_csv('tmp/compartments_stat_sign_polys_readcount.csv', index=False)


""" Random forest analysis"""
table = pd.pivot(data=df, index=['SampleID','SampleTypeGroup'], columns='PolyID', values="AltSupportNorm").fillna(0)
table[table > 0] = 1
from sklearn.preprocessing import LabelEncoder
le=LabelEncoder()
X = table.values
print(table.index)
y = table.index.get_level_values('SampleTypeGroup').values
y =le.fit_transform(y) ## replaces categories with integers - maybe use the OneHotEncoder approach instead?

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.7, random_state=42)
clf = RandomForestClassifier(n_estimators=100, min_samples_split=2)
clf.fit(X_train, y_train)
pred = clf.predict(X_test)
rand_score = clf.score(X_test, y_test)
pre_prob = clf.predict_proba(X_test)
report = classification_report(y_test, pred)
print(rand_score)
print(clf.classes_)
print(report)
print(clf.decision_path(X_train))
importances = list(clf.feature_importances_)
feature_importances = sorted([(feature,round(importance,8)) for feature,importance in zip(table.columns, importances) if importance >= 0] , key= lambda x:x[1],reverse=True)
imp_df = pd.DataFrame({'feature':[x[0] for x in feature_importances], 'importance':[x[1] for x in feature_importances]})
imp_df.to_csv('tmp/feature_importance.csv', index=False)

imp_polys = imp_df.loc[imp_df.importance > 0]

