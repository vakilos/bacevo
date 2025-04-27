
from bacevo_fun import *
from popgen_fun import *


parser = argparse.ArgumentParser(description = '')
parser.add_argument('-i', dest='dnds_file', help='Path to a gene_dNdS.csv file.')
parser.add_argument('-p', dest='polytable', help='Path to a FinalPolyTable.tsv.')
parser.add_argument('-o', dest='out', help='Directory to export files.')
parser.add_argument(-'d', dest='proj_dir', help='Project directory.', default = os.getcwd())
args = parser.parse_args()

dn = pd.read_csv(args.dnds_file)
currentDir = args.proj_dir
os.chdir(currentDir)

if not os.path.isdir(currentDir+'/tmp'): os.mkdir('tmp')

df = pd.read_csv(args.polytable, sep='\t', low_memory=False,
                 usecols=['Day','Mouse','Cage','Dataset','SampleID','PolyID','Pos','Chrom',"Ref",'Alt','TotalCov','Freq', 'SampleType',
                          'GeneName','oldGeneName','Synonymy','Ancestry','Type','MutationType',"product",
                          'cddID','cddName','cddDescription',  "CogID", "CogAnnotation","CogFamily", "DatasetPrevalencePercent", "DatasetPrevalence",
                          'MousePrevalence', "MousePrevalencePercent", 'pseudo','Coverage_mean', "Prevalence",'MeanFreq']).drop_duplicates(['SampleID','PolyID'])
mf = pd.read_csv("tmp/confirmed_invivo.csv")
mask = mf.drop_duplicates("PolyID")['PolyID'].values
df = df.loc[~df.PolyID.isin(mask)]
keep_samples = df.loc[(df.Chrom == 'NC_004663') & (df['Coverage_mean'] >= 100), 'SampleID'].unique()
df = df.loc[(df.SampleID.isin(keep_samples))].drop_duplicates(['SampleID','PolyID'])
#pd.set_option('display.max_colwidth', None)
#print(df.drop_duplicates("GeneName").loc[df.oldGeneName == "BT0436", 'cddDescription'])

"""                              Read dNdS files (after dnds_calculate.py)                       """
dn = pd.read_csv(args.dnds_file, na_values=np.nan)#.fillna(0)
dn = dn.drop_duplicates(['SampleID','GeneName'])  ## remove after fixing





dn = dn.merge(df.drop_duplicates("SampleID")[["SampleID",'Cage','SampleType']], on='SampleID', how ='left') ## merge with polytable
print(dn.shape)
select_genes = []
N_samples = 2

#### pi files all shared the same pi - piN and pS columns (this has been fixed in dndS_calculate,py) - however here need to be differentitate separately
pn = pd.read_csv('invivo_gene_pi.csv')
pn.rename(columns={'pi':'gene_pi','piN':'gene_piN', 'piS':'gene_piS'}, inplace=True)

#pc = pd.read_csv('invivo_codon_pi.csv')
ps = pd.read_csv('invivo_site_pi.csv')

#pc = pd.concat([pd.read_csv('bern_codon_pi.csv'), pd.read_csv("ploen_codon_pi.csv")])
dn = dn.merge(pn, on=['SampleID','GeneStart','GeneLength','GeneName'])
dn['total'] = dn.Gene_N_count + dn.Gene_S_count

""" Export data for compartments..."""
dn_comp = dn.loc[(dn.Dataset == 'ploen') & (dn.Day ==28)].merge(df.drop_duplicates('GeneName')[['GeneName','oldGeneName']], on='GeneName',how='left')
dn_comp.to_csv("tmp/gene_dnds_compartments.csv", index=False)

""" Select genes which have pNpS > 1 in at least 2 samples"""
for name, group in dn.loc[dn.SampleType == 'Feces'].groupby('GeneName'):
    if group.loc[group['Gene_pNpS'] > 1].shape[0] >= N_samples:
        select_genes += [name]

df = df.merge(dn, on=["SampleID",'Mouse','Dataset',"SampleType",'Cage','Day',"GeneName"], how='left')
dn = dn.merge(df.drop_duplicates('GeneName')[['GeneName','oldGeneName']], on='GeneName',how='left')
dn.loc[dn.GeneName.isin(select_genes)].to_csv("tmp/gene_dnds_fig3.csv", index=False)





#dn['Gene_pNpS'] = dn['Gene_pN'] / dn['Gene_pS']

dc = pd.read_csv("invivo_codon_dNdS.csv")#.fillna(0)


### this is per gene stats
print(dn.loc[dn.SampleType == 'Feces'].groupby(["Day"]).agg({'Gene_dNdS':'mean'}))




### pi diversity over genome - sum of pi in polymorphic sites divided by genome size
# grep "LOCUS" genomic.gbff
genome_size = 6260361.0
plasmid_size = 33038.0

chrom_size = {"NC_004663":genome_size, "NC_004703": plasmid_size}

ps = ps.groupby(["SampleID","Chrom"]).agg({'pi':'sum','piS':'sum','piN':'sum'}).reset_index()

ps = ps.merge(df.drop_duplicates("SampleID")[['Mouse',"SampleID",'Cage','Day',"Dataset", 'SampleType']],
                                                                      on='SampleID', how ='left')

ps['ChromSize'] = ps.Chrom.map(chrom_size)
ps[['pi','piN','piS']] = ps.apply(lambda x: x[['pi','piN','piS']] / x['ChromSize'], axis=1)
ps.to_csv('tmp/genome_pi.csv', index=False)

## Number of polys per gene)
gpc = df.groupby(["SampleID",'GeneName']).agg({"PolyID":'size'}).reset_index().rename(columns={'PolyID':"GenePolyCount"})
dn = dn.merge(gpc, on=['SampleID','GeneName'], how='left')
# sns.distplot(dn.loc[dn.GenePolyCount == 1]['Gene_pNpS'])
# plt.show()

df = df.merge(gpc, on=['SampleID','GeneName'], how='left')
#df['Gene_pNpS'] = df['Gene_pNpS'].fillna(-1)


def dNds_calculation_anomalies(df):
    anomalies = pd.DataFrame()

    for name,group in df.groupby(['SampleID','GeneName']):
        not_called_syn = 0
        not_called_nonsyn = 0
        if group['Synonymy'].value_counts().reset_index().shape[0] == 1:
            if group['Synonymy'].value_counts().reset_index()['Synonymy'].values[0] == 'Synonymous' and group['Gene_pS'].values[0] == 0:

                not_called_syn += 1
            if group['Synonymy'].value_counts().reset_index()['Synonymy'].values[0] == 'Nonsynonymous' and group['Gene_pN'].values[0] == 0:
                not_called_nonsyn += 1
        anomalies = pd.concat([anomalies, pd.DataFrame({'SampleID':[name[0]],'GeneName':[name[1]],
                                            'syn_anomalies':[not_called_syn], 'nonsyn_anomalies':[not_called_nonsyn]})])

    return  anomalies



df = add_trajectories_appearing(df)
df.to_csv("tmp/invivo_FinalPolyTable_figure3.tsv", sep='\t', index=False)
""" --> plot_figure3.R"""
sys.exit()








#dn.to_csv('tmp/invivo_fecal_gene_dNdS_pi.csv')
#df.to_csv('tmp/invivo_fecal_FinalPolyTable_plus_dNdSpi.csv',index=False)




pd.set_option('display.max_columns', None)
### Get persistent genes and look at higher dNdS values
from rich.console import Console
from rich_dataframe import prettify

mpr = df.loc[df.MousePrevalencePercent > 0.6].groupby(["Dataset","GeneName"]).agg({'Gene_dNdS':'median'}).sort_values(["Dataset","Gene_dNdS"],ascending=False).reset_index()

## get persistent genes for bern
persistent_genes = list(set(mpr[mpr.Dataset == 'bern_metagenomes'].head(20)["GeneName"].tolist() + mpr[mpr.Dataset == 'ploen_fecal'].head(20)["GeneName"].tolist()))
per = df.loc[df.GeneName.isin(persistent_genes)]
per.to_csv("persistent_invivo_FinalPolyTable.tsv", index=False, sep='\t')


#prettify(mpr[mpr.Dataset == 'bern_metagenomes'].head(20).merge(df.drop_duplicates("GeneName")[["oldGeneName","GeneName",'cddDescription','CogAnnotation','CogFamily']], on="GeneName",how='left'))
#prettify(mpr[mpr.Dataset == 'ploen_fecal'].head(20).merge(df.drop_duplicates("GeneName")[["oldGeneName","GeneName",'cddDescription','CogAnnotation','CogFamily']], on="GeneName",how='left'))

### Get prevalent genes and look at higher dNdS values
pre = df.loc[df.DatasetPrevalencePercent > 0.7].groupby(["Dataset","GeneName"]).agg({'Gene_dNdS':'median'}).sort_values(["Dataset","Gene_dNdS"],ascending=False).reset_index()
prevalent_genes = list(set(pre[pre.Dataset == 'bern_metagenomes'].head(20)["GeneName"].tolist() + pre[pre.Dataset == 'ploen_fecal'].head(20)["GeneName"].tolist()))
pre = df.loc[df.GeneName.isin(prevalent_genes)]
pre.to_csv("prevalent_invivo_FinalPolyTable.tsv", index=False, sep='\t')

#prettify(pre[pre.Dataset == 'bern_metagenomes'].head(20).merge(df.drop_duplicates("GeneName")[["oldGeneName","GeneName",'cddDescription','CogAnnotation','CogFamily']], on="GeneName",how='left'))
#prettify(pre[pre.Dataset == 'ploen_fecal'].head(20).merge(df.drop_duplicates("GeneName")[["oldGeneName","GeneName",'cddDescription','CogAnnotation','CogFamily']], on="GeneName",how='left'))

### hot spots
hot = df.groupby(['Dataset','GeneName']).agg({'GenePolyCount': 'median'}).sort_values(['Dataset','GenePolyCount'],ascending=False).reset_index()
hot_genes = list(set(hot[hot.Dataset == 'bern_metagenomes'].head(20)["GeneName"].tolist() + hot[hot.Dataset == 'ploen_fecal'].head(20)["GeneName"].tolist()))
hot = df.loc[df.GeneName.isin(hot_genes)]
hot.to_csv("hotspots_invivo_FinalPolyTable.tsv", index=False, sep='\t')



#### check which polys are in the isolates (datasets to check per,pre,hot)
pd.set_option('display.max_rows', None)

iso = pd.read_csv("bern_isolates_FinalPolyTable.tsv",sep='\t')




#prettify(iso.loc[iso.GeneName.isin(hot_genes)].groupby(["PolyID",'oldGeneName']).agg({"Freq":'size'}).reset_index())
print(iso.loc[(iso.GeneName.isin(hot_genes)) & (iso.GeneName != "Intergenic")].drop_duplicates(["PolyID",'oldGeneName'])[["PolyID",'oldGeneName','DatasetPrevalence',"DatasetPrevalencePercent"]]
         .sort_values("DatasetPrevalencePercent", ascending=False))
print(iso.loc[(iso.GeneName.isin(prevalent_genes)) & (iso.GeneName != "Intergenic")].drop_duplicates(["PolyID",'oldGeneName'])[["PolyID",'oldGeneName','DatasetPrevalence',"DatasetPrevalencePercent"]]
         .sort_values("DatasetPrevalencePercent", ascending=False))
print(iso.loc[(iso.GeneName.isin(persistent_genes)) & (iso.GeneName != "Intergenic")].drop_duplicates(["PolyID",'oldGeneName'])[["PolyID",'oldGeneName','DatasetPrevalence',"DatasetPrevalencePercent"]]
         .sort_values("DatasetPrevalencePercent", ascending=False))

### genes from isolate polys
print(df.loc[df.oldGeneName.isin(iso['oldGeneName'].unique())]['oldGeneName'].unique())



# prettify(iso.loc[iso.GeneName.isin(prevalent_genes)].groupby(["PolyID", "oldGeneName"]).agg({"Freq": 'size'}).reset_index())
# prettify(iso.loc[iso.GeneName.isin(persistent_genes)].groupby(["PolyID", "oldGeneName"]).agg({"Freq": 'size'}).reset_index())






"""Search for cases of dNdS increase overtime"""
dnds_time = pd.DataFrame()
for name,group in dn.groupby(['Mouse', 'GeneName']):
    if group['Gene_dNdS'].dropna().size > 1:

        R, p = scipy.stats.spearmanr(group['Day'], group["Gene_dNdS"].fillna(0))
        dnds_time = pd.concat([dnds_time,pd.DataFrame({"Mouse":[name[0]],
                                                   "GeneName":[name[1]],
                                                   "PearsonR":[R], "p-value":[p]})])
    else:
        dnds_time = pd.concat([dnds_time,pd.DataFrame({"Mouse":[name[0]],
                                                   "GeneName":[name[1]],
                                                   "PearsonR":[0], 'p-value':[1]})])
#dc = dc.merge(pc, on=['SampleID','GeneStart','GeneLength','GeneName']).merge(df.drop_duplicates("SampleID")[['Mouse',"SampleID",'Cage','Day']],
#                                                                             on='SampleID')

dnds_time = dnds_time.loc[(dnds_time.PearsonR > 0) & (dnds_time['p-value'] <= 0.05)]
### subset the dn dataframe (contains dNdS data) based on the dnds_time (that holds dnds correlation to time)
dnds_increase_overtime = pd.DataFrame()
for i, row in dnds_time.iterrows():
    dnds_increase_overtime = dnds_increase_overtime.append(dn.loc[(dn.Mouse == row.Mouse) & (dn.GeneName == row.GeneName)])
dnds_increase_overtime = dnds_increase_overtime.merge(df.drop_duplicates("GeneName")[["GeneName","oldGeneName"]], on='GeneName',how='left')
dnds_increase_overtime.fillna(0).to_csv("tmp/dnds_increase_overtime.csv", index=False)
print(dnds_increase_overtime.pivot(index=["SampleID","GeneName"], columns="Day",values="Gene_dNdS").fillna(0))

#print(dn.sort_values(['Gene_dNdS', 'pi']).head(50))
