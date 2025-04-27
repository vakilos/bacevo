## download from https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/ the 2020 COG data
## fun-20.tab --> describes COG families
## cog-20.def.tab --> defines COG id to COG functional categories

"""
Comma-delimited plain text file assigning proteins to COGs
Columns:
1.	Gene ID (GenBank or ad hoc)
2.	NCBI Assembly ID
3.	Protein ID (GenBank if conforms to [A-Za-z0-9_]+\\.[0-9]+ regex; ad hoc otherwise)
4.	Protein length
5.	COG footprint coordinates on the protein. "201-400" means "from position 201 to position 400"; "1-100=201-300" indicates a segmented footprint, 1-100 AND 201-300
6.	Length of the COG footprint on the proteins
7.	COG ID
8.	reserved
9.	COG membership class (0: footprint covers most of the protein and most of the COG profile; 1: footprint covers most of the COG profile and part of the protein; 2: footprint covers most of the protein and part of the COG profile; 3: partial match on both protein and COG profile)
10.	PSI-BLAST bit score for the match between the protein and COG profile
11.	PSI-BLAST e-value for the match between the protein and COG profile
12.	COG profile length
13.	Protein footprint coordinates on the COG profile"""
import os

from bacevo_fun import *
import  pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import preprocessing

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
from rich import print
console = Console()
parser = argparse.ArgumentParser(description="Enrichment analysis")
parser.add_argument('-c', dest='cog', help='Directory with all COG annotatoion, dowloaded from https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/'
                                           'By default, current directory i selected', default=os.getcwd())
parser.add_argument('-g', dest='genbank', help='Path to the genome genbank file', required=True)
parser.add_argument('-i', dest='polytable', help='Path to <name>_FinalPolyTable.tsv', required=True)
parser.add_argument('-m', dest='mask', help='File for masking polymorphisms.', default=None)
parser.add_argument('-o', dest='out', help='Directory to export output files. By default current directory', default=os.getcwd())
args = parser.parse_args()

genbank = args.genbank
genome_length = sum([len(rec.seq) for rec in SeqIO.parse(genbank,'genbank')])


if not os.path.isfile(f'{args.cog}/cog-20.btheta.csv'):
    cog = pd.read_csv(f'{args.cog}/cog-20.cog.csv',header=None, usecols=[x for x in range(13)]) ## will associated genename to cogid
    cog = cog.loc[cog[1] == "GCF_000011065.1"] ## my genome assembly ID,
    cog[0] = cog[0].apply(lambda x:x.replace("BT_","BT"))
    cog.to_csv(f'{args.cog}/cog-20.btheta.csv', index=False)
else:
    cog = pd.read_csv(f"{args.cog}/cog-20.btheta.csv").drop_duplicates(['0','7'])

print(f"{cog.shape[0]} COG IDs for my genome")
cog_info = pd.read_csv(f'{args.cog}/cog-20.def.csv', usecols=['CogID','CogCategory', 'CogName', 'Gene_associated', 'pathway_associated'])
category_info = pd.read_csv(f"{args.cog}/fun-20.txt.csv", index_col="CogCategory")['CogCategoryDescription'].to_dict()

print(f"{cog_info.loc[(cog_info['CogCategory'].str.len() > 1) & (cog_info['CogCategory'] != 'Nan')].shape[0]} COG IDs are associated with multiple COG categories")

# create one entry per association - so that one gene might have multiple enties for multiple COG IDs.
for i, row in cog_info.iterrows():
    if len(row.CogCategory) > 1 and row.CogCategory != 'Nan':
        for x in list(row.CogCategory):
            cog_info = pd.concat([cog_info, pd.DataFrame({'CogID':[row.CogID],'CogCategory':[x], 'Gene_associated':[row['Gene_associated']], 'pathway_associated':[row['pathway_associated']]})])
        cog_info.drop(index=i,inplace=True)

if not os.path.isfile(f"{args.cog}/genbank_with_cog.csv"):
    """Associate COG info with genbank file (genes) """
    gf = gbff_to_df(genbank).drop_duplicates("locus_tag")
    cog.rename(columns={'0':'old_locus_tag', '7':'CogID'}, inplace=True)
    gf = gf.merge(cog, on='old_locus_tag', how='left')

    ### IMPORTANT: one gene might be associated with multiple COG IDs therefore COG categories
    gf["Length"] = (gf.start - gf.end).abs()
    gf = gf.merge(cog_info, on="CogID", how="left")
    gf.rename(columns={'7':'CogID'}, inplace=True)
    gf.to_csv(f'{args.cog}/genbank_with_cog.csv', index=False)
    gf = gf[['CogID','CogCategory','CogName','Gene_associated','pathway_associated','old_locus_tag','locus_tag','Length']]
else:
    gf = pd.read_csv("~/bacevo/genbank_with_cog.csv",
                 usecols=['CogID','CogCategory','CogName','Gene_associated','pathway_associated','old_locus_tag','locus_tag','Length'])

##  This code works because gf has one gene - COG association per row (so multiple rows if gene belongs to multiple COGs)
gf[['CogID', 'CogCategory']] = gf[['CogID', 'CogCategory']].fillna('nan')
gf['CogIDSize'] = gf.groupby("CogID")["Length"].transform('sum').fillna(0)
gf['CogCatSize'] = gf.groupby("CogCategory")["Length"].transform('sum').fillna(0)
gf.rename(columns={'old_locus_tag':'oldGeneName'},inplace=True)
#gf['CogCatSizeNorm'] = preprocessing.normalize([gf.CogCatSize], norm='max')[0]
#gf['CogIDSizeNorm'] = preprocessing.normalize([gf.CogIDSize], norm='max')[0]

df = pd.read_csv(args.polytable, sep='\t',
                 usecols=['Chrom','SampleID','Mouse','Day','Dataset','SampleType','TotalCov','PolyID','Ancestry',
                          'oldGeneName','GeneName','GeneLength', 'Coverage_mean', 'PolyCount', 'cds','product','pseudo']).drop_duplicates(['SampleID','PolyID'])

df = df.loc[df.SampleType == 'Feces']
keep_samples = df.loc[(df.Chrom == 'NC_004663') & (df['Coverage_mean'] >= 100), 'SampleID'].unique()
df = df.loc[(df.TotalCov.ge(100)) & (df.SampleID.isin(keep_samples))].drop_duplicates(['SampleID','PolyID'])

if args.mask:
    mask = pd.read_csv(args.mask)
    df = df[~df.PolyID.isin(mask.PolyID.unique()) & (df.Ancestry == 'evolved')]

df = df.merge(gf, on='oldGeneName', how='left')
df[['CogID', 'CogCategory','product']] = df[['CogID', 'CogCategory',"product"]].fillna('nan')


""" There are genes that are not associated with any COG category. We remove those from downstream analysis??"""
print(f"{df.drop_duplicates('GeneName').loc[df.CogID == 'nan'].shape[0]}/{df.drop_duplicates('GeneName').shape[0]} genes not associated to any COG ID or category.")


def pvalue_from_poisson(allcogsize:int, allcogmutations:int, categorysize:int, categorymutations:int):
    from scipy.stats import poisson, binom, hypergeom
    mu = (allcogmutations * categorysize) / allcogsize ### normalize genome wide lamda for cog category size
    #pvalue = 1 - poisson.cdf(categorymutations-1, mu)
    pvalue = 1 - binom.cdf(k = categorymutations-1, n = categorysize, p= allcogmutations/allcogsize)
    #pvalue = 1 - hypergeom.cdf(k = categorymutations-1, N = categorysize, n=allcogmutations, M=allcogsize)
    #pvalue = poisson.sf(categorymutations, mu)
    return pvalue


poly_cutoff = 5
gn = pd.DataFrame()
for sampleid, group in df.groupby("SampleID"):
    Ntotal = group.shape[0]
    for i, row in df.loc[df.GeneLength > 0].drop_duplicates("GeneName").iterrows():
        if not group.loc[group.GeneName == row.GeneName].empty:
            genedf = group.loc[group.GeneName == row.GeneName].drop_duplicates("PolyID")
            length = row.GeneLength
            N = genedf.shape[0]
            #print(N, length)
            score = (N / length) / (Ntotal / genome_length)
            if group.loc[group.GeneName == row.GeneName].shape[0] >= poly_cutoff:
                gn = pd.concat([gn, pd.DataFrame({'SampleID':[sampleid], 'GeneName':[row.GeneName],"oldGeneName":[row.oldGeneName], "Score":[score], "Count":[N], 'PassPolyCutoff':['yes']})])
            else:
                gn = pd.concat([gn, pd.DataFrame(
                    {'SampleID': [sampleid], 'GeneName': [row.GeneName], "oldGeneName": [row.oldGeneName],
                     "Score": [score], "Count": [N], 'PassPolyCutoff': ['no']})])
        else:
            gm = pd.concat([gn, pd.DataFrame({'SampleID':[sampleid], 'GeneName':[row.GeneName],"oldGeneName":[row.oldGeneName], "Score":[0], "Count":[0],'PassPolyCutoff':['no']})])
gn = gn.merge(df.drop_duplicates("SampleID")[['SampleID','Dataset','SampleType','Day','Mouse']], on='SampleID',how='left')



#### Enriched genes are considered only those with enrichment score > 1 and poly count more or equal to poly_cutoff
gene_p = pd.DataFrame()
for gene, group in gn.loc[gn['PassPolyCutoff'] == 'yes'].groupby(['SampleType','GeneName']):
    d = group.Score-1
    res = scipy.stats.wilcoxon(d, alternative='greater')
    if res.pvalue <= 0.05:
        gene_p = pd.concat([gene_p, pd.DataFrame({"SampleType":[gene[0]],"GeneName":[gene[1]],'oldGeneName':[group.oldGeneName.values[0]],'pvalue':[res.pvalue]})])
        print(f'{gene}: pvalue = {res.pvalue}')
print(gene_p.sort_values('pvalue'))

gn.to_csv(f'{args.out}/gene_enrichment.csv',index=False)
gene_p.to_csv(f"{args.out}/gene_pvalues.csv", index=False)




cc = pd.DataFrame()
for sampleid,group in df.groupby("SampleID"):
     Ntotal = group.shape[0] ## total polys for sample

     for cog in [x for x in df.CogCategory.unique() if x != 'nan']:
         if not group.loc[group.CogCategory == cog].empty:
             ## for each cog category count a PolyID only once (remove redundancy caused by multiple COG IDs in the same category
             cogdf = group.loc[group.CogCategory == cog].drop_duplicates("PolyID")
             Ncog = cogdf.shape[0]
             cog_length = cogdf['CogCatSize'].values[0]
             #print(cog,Ntotal, genome_length,Ncog, cog_length)
             score = (Ncog / cog_length) / (Ntotal / genome_length)
             cc = pd.concat([cc, pd.DataFrame({'SampleID':[sampleid], "CogCategory":[cog], "Score":[score]})])
         ### add zeros for no entry in the poly_table
         else:
             cc = pd.concat([cc, pd.DataFrame({'SampleID': [sampleid], "CogCategory": [cog], "Score": [0]})])

cc = cc.merge(df.drop_duplicates("SampleID")[['SampleID','Dataset','SampleType','Day','Mouse']], on='SampleID',how='left')

for cat, group in cc.loc[cc.SampleType == 'Feces'].groupby('CogCategory'):
    d = group.Score-1
    res = scipy.stats.wilcoxon(d, alternative='greater')
    if res.pvalue <= 0.05:
        print(f'{cat}: pvalue = {res.pvalue}')

for cat,group in cc.groupby(['SampleType','CogCategory']):
    d = group.Score - 1
    res = scipy.stats.wilcoxon(d, alternative='greater')
    if res.pvalue <= 0.05:
        print(f'{cat}: pvalue = {res.pvalue}')

cc['CogCategoryDesc'] = cc['CogCategory'].map(category_info)
cc['CogCategoryDesc'] = cc['CogCategoryDesc'] + " ("+cc['CogCategory']+")"
cc.to_csv(f'{args.out}/cog_cat_enrichment.csv',index=False)

ci = pd.DataFrame()
for sampleid,group in df.groupby("SampleID"):
    Ntotal = group.shape[0] ## total polys for sample

    for cog in [x for x in df.CogID.unique() if x != 'nan']:
        if not group.loc[group.CogID == cog].empty:
            ## for each cog id count a PolyID only once (there should be only one anyway)
            cogdf = group.loc[group.CogID == cog].drop_duplicates("PolyID")
            Ncog = cogdf.shape[0]
            cog_length = cogdf['CogIDSize'].values[0]
            #print(cog,Ntotal, genome_length,Ncog, cog_length)
            score = (Ncog / cog_length) / (Ntotal / genome_length)
            if cog == "COG1373":
                print(Ncog, cog_length,score)
            ci = pd.concat([ci, pd.DataFrame({'SampleID':[sampleid], "CogID":[cog], "Score":[score], "Count":[Ncog]})])
        ### add zeros for no entry in the poly_table
        else:
            ci = pd.concat([ci, pd.DataFrame({'SampleID': [sampleid], "CogID": [cog], "Score": [0], 'Count':[0]})])

ci = ci.merge(df.drop_duplicates("SampleID")[['SampleID','Dataset','SampleType','Day','Mouse']], on='SampleID',how='left')



for cat, group in ci.loc[ci.SampleType == "Feces"].groupby('CogID'):
    d = group.Score-1
    res = scipy.stats.wilcoxon(d, alternative='greater')

    if res.pvalue <= 0.05:
        print(f'{cat}: pvalue = {res.pvalue}')
        print(f'Associated genes')
        print(df.loc[df.CogID == cat]['oldGeneName'].unique())


### selected COG IDs with at least 5 observations (samples) of enrichment score > 1
select_cogIDs = []
for cat, group in ci.loc[(ci.SampleType == "Feces") & (ci.Score > 1)].groupby("CogID"):
    if group.shape[0] > 5:
        #print(group.Dataset.value_counts())
        select_cogIDs += [cat]

print(ci.loc[ci.CogID.isin(select_cogIDs)].groupby('CogID')['Score'].median().sort_values(ascending=False))
ci.to_csv(f'{args.out}/cog_id_enrichment.csv',index=False)

print('\n\nDone!')






