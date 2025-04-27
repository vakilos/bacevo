
from bacevo_fun import *
from popgen_fun import *

parser = argparse.ArgumentParser(description = '')
parser.add_argument('-i', dest='dnds_file', help='Path to a gene_dNdS.csv file.')
parser.add_argument('-p', dest='polytable', help='Path to a FinalPolyTable.tsv.')
parser.add_argument('-o', dest='out', help='Directory to export files.')
args = parser.parse_args()
dn = pd.read_csv(args.dnds_file)

print(dn.shape)
dn = dn.loc[~dn['Gene_pNpS'].isna()]
df = pd.read_csv(args.polytable,sep='\t', dtype={'RefCov':int,'AltCov':int})
meta = df.drop_duplicates("SampleID")[["SampleID",'Replicate','Mix']]
print(f'{df.GeneName.unique().size} different mutated gens in poly table.')
print(f'{df.loc[(df.Type == "SNP") & (df.translation != "Nan") & (df.cds.str.len() % 3 == 0) & (df.Ancestry == "evolved")].GeneName.unique().size} different mutated protein coding genes (with SNPs) in poly table.')

miss_counter = 0
detected = 0
out_df = pd.DataFrame()
df_list = []

for name, group in df.loc[(df.Type == "SNP") & (df.translation != 'Nan') & (df.cds.str.len() % 3 == 0) & (df.Ancestry == 'evolved')].drop_duplicates('PolyID').sort_values("Pos").groupby(['GeneName','SampleID']):
    if not group.loc[group.Synonymy == 'Synonymous'].empty and not group.loc[group.Synonymy == 'Nonsynonymous'].empty:
        df_list.append("|".join([str(x) for x in name]))
print(f'{len(df_list)} cases (Sample - Gene) with at least on syn and on nonsyn mutation')
for x in df_list:
    genename, sampleID = x.split("|")
    if dn.loc[(dn.SampleID == int(sampleID)) & (dn.GeneName == genename)].empty:
        miss_counter += 1
    else:
        detected += 1
print(miss_counter, detected)



gpc = df.groupby(["SampleID",'GeneName'], as_index=False).agg({'Freq':'size'}).rename(columns={"Freq":"GenePolyCount"})
df = df.merge(gpc, how='left', on=['SampleID',"GeneName"])


dnds_cutoff = 1

target_genes_count = dn.loc[dn["Gene_pNpS"] >= dnds_cutoff]['GeneName'].value_counts().reset_index()
target_genes = target_genes_count.loc[target_genes_count["count"] >= 1]['GeneName']

print(dn.loc[dn['Gene_pNpS'] > 1]['GeneName'].unique().size)
dnds = dn.loc[(dn['GeneName'].isin(target_genes))].pivot(index='SampleID',columns='GeneName',values="Gene_pNpS").stack().reset_index().rename(columns={0:"Gene_pNpS"})

dnds = dnds.merge(meta,on='SampleID', how='left')

dnds = dnds.merge(df.drop_duplicates('GeneName')[["GeneName",'oldGeneName', 'GenePolyCount']],on="GeneName", how='left')
dnds.to_csv(f'{args.out}/genes_with_high_pNpS.csv', index=False)

dn = dn.merge(df.drop_duplicates('GeneName')[["GeneName",'GenePolyCount']],on="GeneName", how='left')
dn.to_csv(f'{args.out}/dnds_data_to_R.csv', index=False)

print('\n\nDone!')