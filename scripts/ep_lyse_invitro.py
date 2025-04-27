#!/apps/python3/3.9.0/bin/python3
import sys

from bacevo_fun import *
parser = argparse.ArgumentParser(description="Creates a master table from breseq files combined with a meta data file. Run this script after ep_parse.py.")
parser.add_argument('-i', dest='inFile', help='Path to AnnotatedPolyTable.txt (From ep_parse.py)', required=True)
parser.add_argument("-n", dest='name', help="Add a tag to the exported master table.", default=None)
parser.add_argument('-a', dest='alt', help="Value for the AltCov filter. Equal or greater values are selected. DEFAULT=0", default=0)
parser.add_argument('-f', dest='freq', help="Value for the Freq filter. Equal or greater values are selected. DEFAULT=0", default=0)
parser.add_argument('-s', dest='sampletype', help="Use the parameter only if 'SampleType' column is included in the metadata.", default=False, action='store_true')
parser.add_argument('-c', dest='eval_coverage', help="Evaluate the number of polys for a range of genome coverage (optional). DEFAULT=False",default=False, action='store_true')

args = parser.parse_args()
console = Console()

currentDir = os.getcwd()
console.print(f'Working directory set to {currentDir}\n')

##########  set global variables
mode = "vitro"
ancestralFreqThres = 0.98
total_support = 0  ### count statistics --> std error of count is square root
alt_support = int(args.alt)
freq_filter = float(args.freq)
lda_filter = False ## for no filter use False, for filtering use a positive number :)

remove_ancestral = True ## will remove ancestral polys
include_only_denovo = False ## will include only new polymorphic positions not present in the ancestral (Reversed = 0 )
run_filter  = 'VBCF' ## choose between all, VBCF and BSF
run_coverage_threshold = False ## run this only the first time

currentDir = os.path.dirname(args.inFile)
NameTag = args.name+"_" if args.name else ""
BioUnit = 'Replicate'
Condition = 'Mix'

exported_file_name =NameTag+"FinalPolyTable"

os.chdir(currentDir)
if not os.path.isdir("Figures"): os.mkdir("Figures")
if not os.path.isdir("tmp"): os.mkdir("tmp")

myFile = args.inFile
df = pd.read_csv(myFile, sep='\t').rename(columns={'locus_tag':"GeneName",'old_locus_tag':'oldGeneName'})

print(df[df.duplicated(keep=False)]['Type'].value_counts())
assert df.shape[0] == df.drop_duplicates(['PolyID',"SampleID"]).shape[0]
df[['AltCov', 'TotalCov',"RefCov", "Freq"]] = df[['AltCov', 'TotalCov', "RefCov", 'Freq']].astype(float)
df['easy_name'] = "m"+df['Mix'].astype(str)+"_r"+df["Replicate"].astype(str) ## this was specific to invitro added again if needed

df = df.assign(SampleProfile=[None for x in df.index]) ### isolate or metagenome
df.loc[(df["SampleID"] != df["MetagenomeID"]), "SampleProfile"] = 'isolate'; df.loc[(df["SampleID"] == df["MetagenomeID"]), "SampleProfile"] = 'metagenome'


N_isolates = df.drop_duplicates('SampleID').loc[df.SampleProfile == 'isolate'].shape[0]
N_metagenomes = df.drop_duplicates('SampleID').loc[df.SampleProfile == 'metagenome'].shape[0]
initial_N_samples = df.drop_duplicates('SampleID').shape[0]
initial_enties = df.shape[0]
initial_different_polys = len(df['PolyID'].unique())
table1 = Table(title='Parameters and filters')
for x in ['Parameter name', 'Value']:
    table1.add_column(x)
parameters = {x:y for x,y in zip(["Dataset name","BioUnit", "Condition","AncestralFreqThres", "TotalCov filter (>=)",'AltCov filter (>=)','Freq filter (>=)', 'LDA filter', 'Remove ancestral filter', "Sequencing Run filter"],
                                 [NameTag.rstrip("_"), BioUnit, Condition,ancestralFreqThres, total_support,alt_support, freq_filter,lda_filter, remove_ancestral, run_filter])}
for x in parameters:
    if 'filter' in x:
        table1.add_row(f'[#f57089 b]{x}', f'[b]{parameters[x]}')
    else:
        table1.add_row(f'[blue b]{x}', f'[b]{parameters[x]}')
console.print(table1)


if run_coverage_threshold:
    bench = benchmark_coverage_threshold(df)
    bench.to_csv("tmp/"+NameTag+"BenchmarkCoverageThreshold.csv",index=None)


console.print(f'\n[b]Collecting data and estimating parameters...')


"""                            Collect coverage information and estimate other parameters in the same dataframe                         """


df["MaxFreqTraj"] = df.groupby([BioUnit, "PolyID"])["Freq"].transform("max") ### Max frequency per replicate/mouse
df["BioUnitPolyAbundance"] = df.groupby([BioUnit,"PolyID"])["PolyID"].transform(len) ## How many times the poly appears for a replicate/mouse's samples..
df["Prevalence"] = df.groupby("PolyID")["Day"].transform(len) ### How many times each poly is detected in the dataset....

add_group_poly_prevalence_and_mean(df, column="Mix")

def add_MixReplicateNumber_col(df):
    mix_replicates = df.drop_duplicates('SampleID').groupby("Mix").size().to_dict()
    df["MixReplicateNumber"] = df.Mix.map(mix_replicates)

add_MixReplicateNumber_col(df)
df['MixPrevalencePercent'] = df.MixPrevalence / df.MixReplicateNumber
df["CumFreq"] = df.groupby([BioUnit,"PolyID"])["Freq"].transform("sum") ### useless????
df["MeanFreq"] = df.groupby([BioUnit, "PolyID"])["Freq"].transform("mean") ### Mean poly frequency in replicate/mouse's samples
df['GeneName'] = df["GeneName"].replace({"Nan":"Intergenic"})
df['oldGeneName'] = df["oldGeneName"].replace({"Nan": "Intergenic"})
df['GeneLength'] = df['cds'].apply(lambda x:len(x) if isinstance(x, str) else 0)
df["CogFamilySize"] = df.groupby('CogFamily')['GeneLength'].transform("sum")
df["CogSize"] = df.groupby('CogID')['GeneLength'].transform("sum")
df['CogID'] = df["CogID"].replace({"Nan": "Intergenic"})
df['Run'] = df.SampleID.astype(str).apply(lambda x:'VBCF' if x.startswith("6") or x.startswith("512") or x.startswith("3") and x[-1] != "B" else "" "Miseq" if x.startswith("B") else "BSF" ) ### mark samples based on the sequencing workflow
vbcf = df.loc[df['Run'] == "VBCF", "PolyID"].unique()
bsf = df.loc[df["Run"] == "BSF", "PolyID"].unique()
onlyVBCF = [snp for snp in vbcf if snp not in bsf] ### select polys present exclusively in VBCF samples
onlyBSF = [snp for snp in bsf if snp not in vbcf] ### select polys present exclusively in BSF samples
shared = list(set(vbcf).intersection(bsf)) ## select polys shared between sequencing workflows

def seqPatternFunc(x):
    if x in onlyVBCF:
        return "VBCF"
    if x in onlyBSF:
        return "BSF"
    if x in shared:
        return "shared"

df['SequencingPattern'] = df['PolyID'].apply(lambda x:seqPatternFunc(x))
add_group_poly_prevalence_and_mean(df,column='Run')


console.print(f"Writing file with NO filters [u i]{exported_file_name}_nofilter.tsv[/u i].")
df.to_csv(exported_file_name + "_nofilter.tsv", index=None, na_rep="Nan", sep='\t')


console.print(f"\n[b]Filtering data...")


if lda_filter is not False:
    df = df.loc[df.LDAcoef.between(-lda_filter, lda_filter)]
    print("Number of Polymorphisms in the dataset after LDA filtering (-{} <=coef <={}): {}".format(lda_filter,lda_filter, df.drop_duplicates("PolyID").shape[0]))

""" Remove ancestral SNPs..."""
if remove_ancestral is True:
    df = df.loc[df.Ancestry == 'evolved']

"""Hard Filter Polys on Reads supporting a site  """
print(color.BOLD+color.RED+"\nRemove calls that are supported by less than {} reads in total and {} for alterantive state or below {} frequency".format(total_support,alt_support, freq_filter)+color.END)
df = df.loc[df["TotalCov"].astype(int) >= total_support]


"""Hard Filter Polys on Reads supporting the Alternative state  """
df = df.loc[df["AltCov"] >= alt_support]

""" Hard Filter Polys on frequency  """
df = df.loc[df['Freq'] >= freq_filter]


""" Keep only libraries from a Run if needed """
if run_filter != 'all':
    df = df.loc[df.Run == run_filter]
"""                            Collect coverage information and estimate other parameters in the same dataframe                         """
added_columns = {x:y for x,y in zip(["MetagenomeID","SampleProfile",'MaxFreqTraj', f'{BioUnit}PolyAbundance', 'Prevalence', "DatasetPrevalence","DatasetPrevalencePercent",f'MixPrevalence',f'MixPrevalencePercent', 'CumFreq',
                                     'MeanFreq', "DatasetMeanFreq",f'MixMeanFreq',
                                     'PolyCount' ,
                                     "GeneName", "oldGeneName", "GeneLength",  'cds',"cddID", 'cddName','cddDescription', "CogID",  "CogFamily", "CogFamilySize", "CogSize",
                                     "Run", "SequencingPattern", "RunPrevalence"],
                                    ["If isolate the SampleID of the metagenome it derived from else same as SampleID","isolate or metagenome",f"Max Frequency of poly in the {BioUnit}",f'How many times the poly appears in data',
                                     f'Poly prevalence in the dataset', "Poly's prevalence within Dataset", "Poly's prevalence within Dataset (% of samples)",
                                     f"Poly's prevalence within Mix",f"Poly's prevalence within Mix (% of samples)",
                                     "Poly's cumulative frequency over time",
                                     "Poly's mean frequency in data ","Poly's mean frequency within dataset",f"Poly's mean frequency with {BioUnit}",
                                    "Count of polys for each sample.", "Name of the gene locus or intergenic",'Old name for the gene or intergenic',
                                     'Length of the CDS','The actual CDS',"cdd ID",'cdd Name','cdd Description', "COG ID",'COG Family',
                                     'Size of the COG Family (bp)', "Size of COG (bp)", "Sequencing Run",
                                     'Describes if poly is unique to one or both sequencing runs.', "Prevalence of poly within each sequencing run",
                                     ""])}
table2 = Table(title="Added columns in the FinalPolyTable")
for x in ['Column name','Description']:
    table2.add_column(x)
for x in added_columns:
    if 'Freq' in x:
        table2.add_row(f'[b #f1c50d]{x}', added_columns[x])
    elif 'Prevalence' in x:
        table2.add_row(f'[b #d16666]{x}', added_columns[x])
    elif 'Gene' in x:
        table2.add_row(f'[b #3ba3eb]{x}', added_columns[x])
    elif 'cdd' in x:
        table2.add_row(f'[b #3581bd]{x}', added_columns[x])
    elif 'Cog' in x:
        table2.add_row(f'[b #3aaca3]{x}', added_columns[x])
    elif x in ['MetagenomeID','SampleProfile']:
        table2.add_row(f'[b #e68331]{x}', added_columns[x])
    else:
        table2.add_row(f'[b ]{x}', added_columns[x])

console.print(table2)

df["Prevalence"] = df.groupby("PolyID")["Day"].transform(len) ### How many times each poly is detected in the dataset....
add_group_poly_prevalence_and_mean(df, column="Mix")
df["PolyCount"] = df.groupby(["SampleID","Chrom"])["PolyID"].transform(len) ### How many polys a sample has....


print("\nPolys in the dataset after hard filtering "+color.BOLD+color.CYAN+": {}\n".format(len(df["PolyID"].unique()))+color.END )

console.print("Dumping the dataframe as 'masterTable.obj'.")
console.print(f"Exporting [i u]{exported_file_name}.tsv[/i u]...\n")
pickle.dump(df, open("masterTable.obj",'wb'))
df.to_csv(exported_file_name+'.tsv', index=None, na_rep="Nan", sep='\t')
final_N_samples = df.drop_duplicates("SampleID").shape[0]
final_N_isolates = df.drop_duplicates('SampleID').loc[df.SampleProfile == 'isolate'].shape[0]
final_N_metagenomes = df.drop_duplicates('SampleID').loc[df.SampleProfile == 'metagenome'].shape[0]
table3 = Table(title="Poly numbers")
for x in ['Stage','Over Poly Number','Poly Catalog', 'Samples', "Metagenomes", 'Isolates']:
    table3.add_column(x)
table3.add_row("Start", f'{initial_enties}', f'{initial_different_polys}', f'{initial_N_samples}', f'{N_metagenomes}',f'{N_isolates}')
table3.add_row("End", f'{df.shape[0]}', f'{len(df["PolyID"].unique())}', f'{final_N_samples}', f'{final_N_metagenomes}',f'{final_N_isolates}')

console.print(table3)



pf = pd.pivot(df, index="PolyID", columns='SampleID', values="Freq").fillna(0).reset_index()
console.print(f"Exporting [i u]{currentDir}/tmp/PolyFrequencyHeatmap.tsv[/i u]...\n")
console.print(f"Exporting [i u]{currentDir}/tmp/PolyFrequencyHeatmapMetaColumns.tsv[/i u]...\n")

pf.to_csv(currentDir+"/tmp/PolyFrequencyHeatmap.tsv",sep='\t',index=False)
pf.merge(df.drop_duplicates('PolyID'), on='PolyID', how='left').to_csv(currentDir+'/tmp/PolyFrequencyHeatmapMetaRows.tsv',sep='\t',index=False)
df.drop_duplicates('SampleID').to_csv(currentDir+'/tmp/PolyFrequencyHeatmapMetaColumns.tsv',sep='\t',index=False)

print("\n\nDone!")



