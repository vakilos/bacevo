#!/apps/python3/3.9.0/bin/python3
import sys

from bacevo_fun import *

parser = argparse.ArgumentParser(
    description="Creates the FinalPolyTable.tsv file from an AnnotatedPolyTable.txt file. Run this script after ep_parse.py.")
parser.add_argument('-i', dest='inFile', help='Path to an AnnotatedPolyTable.txt (From ep_parse.py)', required=True)
# parser.add_argument('-r', dest='ref', help='directory to reference file. (fasta, genbank, gff3)', required=True)
parser.add_argument("-n", dest='name', help="Add a name tag to the exported master table. i.e. <name>_FinalPolyTable.tsv, DEFAULT=''.", default="")
parser.add_argument('-a', dest='alt', help="Value for the AltCov filter. Equal or greater values are selected. DEFAULT=0", default=0)
parser.add_argument('-f', dest='freq', help="Value for the Freq filter. Equal or greater values are selected. DEFAULT=0", default=0)
parser.add_argument('-t', dest='total', help='Total read support for polymorphic site. DEFAUL=0, -t 100 is recommended for strongly supported sites for deeply sequenced metagenomes.', default=0, type=int)
parser.add_argument('-s', dest='sampletype', help="Use the parameter only if 'SampleType' column is included in the metadata.", default=False, action='store_true')
parser.add_argument('-c', dest='eval_coverage', help="Evaluate the number of polys for a range of genome coverage (optional). DEFAULT=False",default=False, action='store_true')
parser.add_argument('-p', dest='fecal', help='Subset only for fecal metagenomes.', action='store_true')
args = parser.parse_args()
currentDir = os.path.dirname(args.inFile)
console = Console()
console.print(f'Working directory set to {currentDir}\n')
##########  set global variables
#ancestralFreqThres = 0.9 ### fixed muatation in the ancestral population thus it is removed from the analysis is defined in ep_parse.py ( Freq >= 0.9 in the ancestral population)
total_support = args.total  ### count statistics --> std error of count is square root
alt_support = int(args.alt)  ## equal or greater than filter
freq_filter = float(args.freq) ## equal or greater than filter
lda_filter = False  ## for no filter use False, for filtering use a positive number :)
remove_ancestral = True  ## will remove ancestral polys
include_only_denovo = False  ## will include only new polymorphic positions not present in the ancestral (Reversed = 0 )
run_filter = 'all'  ## choose between all, VBCF and BSF
run_coverage_threshold = args.eval_coverage ## run this only the first time
NameTag = args.name
NameTag = args.name+"_" if args.name else ""
subset_fecal = args.fecal
BioUnit = 'Mouse'
Condition = 'Cage'
with_SampleType = args.sampletype
exported_file_name = NameTag + "FinalPolyTable"
myFile = args.inFile

os.chdir(currentDir)
if not os.path.isdir("Figures"): os.mkdir("Figures"), Console().print('fCreating [u i]Figures[/u i] directory.')
if not os.path.isdir("tmp"): os.mkdir("tmp"), Console().print('fCreating [u i]tmp[/u i] directory.')

df = pd.read_csv(myFile, sep='\t', dtype={'AltCov':int,'TotalCov':int,'RefCov':int,'Freq':float}).rename(columns={'locus_tag':"GeneName",'old_locus_tag':'oldGeneName'})
assert df.shape[0] == df.drop_duplicates(['PolyID',"SampleID"]).shape[0]
#df[['AltCov', 'TotalCov',"RefCov", "Freq"]] = df[['AltCov', 'TotalCov', "RefCov", 'Freq']].astype(float)



df = df.assign(SampleProfile=[None for x in df.index]) ### isolate or metagenome
df.loc[(df["SampleID"] != df["MetagenomeID"]), "SampleProfile"] = 'isolate'; df.loc[(df["SampleID"] == df["MetagenomeID"]), "SampleProfile"] = 'metagenome'
df.loc[df.SampleProfile == 'isolate', 'SampleType'] = 'Feces (isolate)' ### quick fix

df['Dataset'] = [NameTag.strip("_") for x in df.index]
df['Mouse'] = df['Dataset'].apply(lambda x: x[0]) + df['Mouse'].astype(str)

if subset_fecal:
    df = df.loc[df.SampleType == 'Feces']


N_isolates = df.drop_duplicates('SampleID').loc[df.SampleProfile == 'isolate'].shape[0]
N_metagenomes = df.drop_duplicates('SampleID').loc[df.SampleProfile == 'metagenome'].shape[0]
initial_N_samples = df.drop_duplicates('SampleID').shape[0]
initial_enties = df.shape[0]
initial_different_polys = len(df['PolyID'].unique())
table1 = Table(title='Parameters and filters')
for x in ['Parameter name', 'Value']:
    table1.add_column(x)
parameters = {x:y for x,y in zip(["Dataset name","BioUnit", "Condition", "TotalCov filter (>=)",'AltCov filter (>=)','Freq filter (>=)', 'LDA filter', 'Remove ancestral filter', "Sequencing Run filter"],
                                 [NameTag.rstrip("_"), BioUnit, Condition, total_support,alt_support, freq_filter,lda_filter, remove_ancestral, run_filter])}
for x in parameters:
    if 'filter' in x:
        table1.add_row(f'[#f57089 b]{x}', f'[b]{parameters[x]}')
    else:
        table1.add_row(f'[blue b]{x}', f'[b]{parameters[x]}')
console.print(table1)

if run_coverage_threshold: benchmark_coverage_threshold(df).to_csv("tmp/" + NameTag + "BenchmarkCoverageThreshold.csv", index=None)


console.print(f'\n[b]Collecting data and estimating parameters...')
"""                            Collect coverage information and estimate other parameters in the same dataframe                         """
#print("Adding columns: 'MaxFreqTraj', 'BioUnitPolyAbundance', 'Prevalence', 'MixPrevalence', 'CumFreq', 'PolyCount', 'MeanFreq',"
#    " GeneName, GeneLength, CogFamilySize, CogSize, CogID, Run, SequencingPattern, RunPrevalence ")

df["MaxFreqTraj"] = df.groupby([BioUnit, "PolyID"])["Freq"].transform('max')  ### Max frequency per replicate/mouse
df["BioUnitPolyAbundance"] = df.groupby([BioUnit, "PolyID"])["PolyID"].transform('size')  ## How many times the poly appears for a replicate/mouse's samples..
df["Prevalence"] = df.groupby("PolyID")["Day"].transform(
    len)  ### How many times each poly is detected in the dataset....

df["CumFreq"] = df.groupby([BioUnit, "PolyID"])["Freq"].transform('sum')  ### useless????
df["PolyCount"] = df.groupby(["SampleID", "Chrom"])["PolyID"].transform(len)  ### How many polys a sample has....

df["MeanFreq"] = df.groupby([BioUnit, "PolyID"])["Freq"].transform('mean')  ### Mean poly frequency in replicate/mouse's samples
df['GeneName'] = df["GeneName"].replace({"Nan": "Intergenic"})
df['oldGeneName'] = df["oldGeneName"].replace({"Nan": "Intergenic"})
df['GeneLength'] = df["cds"].apply(lambda x: len(x) if isinstance(x, str) else 0)
df['Run'] = df.SampleID.astype(str).apply(lambda x: 'VBCF' if x.startswith("6") or x.startswith("512") or x.startswith(
    "3") else "BSF")  ### mark samples based on the sequencing workflow
vbcf = df.loc[df['Run'] == "VBCF", "PolyID"].unique()
bsf = df.loc[df["Run"] == "BSF", "PolyID"].unique()
onlyVBCF = [snp for snp in vbcf if snp not in bsf]  ### select polys present exclusively in VBCF samples
onlyBSF = [snp for snp in bsf if snp not in vbcf]  ### select polys present exclusively in BSF samples
shared = list(set(vbcf).intersection(bsf))  ## select polys shared between sequencing workflows
def seqPatternFunc(x):
    if x in onlyVBCF:
        return "VBCF"
    if x in onlyBSF:
        return "BSF"
    if x in shared:
        return "shared"

df['SequencingPattern'] = df['PolyID'].apply(lambda x: seqPatternFunc(x))


console.print(f"Writing file with NO filters [u i]{exported_file_name}_nofilter.tsv[/u i].")
df.to_csv(exported_file_name + "_nofilter.tsv", index=None, na_rep="Nan", sep='\t')

console.print(f"\n[b]Filtering data...")

if lda_filter:
    df = df.loc[df.LDAcoef.between(-lda_filter, lda_filter)]
    print("Number of Polymorphisms in the dataset after LDA filtering (-{} <=coef <={}): {}".format(lda_filter,
                                                                                                    lda_filter,
                                                                                                    df.drop_duplicates(
                                                                                                        "PolyID").shape[
                                                                                                        0]))

""" Remove ancestral SNPs..."""

if remove_ancestral:
    ### removes fixed polymorphisms in the ancestral population
    df = df.loc[df.Ancestry != 'fixed']
"""Hard Filter Polys on Reads supporting a site  """
df = df.loc[df["TotalCov"] >= total_support]

"""Hard Filter Polys on Reads supporting the Alternative state  """
df = df.loc[df["AltCov"] >= alt_support]

""" Hard Filter Polys on frequency  """
df = df.loc[df['Freq'] >= freq_filter]

""" Keep only libraries from a Run if needed """
if run_filter != 'all': df = df.loc[df.Run == run_filter]

"""                            Collect coverage information and estimate other parameters in the same dataframe                         """

added_columns = {x:y for x,y in zip(["MetagenomeID","SampleProfile",'MaxFreqTraj', f'{BioUnit}PolyAbundance', 'Prevalence', "DayPrevalence","DayPrevalencePercent",f'{BioUnit}Prevalence',f'{BioUnit}PrevalencePercent', 'CumFreq',
                                     'MeanFreq', "DayMeanFreq",f'{BioUnit}MeanFreq',
                                     'PolyCount' ,
                                     "GeneName", "oldGeneName", "GeneLength",  'cds',"cddID", 'cddName','cddDescription', "CogID",  "CogFamily", "CogFamilySize", "CogSize",
                                     "Run", "SequencingPattern", "RunPrevalence"],
                                    ["If isolate the SampleID of the metagenome it derived from else same as SampleID","isolate or metagenome",f"Max Frequency of poly in the {BioUnit}",f'How many times the poly appears in the {BioUnit}',
                                     f'Poly prevalence in the dataset', "Poly's prevalence within Day", "Poly's prevalence within Day (% of samples)",
                                     f"Poly's prevalence within {BioUnit}",f"Poly's prevalence within {BioUnit} (% of samples)",
                                     "Poly's cumulative frequency over time",
                                     "Poly's mean frequency in the dataset ","Poly's mean frequency within timepoint",f"Poly's mean frequency with {BioUnit}",
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

df["MaxFreqTraj"] = df.groupby([BioUnit, "PolyID"])["Freq"].transform('max')  ### Max frequency per replicate/mouse
df[f"{BioUnit}PolyAbundance"] = df.groupby([BioUnit, "PolyID"])["PolyID"].transform(
    'size')  ## How many times the poly appears for a replicate/mouse's samples..
df["CumFreq"] = df.groupby([BioUnit, "PolyID"])["Freq"].transform('sum')  ### useless????
df["DatasetPrevalence"] = df.groupby("PolyID")["Freq"].transform(len)  ### How many times each poly is detected in the dataset....
df["DatasetPrevalencePercent"] = df['DatasetPrevalence'] / len(df.SampleID.unique())

## make functions that create mean and prevalence columns per grouping

add_group_poly_prevalence_and_mean(df, column="Day")
add_group_poly_prevalence_and_mean(df, column='Mouse')
add_group_poly_prevalence_and_mean(df,column='Run')


df["PolyCount"] = df.groupby(["SampleID", "Chrom"])["PolyID"].transform(len)  ### How many polys a sample has....
df["MeanFreq"] = df.groupby(["PolyID"])["Freq"].transform("mean")  ### Mean poly frequency in the dataset



print(f"Dumping the dataframe as '{exported_file_name}_masterTable.obj'.")
console.print(f"Exporting [i u]{exported_file_name}.tsv[/i u]...\n")
pickle.dump(df, open(f"{exported_file_name}_masterTable.obj", 'wb'))
df.to_csv(exported_file_name + '.tsv', index=None, na_rep="Nan", sep='\t')

final_N_samples = df.drop_duplicates("SampleID").shape[0]
final_N_isolates = df.drop_duplicates('SampleID').loc[df.SampleProfile == 'isolate'].shape[0]
final_N_metagenomes = df.drop_duplicates('SampleID').loc[df.SampleProfile == 'metagenome'].shape[0]
table3 = Table(title="Poly numbers")
for x in ['Stage','Over Poly Number','Poly Catalog', 'Samples', "Metagenomes", 'Isolates']:
    table3.add_column(x)
table3.add_row("Start", f'{initial_enties}', f'{initial_different_polys}', f'{initial_N_samples}', f'{N_metagenomes}',f'{N_isolates}')
table3.add_row("End", f'{df.shape[0]}', f'{len(df["PolyID"].unique())}', f'{final_N_samples}', f'{final_N_metagenomes}',f'{final_N_isolates}')

console.print(table3)
print('\n\nDone!')