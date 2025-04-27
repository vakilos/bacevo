import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import rich

from bacevo_fun import *
from ai_fun import *
console = rich.console.Console(color_system='windows')
tprint("\n\nGene Order (I)")
parser = argparse.ArgumentParser(prog='ai_gene_order.py',description='Compare gene order of homologous genomic regions')
parser.add_argument('-d', dest='basedir',help='Project directory. Default = current working directory', default=os.getcwd())
parser.add_argument('-o', dest='outDir',help='Directory to export files..', default=os.getcwd()+'/gene_order')
parser.add_argument('-e', dest='experiment', help='Choose between invivo and invitro experiment', choices=['invivo','invitro'], default='invivo')
parser.add_argument('-g', dest='anno', help='Path to AnnotationTable.txt. Default = CWD', default = os.getcwd()+"/AnnotationTable.txt")
parser.add_argument('-a', dest='ancestor', help='ID for ancestral library', default='35254', type=str)
args = parser.parse_args()


outdir = args.outDir



"""

  Read the all_coords.csv file from ai_reverse_nodes.py located in /gene_order folder.


1. ai_gene_order.py
2. assembly-isolates-gene_alignment.py 
"""
####################### parameters ##########################################
oneNode = False
targetMode = True
n50_filter = 100000
experiment = args.experiment
ancestorID = "S"+args.ancestor
#############################################################################
df = pd.read_csv(f"{args.basedir}/gene_order/all_coords.csv")
meta = pd.concat([pd.read_csv(x, sep='\t', dtype='str') for x in glob.glob(f'/home/christos/bacevo/{experiment}/*exp*')]) ### this should be changed according to project
meta['SampleID'] = ("S"+meta['SampleID'])
df = df.merge(meta, on='SampleID', how='left')
sns.histplot(df.N50)
plt.show()
all_libraries = len(df.SampleID.unique())
df = df.loc[((df.N50.ge(n50_filter)) & (df.SampleID != 'RefSeq')) | (df.SampleID == 'RefSeq')]
good_assemblies = len(df.SampleID.unique())
console.print(f'{good_assemblies}/{all_libraries} filtered as good quality assemblies.')
############################################################################################################

tf = pd.read_csv(f'{outdir}/target_genes.csv')
annotation = args.anno
def compare_homologous_target_genes(df:pd.DataFrame, targets:list, geneGroup:str, annotation:str, outdir:str, ancestorID=ancestorID):
    """

    :param df: output from ai_reverse_nodes.py
    :param targets: list of RefSeq gene IDs
    :param annotation: AnnotationTable.txt
    :return:
    """

    df = df.loc[df.RefBestHit.isin(target_genes)]

    oneNodeSamples = []

    for name, group in df.groupby("SampleID"): ## remove samples with targets in more than one contig
        if len(group.nodeID.unique()) == 1:
            oneNodeSamples += [name]
    if oneNode:
        df = df.loc[df.SampleID.isin(oneNodeSamples)]

    df = extract_relative_gene_coords(df, mode=experiment)
    #print(df.loc[(df.SampleID == 'S51015'), ['RefBestHit', 'generelstart','generelend']])
    ano = pd.read_csv(annotation,
                      sep='\t', usecols=['old_locus_tag','locus_tag','product']).set_index('locus_tag')
    ano_names  = ano.to_dict()['old_locus_tag']
    ano_info = ano.to_dict()['product']


    df['string'] = ["" for x in df.index]
    for name,group in df.groupby('SampleID'):
        group['string'] = group.RefBestHit +"["+group.nodeN.astype(str)+"]"+"{"+group.generelstart.astype(str)+":"+group.generelend.astype(str)+"}"+\
                          "("+group.genelength.astype(str)+"bp)"+"{"+group.genestrand+"}"
        string ="| ".join(group['string'].values)
        df.loc[df.SampleID == name, 'string'] = string


    #### remove redundant variants (duplicates)
    keepSamples = []
    strings = []
    for name, group in df.loc[df.SampleID != "RefSeq"].groupby('SampleID'):
        if group.string.values[0] not in strings:
            strings += [group.string.values[0]]
            keepSamples += [name]
    keepSamples += ["RefSeq"] ## always add the reference
    #df['string'] = df['SampleID'] + "::: "+ df['string']
    ## create labels for variants
    df['VariantName'] = ['' for x in df.index]
    vindex = 1
    for name, group in df.groupby('string'):
        if 'RefSeq' in group.SampleID.values and ancestorID in group.SampleID.values:
            df.loc[df.string == name, 'VariantName'] = "RefSeq / Ancestral variant"
        elif 'RefSeq' in group.SampleID.values:
            df.loc[df.string == name, 'VariantName'] = "RefSeq variant"

        elif ancestorID in group.SampleID.values:
            df.loc[df.string == name, 'VariantName'] = "Ancestral variant"
        else:
            df.loc[df.string == name, 'VariantName'] = "V" + str(vindex)
            vindex += 1

    freq_dict = df.loc[df.SampleID != "RefSeq"].drop_duplicates("SampleID").value_counts('string',normalize=True).to_dict()

    #print(mouse_sizes)
    if experiment == 'invivo':
        mouse_sizes = df.drop_duplicates("SampleID").groupby("Mouse").size().to_dict()
        freq_mouse = df.drop_duplicates("SampleID").groupby("Mouse").agg({"VariantName":'value_counts'}).rename(columns={'VariantName':"Count"}).reset_index()
        freq_mouse['size'] = freq_mouse['Mouse'].map(mouse_sizes)
        freq_mouse['Freq'] = freq_mouse['Count'] / freq_mouse['size']
        freq_mouse.to_csv(f'{outdir}/{geneGroup}_variants_frequency_per_mouse.csv')
    else:
        mix_sizes = df.drop_duplicates("SampleID").groupby("Mix").size().to_dict()
        freq_mix = df.drop_duplicates("SampleID").groupby("Mix").agg({"VariantName":'value_counts'}).rename(columns={"VariantName":'Count'}).reset_index()
        freq_mix['size'] = freq_mix['Mix'].map(mix_sizes)
        freq_mix['Freq'] = freq_mix["Count"] / freq_mix['size']
        freq_mix.to_csv(f'{outdir}/{geneGroup}_variants_frequency_per_mix.csv')

    df['VariantFreq'] = df['string'].map(freq_dict).fillna(0)

    df = df.loc[df.SampleID.isin(keepSamples)] ### this will keep representative (for variants) samples only


    df['oldName'] = df.RefBestHit.map(ano_names)
    df['info'] = df.RefBestHit.map(ano_info).fillna('no info')

    df.to_csv(f'{outdir}/{geneGroup}_relative_coords.csv', index=False)
    return df

for name, group in tf.groupby("GeneGroup"):
    print(name)
    target_genes = group['GeneName'].values
    gf = compare_homologous_target_genes(df, targets = target_genes, geneGroup=name,annotation = annotation, outdir=outdir)
    #print(df.info())










