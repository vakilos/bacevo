import os
import re
import sys

import pandas as pd

from bacevo_fun import *

"""
Directory --> ~/Downloads/isolate_assemblies/annotations
contains folders with annotations for each isolate plus fasta file with the assembly.
"""

parser = argparse.ArgumentParser(description='Starts pipeline for identification of phase variable regions'
                                             'by comparing gene orientation between genome assemblies and reference genome.')
parser.add_argument('-d', dest='folder', help='Path to directory with genome assemblies and annotation files (per sample subdirectories).')
parser.add_argument('-sv', dest='invivo_samples', help='Path to file with in vivo SampleID - sequencing files associations.')
parser.add_argument('-st', dest='invitro_samples', help='Path to file with in vitro SampleID - sequencing files associations.')
parser.add_argument('-a', dest='anno', help='Full path to AnnotationTable.txt',  default='AnnotationTable.txt')
parser.add_argument('-g', dest='genes', help='Path to name of a blast database for gene sequences.', default='genomic_genes.fa')
parser.add_argument('-o', dest='out', help='Directory to export files. Default is CWD', default=os.getcwd())
args = parser.parse_args()

pd.set_option("display.max_columns", None)


sv_samples = pd.read_csv(args.invivo_samples, dtype={'SampleID':str})
st_samples = pd.read_csv(args.invitro_samples)

exp1 = re.compile(r'BSF_\d+_\w{9}_\d_(.+).bam')
exp2 = re.compile(r'(.+)(?=_spades)')

## maybe needs to be updated
sv_samples['LiSC_path'] = sv_samples['LiSC_path'].apply(lambda x: re.match(exp1, x.split("/")[-1]).group(1) if re.match(exp1, x.split("/")[-1]) else x)
st_samples['LiSC_path'] = st_samples['LiSC_path'].apply(lambda x: os.path.basename(x))

os.chdir(args.folder)
idx = 0
ref_db  = args.genes

def blastn_to_ref_genes(ref_genes_db = '', query='', out= ''):
    os.system(f"blastn -db {ref_db} -query {query} -out {out}  -outfmt 7")
    blastDf = pd.read_csv(f"{folder}/blastnResults.out", sep="\t", comment="#", header=None)
    blastDf.rename(columns={0:"GeneID",1:"RefBestHit",2:"RefIdentity",3:"Alignment Length",4:"mismatches",5:"gap opens", 6:'qstart',7:'qend',8:'sstart',9:'send',10:"Blast_E-value",11:"Blast_Bit-score"},inplace=True)
    blastDf = pd.DataFrame(blastDf.sort_values(["GeneID","Blast_E-value"]).groupby("GeneID",as_index=False).first()) ## extracts the best hit (based on E-value) for each query
    blastDf["GeneID"] = blastDf["GeneID"].apply(lambda x: x.split("#")[0].strip())
    blastDf = blastDf[["GeneID","RefBestHit", "RefIdentity","Blast_E-value"]]
    return blastDf
newdir_v = f'{args.out}/invivo_isolates_assemblies'
newdir_t = f'{args.out}/invitro_isolates_assemblies'

for i,folder in zip(track(range(len(glob.glob("*")))),glob.glob("*")):
    newdir = None
    new_id = None
    try:
        #print(folder)
        if re.match(exp2, folder).group(1) in sv_samples['LiSC_path'].values:  ## selects each folder based on their bamFileName - switch or rename with NewFileName
            idx += 1
            print(re.match(exp2, folder).group(1))
            new_id = sv_samples.loc[sv_samples['LiSC_path'] == re.match(exp2, folder).group(1), 'SampleID'].values[0]
            newdir = newdir_v
    except:
        pass
    try:
        if re.match(re.compile(r'(\d+)_.*'), folder):
            if re.match(re.compile(r'(\d+)_.*'), folder).group(1) == '35254':
                new_id = "35254"
                newdir = newdir_v
                idx += 1
    except:
        pass
    if re.match(re.compile(r'(\d+)_.*'), folder):
        if re.match(re.compile(r'(\d+)_.*'), folder).group(1) == '69383':
            new_id = "69383"
            newdir = newdir_t
            idx += 1

    try:
        if re.match(exp2, folder).group(1) in st_samples['LiSC_path'].values:
            idx += 1
            new_id = st_samples.loc[st_samples["LiSC_path"] == re.match(exp2, folder).group(1), 'SampleID'].values[0]
            newdir = newdir_t
    except:
        pass
    ##if os.path.exists(f'{folder}/blastnResults.out'):
    if new_id: ## run this only for valid samples
        if glob.glob(f"{folder}/*blastnResults.out"):
            os.system(f"rm {folder}/*blastnResults.out")
        blastDf = blastn_to_ref_genes(ref_genes_db=ref_db, query=f'{folder}/{folder}.ffn', out=f'{folder}/blastnResults.out')
        blastDf['SampleID'] =  [new_id for x in blastDf.index]

        ## add contig info from respective gff files
        gff = pd.read_csv(f'{folder}/{folder}.gff', comment="#", header=None, sep='\t', dtype=str).rename(
            columns={0: "nodeID", 8: "GeneID"})
        gff = gff.loc[gff.nodeID.str.startswith("NODE")]  ## removes lines with DNA sequences in the end of file.
        gff = gff[["nodeID", 'GeneID']]
        expr = re.compile(r'NODE_(\d+)_length_\d+_cov_\d+\.*\d+')
        gff.nodeID = gff.nodeID.apply(lambda x: f"S{new_id}n{re.match(expr, x).groups()[0]}")  ### defines the new node_ids {SAMPLEID}nNODE_NUMBER
        gff.GeneID = gff.GeneID.apply(lambda x: re.match(re.compile(r"ID=(\w+);.*"), x).groups()[0])

        blastDf = blastDf.merge(gff, on ='GeneID', how='left')

        blastDf.to_csv(f'{folder}/{new_id}_blastnResults.csv', index=False)
        os.system(f'rm {folder}/blastnResults.out')
        if not os.path.exists(f'{newdir}/{new_id}'):
            os.system(f'mkdir -p {newdir}/{new_id}')
        os.system(f'rsync -a {folder}/* {newdir}/{new_id}')
        os.system(f"rename '{folder}' {new_id} {newdir}/{new_id}/*")



