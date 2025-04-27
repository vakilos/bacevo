
import os
import sys

import pandas as pd
from debugpy.common.json import default
from numpy.core.defchararray import endswith
from polars import read_excel

from bacevo_fun import *
parser = argparse.ArgumentParser(description="""Call from the dNdS directory to concatenate files from metagenomic populations (after running dnds_calculate.py).
 It adds SampleID column to the final file.""")
parser.add_argument('-n', dest='name', help='Name tag to used in the exported files', default=False)
parser.add_argument('-o', dest='outdir', help='Name of the output directory. Default is parent directort of dNdS', default=os.path.dirname(os.getcwd()))
parser.add_argument('-a', dest ='annotable', help='FinalPolyTable.tsv used in the project', default=os.path.dirname(os.getcwd())+"/AnnotationTable.txt")
parser.add_argument('-m', dest='meta', help='Meta table for extra annotation columns for genes.', default=os.path.dirname(os.getcwd())+"/All_genes_polycount_norm_incl_PUL_table_forOlga.csv")
args = parser.parse_args()
## run in the dNdS directory
suffixes = ["codon_dNdS.csv", 'gene_dNdS.csv', 'site_pi.csv'] #'codon_pi.csv' 'site_pi.csv','gene_pi.csv'

meta = pd.read_csv(args.meta, sep='\t',usecols=['GeneName', 'oldGeneName', "PUL","Function",'EC','Overall function','YesNoPUL','Substrate general','Substrate specific','product', 'cddID','cddName','cddDescription']).drop_duplicates('GeneName')
ano = pd.read_csv(args.annotable, sep='\t', usecols=['locus_tag', 'cds','CogID','CogAnnotation', 'CogFamily']).drop_duplicates(['locus_tag','CogID','CogFamily']).rename(columns={'locus_tag':'GeneName'})


ano = ano.merge(meta, on='GeneName', how='left')


if not args.outdir.endswith("/"):
    args.outdir += '/'

if args.name:
    if not args.name.endswith("_"):
        args.name += "_"
else:
    args.name = ""
for suf in suffixes:
    cnt = 0
    cnt_d = 0
    out_df = pd.DataFrame()
    os.chdir(os.getcwd())
    for x in glob.glob(f'*{suf}'):
        this_df = pd.read_csv(x, dtype={'SampleID':str, "GeneName":str})
        out_df = pd.concat([out_df, this_df])

    out_df.merge(ano, on='GeneName',how='left').to_csv(args.outdir+args.name+"full_annotation_"+suf, index=False)
    out_df.merge(ano.drop_duplicates("GeneName"), how='left', on='GeneName').to_csv(args.outdir+args.name+suf, index=False)



