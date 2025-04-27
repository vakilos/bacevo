import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import rich

from bacevo_fun import *
from ai_fun import *
console = rich.console.Console(color_system='windows')
parser = argparse.ArgumentParser(prog='ai_reverse_nodes.py',description='Reverse contig orientation according to portion of genes that share same orientation in the reference')
parser.add_argument('-d', dest='samplesDir',help='Directory with sample folders with all blast hits. Default is CWD', default=os.getcwd())
parser.add_argument('-o', dest='outDir',help='Directory to export files..', default=f'{os.getcwd()}/gene_order')

args = parser.parse_args()
sampledir = args.samplesDir
outdir = args.outDir
os.chdir(sampledir)
#binStatsFile = args.binstatFile


if not os.path.exists(f'{outdir}/contigs_cov_length.csv'):
    bin_stats = pd.DataFrame()
    for _,folder in zip(track(range(len(glob.glob("*"))),description="## Preparing files..."),glob.glob("*")):

            if not folder.startswith("NC"): ## runs only for isolates (not for ref)
                bin_stats = pd.concat([bin_stats, create_bin_stats(SeqIO.parse(f"{folder}/{folder}.fna", 'fasta'), folder)])
            bin_stats.to_csv(f'{outdir}/contigs_cov_length.csv',index=False)
else:
    bin_stats = pd.read_csv(f'{outdir}/contigs_cov_length.csv')
"""
ALL .ffn files have genes encoded as if in + strand!!!!!!!
"""

console.print('\n\n## Parsing best hits from blastn results...')
df =pd.concat([pd.read_csv(f'{f}/{f}_blastnResults.csv') for f in glob.glob("*") if not f.startswith('NC')])
target_genes_best_hits = df[['GeneID', 'RefBestHit']].set_index("GeneID").to_dict()['RefBestHit']
lc = pd.DataFrame()
for i,folder in zip(track(range(len(glob.glob("*"))),description="## Locating genes on contigs..."),glob.glob("*")):
    if not folder.startswith("NC"):
        lf = parse_gff(f'{folder}/{folder}.gff')
        if not lf.empty:
            fnn = fasta_to_df(f'{folder}/{folder}.ffn')
            lf['DNA'] = lf['GeneID'].apply(lambda x: fnn.loc[fnn.Chrom == x,'seq'].item())
            lf['RefBestHit'] = lf['GeneID'].map(target_genes_best_hits)
            lf.sort_values('nodeID',inplace=True)
            #node_dict = {k:v for k,v in zip(lf.nodeID.values, range(len(lf.nodeID.values)))}
            #lf['nodeN'] = lf.nodeID.map(node_dict) ### stores indexes for nodes associated with target genes
            ## check if there is only one nodeID - target genes should be positioned on the same contig/node
            lf['GeneLength'] = lf.apply(lambda x: np.abs(int(x.start) - int(x.end)) +1,axis=1)
            lc = pd.concat([lc,lf])

    else:
        """The AnnotationTable.txt genomic coordinates are used for slicing the refseq. Thus
        Add +1 to the gene end coord and 0 to gene-start difference when calculating gene length"""
        lf = pd.read_csv(f'{folder}/AnnotationTable.txt',sep='\t').rename(columns={"Chrom":'nodeID','locus_tag':"GeneID"})
        lf = lf.drop_duplicates("GeneID")
        lf['SampleID'] = ['RefSeq' for x in lf.index]
        lf = lf.rename(columns ={'cds':'DNA'})
        lf['RefBestHit'] = lf.GeneID.values
        #lf['nodeN'] = [0 for x in lf.index]
        lf['GeneLength'] = lf.apply(lambda x: np.abs(int(x.start) - int(x.end)),axis=1)
        lf['strand'] = lf['strand'].apply(lambda x: "-" if x == -1 else "+")
        lf['end'] = lf['end'].astype(int) - 1
        lf = lf[['nodeID','GeneID','start','end','strand', 'GeneLength', "SampleID", "RefBestHit",'DNA']] ##"DNA" 'nodeN'
        lc = pd.concat([lc,lf])




lc = check_gene_order(lc)
"""
the gene coordinates for 
genes with OriWithRef == '-' (reversed contigs) should be 
LengthOfContig - coord + 1
"""


lc = lc.merge(bin_stats.drop(columns='SampleID'), on='nodeID', how='left').fillna(0)
lc['start'] = lc.apply(lambda x: int(x.length) - x.start + 1 if x.OriScore < 0.5 else x.start, axis=1)
lc['end'] = lc.apply(lambda x: int(x.length) - x.end + 1 if x.OriScore < 0.5 else x.end, axis=1)
lc.to_csv(f"{outdir}/all_coords.csv", index=False)


