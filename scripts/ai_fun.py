import pandas as pd
import numpy as np
import os,sys,re, logging
import matplotlib.pyplot as plt
from rich.console import Console
from rich.progress import track
from art import *
import argparse
from Bio import SeqIO
import seaborn as sns
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastpCommandline
def create_bin_stats(fasta, sampleid:str):
    """
    :param fasta is a SeqIO parsed iterator
    :param sampleid is a string for sample ID of the fasta file
    :return a dataframe
    """
    #!/usr/bin/env python3
# -*- coding: utf-8 -*-

    def calculate_N50(list_of_lengths):
        """Calculate N50 for a sequence of numbers.

        Args:
            list_of_lengths (list): List of numbers.

        Returns:
            float: N50 value.

        """
        tmp = []
        for tmp_number in set(list_of_lengths):
                tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
        tmp.sort()

        if (len(tmp) % 2) == 0:
            median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
        else:
            median = tmp[int(len(tmp) / 2)]

        return median


    df = pd.DataFrame()
    for rec in fasta:
        expr = re.compile(r'NODE_(\d+)_length_(\d+)_cov_(\d+\.*\d+)')
        node, length, cov = re.match(expr, rec.id).groups()
        df = pd.concat([df,pd.DataFrame({"SampleID":["S"+sampleid], "length":[length], 'cov':[cov],
                                         'nodeID':["S"+sampleid+"n"+node]})])

    n50  = calculate_N50(df['length'].astype(int).tolist())
    df['N50'] = [n50 for x in df.index]
    return df
def rename_fasta_bins(fasta, folder):
    """
    :param fasta is a SeqIO parsed iterator
    """
    sequences = []
    for rec in SeqIO.parse(fasta, 'fasta'):
        expr = re.compile(r'NODE_(\d+)_length_(\d+)_cov_(\d+\.*\d+)')
        node, length, cov = re.match(expr, rec.id).groups()
        rec.id = "S"+folder+"n"+node
        rec.name = ''
        rec.description = ''
        sequences.append(rec)
    with open(f'{fasta.replace(".fna","_sib.fna")}', 'w') as export:
        SeqIO.write(sequences, export, 'fasta')

def parse_gff(filepath):
    gff = pd.read_csv(f"{filepath}", comment="#", sep='\t', header=None, dtype=str).rename(columns={0:'nodeID',3:"start",4:"end",6:"strand",8:'GeneID'})
    gff = gff.loc[gff['nodeID'].str.startswith("NODE")]  ## removes lines with DNA sequences in the end of file.
    expr = re.compile(r'NODE_(\d+)_length_\d+_cov_\d+\.*\d+')
    gff['nodeID'] = gff['nodeID'].apply(lambda x: f"S{os.path.basename(filepath).replace('.gff','')}n{re.match(expr, x).group(1)}")  ### defines the new node_ids {SAMPLEID}nNODE_NUMBER
    gff['GeneID'] = gff['GeneID'].apply(lambda x: re.match(re.compile(r"ID=(\w+);.*"), x).group(1))
    gff['SampleID'] = [f"S{os.path.basename(filepath).replace('.gff','')}" for x in gff.index]
    gff = gff[['nodeID','GeneID','start','end','strand','SampleID']]
    return gff



def relative_coords(df: pd.DataFrame()):
    """
    :param df: A dataframe with Block coords and gene coords - Run this function after locate_genes_in_synteny() function
    In the input dataframe BlockStart represents the leftmost (5'-end) compared to BlockEnd - The original block orientation is stored as BlockOrientation
    :return: adds RelStart and RelEnd columns
    """
    of = pd.DataFrame()
    idx = 1
    #df['BlockOrientation'] = ["+" for x in df.index]
    for name, group in df.groupby("nodeID"):
        group.sort_values('BlockStart', inplace=True)
        leftmost_coord = group.BlockStart.values[0] ## creates 0-index coords based on the leftmost (5'-end) of the block on the contig -- +1 if need 1-index
        ### create index based on reference gene coordinate..
        group["RelStart"] = group.BlockStart - leftmost_coord
        group['RelEnd'] = group.BlockEnd - leftmost_coord
        group["nodeID"] = [name for x in group.index] ## I need this - because it is dropped from the original dataframe..(with groupby)
        #group['Yval'] = [idx for x in group.index]
        of = pd.concat([of, group])
        idx += 1

    return of
def relative_coords_2(df: pd.DataFrame):
    '''
    length columns corresponds to node length
    :param df:
    :return:
    '''

    of = pd.DataFrame()
    sample_idx = 0
    df['RelStart'] = [None for x in df.index]
    df['RelEnd'] = [None for x in df.index]
    for sampleID, group in df.groupby("SampleID"):
        sample_idx += 1
        node_idx = 0
        group.sort_values(['nodeID','BlockStart'], inplace=True)
        distance_from_right_edge = 0
        for nodeID, node in group.groupby('nodeID'):
            node.sort_values('BlockStart', inplace=True)
            node_idx += 1
            leftmost_coord = node.BlockStart.values[0] ## creates 0-index coords based on the leftmost (5'-end) of the block on the contig -- +1 if need 1-index
            for i, block in node.iterrows():
                ### create index based on reference gene coordinate..
                ### distance_from_right_edge only affects output when genes span more than one contig
                ### because blocks can be associated with more than one gene / regions add the Block BlockStart and BlockEnd slicing df
                df.loc[(df.SampleID == sampleID) & (df.nodeID == nodeID) &
                       (df.Block == block.Block) &
                       (df.BlockStart == block.BlockStart) & (df.BlockEnd == block.BlockEnd)
                ,['RelStart','RelEnd']] = [distance_from_right_edge + block.BlockStart - leftmost_coord,
                                                                 distance_from_right_edge + block.BlockEnd - leftmost_coord]
                #block["RelStart"] = distance_from_right_edge + block.BlockStart - leftmost_coord
                #block['RelEnd'] = distance_from_right_edge + block.BlockEnd - leftmost_coord
                #block["nodeID"] = [nodeID for x in block.index]
                #block['SampleID'] = [sampleID for x in block.index]
                #group['Yval'] = [idx for x in group.index]
                #of = pd.concat([of, block])
            distance_from_right_edge = node.nodelength.astype(int).values[-1] - node.BlockEnd.values[-1]

    return df
def add_gene_relative_coords(df: pd.DataFrame):
    """
    For all gene entries, Yval = -10

    :param df: A dataframe which already has RelStart RelEnd columns
    ( from relative_coords function)
    :return: adds an entry for gene relative coordinates to the input dataframe
    """

    for name, group in df.groupby('nodeID'):
        if name == 'NC_004663.1':
            group.sort_values('BlockStart',inplace=True)
            leftmost_coord = group.GeneStart.values[0] ##### should this be block start or gene start
            for geneID in group.drop_duplicates("GeneID")['GeneID'].values:
                gene_start = group.loc[group.GeneID == geneID, 'GeneStart'].values[0]
                gene_end = group.loc[group.GeneID == geneID, 'GeneEnd'].values[0]

                gene_strand = group.loc[group.GeneID == geneID, 'Strand'].values[0]
                # create relative coordinates also for target gene..
                gene_rel_start = gene_start - leftmost_coord
                gene_rel_end = gene_end - leftmost_coord
                gene = group.head(1).copy()
                gene['RelStart'] = gene_rel_start
                gene['RelEnd'] = gene_rel_end
                gene["Yval"] = -10
                gene['nodeID'] = geneID
                gene['Block'] = geneID
                gene['GeneID'] = geneID
                gene['Strand'] = gene_strand
                gene['SampleID'] = 'genecoord'
                gene['RefBestHit'] = geneID

                df = pd.concat([df, gene])
    return df

def add_ref_coords_of_block(gff):
    gff['RefStart'] = [np.nan for x in gff.index]
    gff['RefEnd'] = [np.nan for x in gff.index]
    for _,x in zip(track(range(len(gff.groupby("Block").groups))),gff.groupby("Block")):
        name, group = x
        if 'NC_004663.1' in group.nodeID.values:
            refStart = group.loc[group.nodeID == 'NC_004663.1','Start'].values[0]
            refEnd = group.loc[group.nodeID == 'NC_004663.1','End'].values[0]
            for i, row in group.iterrows():
                gff.loc[i, "RefStart"] = refStart
                gff.loc[i, "RefEnd"] = refEnd
def add_blockstring(df:pd.DataFrame):
    df['BlockString'] = ["" for x in df.index]
    for name, group in df.groupby('nodeID'):
        group['blockString'] = "L"+group.BpInBlock.astype(str) + "B" + group.Block.astype(str) + group.BlockOrientation
        string = "".join(group['blockString'].values)
        df.loc[df.nodeID == name, 'BlockString'] = string
    return  df

def hclust_samples(df, targets:list):
    """
    removes any target genes from the Yval calculations -
    they should be calculated separately downstream anyway

    input dataframe has duplicated enties for blocks that span more than one gene - so I need to sum column bpInBlock to include all entries
    this is stored in TotalBpInBlock column - must iterate over target genes -> over blocks for each node
    :param df:
    :param targets:
    :return:
    """

    df['TotalBpInBlock'] = df.groupby(['nodeID','Block'])['BpInBlock'].transform(sum) ### sums up bp contribution from all genes spanning the block

    #print(df.drop_duplicates(['nodeID','Block']).shape)
    #print(df[df.duplicated(['nodeID','Block'], keep=False)].sort_values(['nodeID','Block']))
    blocks = [gene+"B"+str(block) for block in df.Block.unique() for gene in targets if block not in targets] ## exclude the gene entries added for reference visualization
    nodes = [x for x in df.nodeID.unique() if x not in targets]
    #bl_matrix = np.zeros(shape=[len(nodes),len(blocks)])
    bl = pd.DataFrame(columns=blocks, index=nodes)
    for node in nodes:
        for gene in targets:
            for blk in df.Block.unique():
                if blk in df.loc[(df.nodeID == node) & (df.GeneID == gene),"Block"].values:
                    bl.loc[node, gene+"B"+str(blk)] = df.loc[(df.nodeID == node) & (df.GeneID == gene) & (df.Block == blk), 'TotalBpInBlock'].item()
                else:
                    bl.loc[node, gene+"B"+str(blk) ] = 0
    from skbio.stats.distance import DistanceMatrix
    from sklearn.preprocessing import StandardScaler
    data_scaler = StandardScaler()
    scaled_data = data_scaler.fit_transform(bl)
    from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
    #dm = DistanceMatrix(scipy.spatial.distance.pdist(bl)).data
    average_clustering = linkage(scaled_data, method="complete", metric="euclidean")
    #dm = pd.DataFrame(index=nodes, columns=nodes, data=scaled_data)
    #cluster_labels = cut_tree(average_clustering, n_clusters=2).reshape(-1, )
    #print(cluster_labels)
    #bl["Cluster"] = cluster_labels
    dendro_idx = dendrogram(average_clustering, orientation='left')['leaves']
    ordered_samples = [nodes[i] for i in dendro_idx]
    dendro_order = {k:v for k,v in zip(ordered_samples,range(len(ordered_samples)))}
    df['Yval'] = [0 for x in df.index]
    df['Yval'] = df['nodeID'].map(dendro_order)
    return df

def locate_genes_in_synteny(lc:pd.DataFrame, sb:pd.DataFrame) -> pd.DataFrame:
    """

    :param lc: dataframe with parsed data from gff / AnnotationTable
    :param sb: dataframe with block info from sibeliaz output (synteny_blocks.csv from OLDassembly-isolates-sibeliaz.py) or maf2synteny output parsed as
    analysis_out.csv from assembly-isolates-stat.py
    :return: sb_hits
    """

    sb_hits = pd.DataFrame()
    for _,x in zip(track(range(len(sb.index)), description='## Locating target gene(s) in synteny blocks'), sb.iterrows()):
        i,block = x

        block_start_on_contig = int(block.start)
        block_end_on_contig = int(block.end)
        block_orientation = "+"
        block_strand = block.strand
        if not block_start_on_contig < block_end_on_contig:
            block_start_on_contig = int(block.end)
            block_end_on_contig = int(block.start)
            block_orientation = '-'

        this_lc = lc.loc[(lc.nodeID == block.nodeID)] #### subset based on block's contig
        #if this_lc.shape[0] > 1: logging.warning(f'This block spans {this_lc.shape[0]} genes!!')

        ## what if block is associated with more than one gene??? need to iterate over all spanning genes = This creates  duplicated block entries - each for each gene it spans - for plottin remove those
        for ii, gene in this_lc.iterrows(): ## iterate over target genes on this node

            gene_start = int(gene.start)
            gene_end = int(gene.end)
            gene_strand = gene.strand
            gene_id = gene.GeneID
            gene_refbest = gene.RefBestHit
            gene_length = np.abs(gene_end - gene_start)
            sampleID = gene.SampleID
            """Start and end in this context are defined based on the site index or 5' end NOT based on CDS"""
            if not gene_start < gene_end: ## this should not happen
                gene_start = lc.loc[(lc.nodeID == block.nodeID)]['end'].astype(int).item()
                gene_end = lc.loc[(lc.nodeID == block.nodeID)]['start'].astype(int).item()

            if gene_start >=  block_start_on_contig and gene_end <= block_end_on_contig:
                #print(gene_length, dif, bp_in_block , gene_start,gene_end, row.start, row.end, row.length)
                #assert bp_in_block == row.length
                logging.info(f"\tFound complete target gene {gene_id} in block {block.Block}")
                h = pd.DataFrame({'GeneID':[gene_id], "RefBestHit":[gene_refbest],'GeneStart':[gene_start],'GeneEnd':[gene_end],'Strand':[gene_strand], "GeneLength":[gene.GeneLength],
                                   "nodeID":[block.nodeID], 'nodeN':[gene.nodeN],"Block":[block.Block],'BlockStart':[block_start_on_contig],'BlockEnd':[block_end_on_contig],
                                  "BlockSize":[block.length], "BlockOrientation":[block_orientation], "BlockStrand":[block_strand],
                                   "FractionInBlock":[100.0],'BpInBlock':[gene_length],'Completeness':['complete gene in block'], "SampleID":[sampleID]})
                sb_hits = pd.concat([sb_hits, h])

            elif gene_start < block_start_on_contig and gene_end > block_end_on_contig: ### block is within gene
                dif = gene_end - block_end_on_contig + block_start_on_contig - gene_start
                bp_in_block = abs(gene_length - dif)
                portion_in_block = round((bp_in_block)*100 / gene_length, 1)

                #print(gene_length, dif, bp_in_block, gene_start,gene_end, row.start, row.end , row.length)
                h = pd.DataFrame({'GeneID':[gene_id], "RefBestHit":[gene_refbest],'GeneStart':[gene_start],'GeneEnd':[gene_end],'Strand':[gene_strand],"GeneLength":[gene.GeneLength],
                                   "nodeID":[block.nodeID], 'nodeN':[gene.nodeN], "Block":[block.Block],'BlockStart':[block_start_on_contig],'BlockEnd':[block_end_on_contig],
                                  "BlockSize":[block.length], "BlockOrientation":[block_orientation],"BlockStrand":[block_strand],
                                   "FractionInBlock":[portion_in_block],'BpInBlock':[bp_in_block],'Completeness':['block within gene'], "SampleID":[sampleID]})
                sb_hits = pd.concat([sb_hits, h])
            elif gene_start <= block_end_on_contig and gene_end > block_end_on_contig:
                dif = gene_end - block_end_on_contig
                portion_in_block = round((gene_length - dif)*100 / gene_length, 1)
                bp_in_block = gene_length - dif
                #assert bp_in_block == row.length
                h = pd.DataFrame({'GeneID':[gene_id],"RefBestHit":[gene_refbest],'GeneStart':[gene_start],'GeneEnd':[gene_end],'Strand':[gene_strand],"GeneLength":[gene.GeneLength],
                                   "nodeID":[block.nodeID], 'nodeN':[gene.nodeN],"Block":[block.Block],'BlockStart':[block_start_on_contig],'BlockEnd':[block_end_on_contig],
                                  "BlockSize":[block.length], "BlockOrientation":[block_orientation],"BlockStrand":[block_strand],
                                   "FractionInBlock":[portion_in_block],'BpInBlock':[bp_in_block],'Completeness':['only gene start in block'], "SampleID":[sampleID]})
                sb_hits = pd.concat([sb_hits, h])

                logging.info(f"\tFound incomplete target gene {gene_id} in block {block.Block} | truncated gene end | {portion_in_block} % in block ({bp_in_block} bp) ")

            elif gene_start < block_start_on_contig and gene_end >= block_start_on_contig:
                #print(gene_length,dif, bp_in_block,  row.start, row.end, row.length)
                dif = block_start_on_contig - gene_start
                portion_in_block = round((gene_length - dif)*100 / gene_length,1)
                bp_in_block = gene_length - dif

                h = pd.DataFrame({'GeneID':[gene_id],"RefBestHit":[gene_refbest],'GeneStart':[gene_start],'GeneEnd':[gene_end],'Strand':[gene_strand],"GeneLength":[gene.GeneLength],
                                   "nodeID":[block.nodeID], 'nodeN':[gene.nodeN],"Block":[block.Block],'BlockStart':[block_start_on_contig],'BlockEnd':[block_end_on_contig],
                                  "BlockSize":[block.length], "BlockOrientation":[block_orientation],"BlockStrand":[block_strand],
                                   "FractionInBlock":[portion_in_block],'BpInBlock':[bp_in_block],'Completeness':['only gene end in block'], "SampleID":[sampleID]})
                sb_hits = pd.concat([sb_hits, h])
                logging.info(f"\tFound incomplete target gene {gene_id} in block {block.Block} | truncated gene start | {portion_in_block} % in block ({bp_in_block}bp)")
    return sb_hits
def how_close_to_the_edge(df:pd.DataFrame, cf:pd.DataFrame, plot=False):
    """
    Determine how close to the end of a node/contig a gene is.
    Group block into informative categories: split --> if block within a gene which contig is smaller than reference gene
    Adds BlockShare column that stores info for how many genes a block is associated with
    Adds fromEdge column that is the min dinstance from nodes's edges
    Adds nodelength that is the length of each node
    :param df: input dataframe with gene coords
    :param cf: input dataframe with node/contig info
    :return: adds column fromEdge column
    """
    #cf['SampleID'] = cf["SampleID"].apply(lambda x: "S"+x)
    df = df.merge(cf.rename(columns={'length':'nodelength'}), on=['SampleID','nodeID'], how='left')
    df['nodelength'] = df['nodelength'].fillna(0)
    df['fromEdge'] = df.apply(lambda x: np.min([int(x.nodelength) - int(x.GeneEnd), int(x.GeneStart)]), axis=1)
    df['BlockShare'] = [0 for x in df.index] ### counts in how many genes a block is involved
    for name, group in df.groupby(['SampleID',"Block"]):
        df.loc[(df.SampleID == name[0]) & (df.Block == name[1]), 'BlockShare'] = len(group.RefBestHit.unique())

    if plot:
        sns.histplot(df['fromEdge'])
        plt.ticklabel_format(style='plain', axis='x')
        plt.show()
    return df

def extract_relative_gene_coords(lc:pd.DataFrame, mode='invivo'):

    gene_coords = pd.DataFrame()
    sample_idx = 0
    lc.reset_index(inplace=True)
    lc = lc.astype({'start':int, 'end':int})
    #lc['start'] = lc['start'].astype(int)
    #lc['end'] = lc['end'].astype(int)
    for i, row in lc.iterrows():  ## fix so that start = leftmost coord (noticed that always start < end in the assemblies gff
        if row.start > row.end:
            lc.loc[i, 'start'] = row.end
            lc.loc[i, 'end'] = row.start

    lc.sort_values(['nodeID', 'start'], inplace=True)
    idx = 0
    for sampleID, group in lc.groupby("SampleID"):
        sample_idx += 1
        node_idx = 0

        node_end_coord = 0
        for nodeID, node in group.groupby('nodeID', sort=False):
            node_idx += 1
            node.sort_values('start', inplace=True)
            node_start_coord = node.start.head(1).item()  ## leftmost coord on node
            # node_end_coord = node.tail(1).end.astype(int).item()

            for iii, gene in node.iterrows():
                assert gene.start < gene.end
                generelstart = gene.start - node_start_coord + node_end_coord + 1 ### node_end_coord is 0 for one node samples
                generelend = gene.end - node_start_coord + node_end_coord + 1
                if mode == 'invivo':
                    v1  = "Mouse"
                    v2 = "Cage"
                    v1_v = gene.Mouse
                    v2_v = gene.Cage
                else:
                    v1  = "Mix"
                    v2 = "Replicate"
                    v1_v = gene.Mix
                    v2_v = gene.Replicate

                gene_coords = pd.concat([gene_coords, pd.DataFrame({'SampleID': [sampleID], v1:[v1_v], v2:[v2_v],'SampleType':[gene.SampleType], 'Day':[gene.Day],
                                                                    "nodeID": [nodeID], "N50":[gene.N50],"nodeLength":[gene.length], "OriScore":[gene.OriScore],
                                                                    'geneid':[gene.GeneID], "geneseq":[gene.DNA],
                                                                    'generelstart': [generelstart], 'genestart':[gene.start], 'geneend':[gene.end],
                                                                    'generelend': [generelend],
                                                                    'genestrand': [gene.strand],
                                                                    'genelength': [gene.GeneLength],
                                                                    'yval': [idx], 'nodeN': [node_idx],
                                                                    'RefBestHit': [
                                                                        gene.RefBestHit]})])
            node_end_coord = gene_coords.loc[gene_coords.nodeID == nodeID, 'generelend'].max()
            idx += 1
    gene_coords["GeneIDwithBestHit"] = gene_coords.geneid  + gene_coords['RefBestHit']
    dm = pd.pivot(gene_coords, index='nodeID', columns="GeneIDwithBestHit",values='genelength').fillna(0)

    from skbio.stats.distance import DistanceMatrix
    from sklearn.preprocessing import StandardScaler
    data_scaler = StandardScaler()
    #scaled_data = data_scaler.fit_transform(dm)
    from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
    #dm = DistanceMatrix(scipy.spatial.distance.pdist(bl)).data
    average_clustering = linkage(dm.values, method="average", metric="euclidean")
    dendro_idx = dendrogram(average_clustering, orientation='left')['leaves']
    ordered_samples = [dm.index.values[i] for i in dendro_idx]
    dendro_order = {k:v for k,v in zip(ordered_samples,range(len(ordered_samples)))}
    gene_coords['yval'] = [0 for x in gene_coords.index]
    gene_coords['yval'] = gene_coords['nodeID'].map(dendro_order)
    gene_coords['yval'] = gene_coords['yval'].apply(lambda x: 2*x)
    gene_coords.loc[gene_coords.SampleID == 'RefSeq', 'yval'] = -10

    return gene_coords

def check_gene_order(lc:pd.DataFrame, reference=["NC_004663.1","NC_004703.1"]):
    """
    :param lc: a dataframe that stores gene coords and strand
    :return:
    """
    #lc['libraryOrientation'] = ["+" for x in lc.index]
    lc['OriWithRef'] = [None for x in lc.index]
    lc['OriScore'] = [None for x  in lc.index]
    lc['start'] = lc['start'].astype(int)
    lc['end'] = lc['end'].astype(int)
    ref = lc.loc[lc.nodeID.isin(reference)]
    total = ref.shape[0] ## each entry is a reference gene..

    for i, x in zip(track(range(len(lc.groupby(['SampleID','nodeID']).groups))),lc.groupby(["SampleID",'nodeID'])):
        name, group = x
        score = 0
        this_total = 0
        for i, row in group.iterrows():
            #print(row.nodeID, row.GeneID, row.strand, row.start, row.end)
            assert row.start < row.end
            if not ref.loc[ref.GeneID == row.RefBestHit].empty:
                if row.strand == ref.loc[ref.GeneID == row.RefBestHit,'strand'].item():
                    score += 1
                    this_total += 1
                    lc.loc[(lc.SampleID == name[0]) &( lc.nodeID == name[1]) & (lc.GeneID == row.GeneID), 'OriWithRef'] = '+' ### this should be determined after final score remove!!
                else:
                    lc.loc[(lc.SampleID == name[0]) & (lc.nodeID == name[1]) & (lc.GeneID == row.GeneID), 'OriWithRef'] = '-'
                    this_total += 1
            else:
                this_total += 1
        #print(lc[lc.SampleID == name]['OriWithRef'].value_counts(dropna=False))

        lc.loc[(lc.SampleID == name[0]) &( lc.nodeID == name[1]), "OriScore"] = score /this_total
        #print(f'{name}| {score/this_total:.2%} ({score}/{this_total} | {this_total}/{total} of reference genes detected.)')

        if score / this_total < 0.5:
            #lc.loc[(lc.SampleID == name[0]) & (lc.nodeID == name[1]),'libraryOrientation'] = "-"
            lc.loc[(lc.SampleID == name[0]) & (lc.nodeID == name[1]),'DNA'] = lc.loc[(lc.SampleID == name[0]) & (lc.nodeID == name[1])]['DNA'].apply(lambda x: str(Seq(x).reverse_complement()))
            lc.loc[(lc.SampleID == name[0]) & (lc.nodeID == name[1]),'strand'] = lc.loc[(lc.SampleID == name[0]) & (lc.nodeID == name[1]),'strand'].apply(lambda x: "-" if x == "+" else "+")


    return lc
