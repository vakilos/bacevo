import itertools
import os
import re
import sys

import pandas as pd
import rich.console
from rich import print
from ai_fun import *
from bacevo_fun import *

from Bio import AlignIO
tprint("\n\nGeneAligner (II)")
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
console = rich.console.Console(color_system='windows')

parser = argparse.ArgumentParser(prog='ai_gene_alignment.py',description='Compare gene order of homologous genomic regions')
parser.add_argument('-d', dest='basedir',help='Project directory. Default = current working directory', default=os.getcwd())
parser.add_argument('-o', dest='outDir',help='Directory to export files.. Default = ./gene_order', default=os.getcwd()+'/gene_order')
parser.add_argument('-r', dest='ref', help='Path to Reference fasta', default='/home/christos/bacevo/reference/Btheta-VPI-5482/GCF_000011065.1/GCF_000011065.1_ASM1106v1_genomic.fna')
parser.add_argument('-a', dest='anno', help='Path to AnnotationTable.txt', default="/home/christos/bacevo/reference/Btheta-VPI-5482/GCF_000011065.1/AnnotationTable.txt")
parser.add_argument('-s', dest='sampleDir', help='Directory with assemblies for each sample.')
parser.add_argument('-n', dest='rename', help='Renames fasta files for sibeliaz, creates _sib.fa', default=False, action='store_true')
parser.add_argument('-g', dest='gene', help='Gene group name', default='', type=str, required=True)
parser.add_argument('-m', dest='mode',help='Invivo or invitro dataset. Default = invivo.', choices=['invivo','invitro'], default='invivo')
args = parser.parse_args()

os.chdir(args.basedir)
gene = args.gene
df = pd.read_csv(f"gene_order/{gene}_relative_coords.csv")

if args.rename:
    os.chdir(args.sampleDir)
    for _,folder in zip(track(range(len(glob.glob("*"))),description="## Preparing files..."),glob.glob("*")):
        if not folder.startswith('NC'):
            rename_fasta_bins(f'{folder}/{folder}.fna', folder)
    os.chdir(args.basedir)

#print(df.loc[(df.SampleID.isin(['S50963', 'RefSeq'])) & (df.oldName == 'BT2260'),'geneseq'].values)
#print(sys.exit())
ano = pd.read_csv(args.anno,sep='\t')
ano['old_locus_tag'] = ano['old_locus_tag'].fillna('')
onenodesamples = []
for name, group in df.groupby("SampleID"):
    if len(group.nodeID.unique()) == 1:
        onenodesamples += [name]
os.chdir(args.outDir)


def get_genes_within_range(ano: pd.DataFrame, limits: list):
    boundaries = ''
    limits = sorted(limits)
    ano[['start', 'end']] = ano[['start', 'end']].astype(int)
    ingene = ''

    for i, row in ano.drop_duplicates('locus_tag').iterrows():
        for x in range(*limits):
            if row.start < x and row.end > x:
                if len(row['old_locus_tag']) > 0:
                    ingene +=row['old_locus_tag']
                else:
                    ingene += row['locus_tag']
                boundaries = f"{limits[0] - row.start + 1}-{limits[1] - row.start + 1}"
                break
    return ingene, boundaries
def compare_homologue_sequences(df:pd.DataFrame, ref_fna = '', anotable=''):

    rf = df.loc[df.SampleID == 'RefSeq']
    bl = pd.DataFrame()
    ## use only genes present in reference (otherwirse comparison is impossible)
    """
    All homologous sequences fed to blastn are aligned (start --> end) despite different orientation on contigs!!!
    
    """
    for i, x in zip(track(range(len(df.groupby('RefBestHit').groups)), description='## Comparing homologs...'),df.loc[df.RefBestHit.isin(rf.geneid.values)].groupby('RefBestHit')):
        name, group = x
        print(f'\n\n{name}')
        #hom_seqs = []
        ref_seq = [SeqRecord(description=name,id='RefSeq', seq=Seq(rf.loc[rf.geneid == name,'geneseq'].item()))]
        SeqIO.write(ref_seq, 'sub.fa', 'fasta')

        for i, gene in group.iterrows():
            if gene.SampleID != 'RefSeq':
                if gene.genestrand == rf.loc[rf.geneid == name, 'genestrand'].item():
                    flip = False
                else:
                    flip = True
                if gene.OriScore >= 0.5: ## query and sbjct have the same orientation
                    hom_seqs = [SeqRecord(description=str(gene.oldName),id=f'{gene.SampleID}-{gene.VariantName}', seq=Seq(gene.geneseq))]
                else:
                    hom_seqs = [SeqRecord(description=str(gene.oldName), id=f'{gene.SampleID}-{gene.VariantName}', seq=Seq(gene.geneseq).reverse_complement())]
                SeqIO.write(hom_seqs, 'query.fa', 'fasta')
                #output = NcbiblastnCommandline(query='query.fa',subject='sub.fa' ,outfmt=1, out='blastn.xml')()[0]
                os.system(str(NcbiblastnCommandline(query='query.fa',subject='sub.fa', outfmt=7, out='oneblastn.csv')))
                out = subprocess.getoutput('grep -v "#" oneblastn.csv').split('\t')

                thisbl = pd.DataFrame([out+[name, gene.genelength, gene.genestrand, flip, gene.OriScore, gene.geneseq, rf.loc[rf.geneid == name, 'geneseq'].item()]],
                                      columns=["query",'subject','identity','alignmentlength',
                              'mismatches','gapopens','qstart','qend','sstart','send','evalue',
                                               'bitscore','genename','genelength','orientation', 'flip', 'oriscore','qseq','sseq'])
                thisbl['startdif'] = thisbl.qstart.astype(int) - thisbl.sstart.astype(int)
                thisbl['enddif'] = thisbl.qend.astype(int) - thisbl.send.astype(int)
                thisbl['SampleID'] = [gene.SampleID for x in thisbl.index]
                thisbl['nodeID'] = [gene.nodeID for x in thisbl.index]
                thisbl['RefBestHit'] = [gene.RefBestHit for x in thisbl.index]
                thisbl['VariantName'] = [gene.VariantName for x in thisbl.index]
                #thisbl['VariantFreq'] = [gene.VariantFreq for x in thisbl.index]
                thisbl['generelstart'] =[gene.generelstart for x in thisbl.index]
                thisbl['generelend'] = [gene.generelend for x in thisbl.index]
                thisbl['oldName'] = [gene.oldName for x in thisbl.index]
                thisbl['genestrand'] = [gene.genestrand for x in thisbl.index]
                bl = pd.concat([bl,thisbl])

    bl.reset_index(inplace=True)
    bl = bl.drop(columns='index')

    bl['fragseq'] = ['' for x in bl.index]
    bl['geneplace'] = ['' for x in bl.index]
    bl['diftype'] = ['' for x in bl.index]
    bl['description'] = ['' for x in bl.index]
    bl['fragrelstart'] = [None for x in bl.index]
    bl['fragrelend'] = [None for x in bl.index]

    #print(bl[['query','genename','genelength','qstart','qend','sstart','send','startdif','enddif']])
    for _,x  in zip(track(range(bl.shape[0]), description='## Labelling genomic rearrangements..'),bl.iterrows()):
        """
        generelcoords are by definition start < end and genestrand gives the orientation, so i need to add insertions
        or deletions accordingly to each end based on the genestrand. Cause insertstart will correspond to generelstart only if
        genestrand == + !!!!!
        """
        i,row = x
        if row.startdif == row.enddif:
            if row.startdif > 0:
                #print(f'There is an insertion at the start of the gene!')
                bl.loc[i, 'fragseq'] = row.qseq[:row.startdif]
                bl.loc[i, 'geneplace'] = 'start'
                bl.loc[i, 'diftype'] = 'INS'
                bl.loc[i, 'description'] += f'{row.startdif} bp insertion at the gene start.'
                if row.genestrand == "+":
                    bl.loc[i, ['fragrelstart','fragrelend']] = [row.generelstart, row.generelstart + row.startdif-1]
                else:
                    bl.loc[i, ['fragrelstart', 'fragrelend']] = [row.generelend, row.generelend - row.startdif]
                assert len(bl.loc[i, 'fragseq']) ==  row.startdif
            elif row.startdif < 0:
                #print(f'There is an deletion at the start of the gene!')
                bl.loc[i, 'fragseq'] = row.sseq[int(row.qstart)-1:int(row.sstart)-1]
                bl.loc[i, 'geneplace'] = 'start'
                bl.loc[i, 'description'] += f'{row.startdif} bp deletion at the gene start.'
                bl.loc[i, "diftype" ] = 'DEL'
                if row.genestrand == "+":
                    bl.loc[i, ['fragrelstart','fragrelend']] = [row.generelstart, row.generelstart + row.startdif-1]
                else:
                    bl.loc[i, ['fragrelstart', 'fragrelend']] = [row.generelend - row.startdif, row.generelend - 1]

                assert len(bl.loc[i, 'fragseq']) == -row.startdif
            else:
                #print(f'Gene is identical to RefSeq!')
                if row.mismatches == "0" and row.gapopens == "0":
                    bl.loc[i, 'description'] += f'Gene identical to RefSeq'
                else:
                    bl.loc[i, 'description'] += f'Mutation to RefSeq'
        elif row.startdif == 0 and row.enddif > 0:
            print(f'There is an insertion in the end of the gene!!')
            bl.loc[i, 'fragseq'] = row.qseq[int(row.send): int(row.qend)]
            bl.loc[i, 'geneplace'] = 'end'
            bl.loc[i, 'diftype'] = 'INS'
            bl.loc[i, 'description'] += f'{row.enddif} bp insertion at the gene end'
            if row.genestrand == "+":
                bl.loc[i, ['fragrelstart','fragrelend']] = [row.generelend+1, row.generelend + row.enddif]
            else:
                bl.loc[i, ['fragrelstart', 'fragrelend']] = [row.generelstart-1,row.generelstart - row.enddif]
            print(f"{len(bl.loc[i, 'fragseq'])} - {row.enddif}")
            assert len(bl.loc[i, 'fragseq']) == row.enddif
            #print(len(bl.loc[i, 'insertend']), row.enddif)

        elif row.startdif == 0 and row.enddif < 0:
            print(f'There is an deletion in the end of the gene!!')
            bl.loc[i, 'fragseq'] = row.sseq[int(row.qend): int(row.send)]
            bl.loc[i, 'geneplace'] = 'end'
            bl.loc[i, 'diftype'] = 'DEL'
            bl.loc[i, 'description'] += f'{row.enddif} bp deletion at the gene end'
            if row.genestrand == '+':
                bl.loc[i, ['fragrelstart','fragrelend']] = [row.generelend + 1 , row.generelend - row.enddif]
            else:
                bl.loc[i, ['fragrelstart', 'fragrelend']] = [row.generelstart + 1, row.generelstart - row.enddif]
            assert len(bl.loc[i, 'fragseq']) == -row.enddif
        else:
            print(row[['query','genename','genelength','qstart','qend','sstart','send','startdif','enddif']])
            logging.warning('Need to add this mutational scenario!!!!')

    os.system('rm oneblastn.csv')

    bl['fragName'] = ['' for x in bl.index]
    bl['fragorigin'] = ['' for x in bl.index]
    for i, row in bl.loc[bl.description != 'Gene identical to RefSeq'].iterrows():
        fragname = ''
        origin = ''
        if row.diftype == "INS":
            for ii, deletion in bl.loc[(bl.diftype == 'DEL') & (bl.SampleID == row.SampleID)].drop_duplicates('fragseq').iterrows():
                if row.fragseq in deletion.fragseq or deletion.fragseq in row.fragseq:
                    if len(row.fragseq) == len(deletion.fragseq):
                        fstart, fend = sorted([int(x) for x in [deletion.fragrelstart,deletion.fragrelend]])
                        fragname += deletion.oldName + f' ({rf.loc[rf.RefBestHit == deletion.RefBestHit,"generelstart"].item()-fstart}-{rf.loc[rf.RefBestHit == deletion.RefBestHit,"generelstart"].item()-fend}\n{len(row.fragseq)} bp)'
                        origin += deletion.oldName
                    # else:
                    #     fragname += deletion.oldName + f' (?-?\n{len(row.fragseq)} bp)'
                    #     origin += deletion.oldName

            if len(fragname) == 0: ## if there is no exact much for the insertion (in the deletions) blastn against the genome)
                #ref = fasta_to_df(ref_fna)
                hom_seqs = [SeqRecord(description=row.RefBestHit, id=f'{row.SampleID}-{row.VariantName}', seq=Seq(row.fragseq))]
                SeqIO.write(hom_seqs, 'fragquery.fa', 'fasta')
                os.system(str(NcbiblastnCommandline(query='fragquery.fa',subject=ref_fna, outfmt=7, out='oneblastn.csv')))
                #os.system(f"blastn -query fragquery.fa -subject {ref_fna} -outfmt 7 -out oneblastn.csv ")
                #out = subprocess.getoutput('grep -v "#" oneblastn.csv').split('\t')
                #os.system(str(NcbiblastnCommandline(query='fragquery.fa',subject=ref_fna ,outfmt=7, out='oneblastn.csv')))
                out = [x.split('\t') for x in subprocess.getoutput('grep -v "#" oneblastn.csv').splitlines()]
                print(out)
                if len(out) > 0:
                    hit = out[0] ## get the first hit

                    print(f'query length: {len(row.fragseq)}, Alignment length: {hit[3]}, mismatch: {hit[4]}, hit_coords: {hit[8]}- {hit[9]}')
                    geneoriginname, boundaries = get_genes_within_range(ano, [int(hit[8]),int(hit[9])])
                    limits = sorted([int(hit[8]),int(hit[9])])
                    fragname = f'{geneoriginname} ({boundaries}\n{limits[1]-limits[0]} bp)'
                    origin += geneoriginname
                    print(fragname)
                    print('\n\n')
                else:
                    print(f'{row.SampleID} - {row.VariantName} | {row.genename}:{row.description}| No hit found')
                    origin += 'nohit'
                    fragname += f'({len(row.fragseq)} bp)'
                    print('\n\n')
                os.system('rm fragquery.fa oneblastn.csv')

            bl.loc[i, 'fragName'] = fragname
            bl.loc[i, 'fragorigin'] = origin
        elif row.diftype == "DEL":
            for ii, insertion in bl.loc[(bl.diftype == "INS") & (bl.SampleID == row.SampleID)].drop_duplicates('fragseq').iterrows():
                if row.fragseq in insertion.fragseq or insertion.fragseq in row.fragseq:
                    fragname += ""
                    origin = ""
                    break
            else:
                #fragname += row.oldName
                origin = ""


            bl.loc[i, 'fragName'] = fragname + f' ({len(row.fragseq)} bp)'
            bl.loc[i, 'fragorigin'] = origin

    if not bl.empty:
        out = df.merge(bl, on=['SampleID','nodeID','RefBestHit','oldName','VariantName', 'generelstart','generelend', 'genestrand','genelength'], how='outer')
    else:
        out = pd.DataFrame()
    return out


def locate_motifs_in_genomic_region_of_interest(df:pd.DataFrame, samplesDir, motifs:dict):
    """
    This is relevant only for one node samples
    :param df:
    :param samplesDir:
    :return:
    """
    motif_df = pd.DataFrame()
    motifs_height = {x:v for x,v in zip(motifs,range(len(motifs)*2))}

    oneNodeSamples = []
    refname = 'NC_004663.1'
    for name,group in df.groupby("SampleID"):## remove samples with targets in more than one contig
        if len(group.nodeID.unique()) == 1:
            oneNodeSamples += [name]
    for name, group in df.loc[df.SampleID.isin(oneNodeSamples)].groupby(["SampleID",'nodeID']):
        samplename,nodename = name

        region_start = group.genestart.astype(int).min()-1
        region_end = group.geneend.astype(int).max()+2
        if samplename != "RefSeq":
            fna = fasta_to_df(f'{samplesDir}/{samplename[1:]}/{samplename[1:]}_sib.fna')
        else:
            fna = fasta_to_df(f'{samplesDir}/{refname}/{refname}.fna')
        chrom = Seq(str(fna.loc[fna.Chrom == nodename]['seq'].item()))

        if group.OriScore.astype(float).values[0] < 0.5:
            chrom = chrom.reverse_complement()
            region = chrom[region_start:region_end]
            reg = str(region)

            for i ,row in group.iterrows():
                if row.genestrand  == "-":

                    #print(f'{samplename}:{nodename}-{row.oldName}|{row.OriScore}| {row.genestrand}\n{chrom[int(row.genestart) : int(row.geneend)+1].reverse_complement()}\n{row.geneseq}')
                    print(f'{samplename}:{nodename}-{row.oldName}|{row.OriScore}| {row.genestrand}\n{region[int(row.generelstart)-1:int(row.generelend)].reverse_complement()}\n--\n{Seq(row.geneseq).reverse_complement()}')

                else:
                    #print(f'{samplename}:{nodename}-{row.oldName}|{row.OriScore}| {row.genestrand}\n{chrom[int(row.genestart) : int(row.geneend)+1]}\n{row.geneseq}')
                    print(f'{samplename}:{nodename}-{row.oldName}|{row.OriScore}| {row.genestrand}\n{region[int(row.generelstart)-1:int(row.generelend)]}\n{Seq(row.geneseq).reverse_complement()}')
        else:
            region = chrom[region_start:region_end]
            reg = str(region)
            for i ,row in group.iterrows():
                if row.genestrand  == "-":

                    #print(f'{samplename}:{nodename}-{row.oldName}|{row.OriScore}| {row.genestrand}\n{chrom[int(row.genestart) : int(row.geneend)+1].reverse_complement()}\n{row.geneseq}')
                    print(f'{samplename}:{nodename}-{row.oldName}|{row.OriScore}| {row.genestrand}\n{region[int(row.generelstart)-1:int(row.generelend)].reverse_complement()}\n{row.geneseq}')
                else:
                    #print(f'{samplename}:{nodename}-{row.oldName}|{row.OriScore}| {row.genestrand}\n{chrom[int(row.genestart) : int(row.geneend)+1]}\n{row.geneseq}')
                    print(f'{samplename}:{nodename}-{row.oldName}|{row.OriScore}| {row.genestrand}\n{region[int(row.generelstart)-1:int(row.generelend)]}\n{row.geneseq}')
        print(f"\n\n{nodename} - {group.VariantName.values[0]} | {region_start} - {region_end} | ({region_end-region_start} bp)")

        for motif in motifs:

            motif_seq = motifs[motif]
            if len([(x.start(),x.end()) for x in re.finditer(motif_seq, str(region))]) > 0:
                start, end = zip(*[(x.start(),x.end()) for x in re.finditer(motif_seq, str(region))])
                for s,e in zip(start,end):
                    motif_df = pd.concat([motif_df,pd.DataFrame({'SampleID': [samplename], 'VariantName': [group.VariantName.values[0]],'nodeID':[nodename],
                                  'motifrelstart':[s], 'motifrelend':[e], 'motif':[motif+"(+)"],'motifseq':[motif_seq],'motifori':[1]})])
            else:
                start,end =['','']

            motif_seq_r = motif_seq[::-1]
            if len([(x.start(),x.end()) for x in re.finditer(motif_seq_r, str(region))]) > 0:
                start_r,end_r = zip(*[(x.start(),x.end()) for x in re.finditer(motif_seq_r, str(region))])
                for s,e in zip(start_r,end_r):
                    motif_df = pd.concat([motif_df,pd.DataFrame({'SampleID': [samplename], 'VariantName': [group.VariantName.values[0]],'nodeID':[nodename],
                                  'motifrelstart':[s], 'motifrelend':[e], 'motif':[motif+"(r)"],'motifseq':[motif_seq_r],'motifori':[0]})])
            else:
                start_r,end_r =['','']
            motif_seq_rc =  str(Seq(motif_seq).reverse_complement())
            if len([(x.start(),x.end()) for x in re.finditer(motif_seq_rc, str(region))]) > 0:
                start_rc,end_rc = zip(*[(x.start(),x.end()) for x in re.finditer(motif_seq_rc, str(region))])
                for s,e in zip(start_rc,end_rc):
                    motif_df = pd.concat([motif_df,pd.DataFrame({'SampleID': [samplename], 'VariantName': [group.VariantName.values[0]],
                                                                 'nodeID':[nodename],
                                  'motifrelstart':[s], 'motifrelend':[e], 'motif':[motif+"(rc)"],'motifseq':[motif_seq_rc],'motifori':[0]})])
            else:
                start_rc,end_rc =['','']

            print(f"Motif: {motif}-{motif_seq}(+)  found: {Seq(str(region)).count(motif_seq)} [{start} - {end}]") ## find returns index
            print(f"Motif: {motif}-{motif_seq_r}(r)  found: {Seq(str(region)).count(motif_seq_r)} [{start_r} - {end_r}]")
            print(f'Motif: {motif}-{motif_seq_rc}(rc) found: {Seq(str(region)).count(motif_seq_rc)} [{start_rc} - {end_rc}]')

    motif_df.reset_index(inplace=True)
    motif_df = motif_df.drop(columns='index')

    motif_df['RefBestHit'] = ['' for x in motif_df.index] ## will store data for genes associated with motifs
    motif_df['overlap'] = ['Intergenic' for x in motif_df.index]## will store data for type of overlap with associated gene

    for i, motif in motif_df.iterrows():
        print(f'\n\n{motif.motif}: {motif.motifrelstart} - {motif.motifrelend}')
        for ii, gene in df.loc[df.nodeID == motif.nodeID].drop_duplicates(['RefBestHit','generelstart']).iterrows():

            print(f'{gene.generelstart} - {gene.generelend}')
            if gene.generelstart >= motif.motifrelstart and gene.generelend <= motif.motifrelend:
                motif_df.loc[i, 'overlap'] = 'Gene within motif! lol!'
                motif_df.loc[i, 'RefBestHit'] = gene.RefBestHit



            elif gene.generelstart < motif.motifrelstart and gene.generelend > motif.motifrelend:
                motif_df.loc[i, 'overlap'] = 'Motif within gene'
                motif_df.loc[i, 'RefBestHit'] = gene.RefBestHit



            elif motif.motifrelstart < gene.generelstart and motif.motifrelend >= gene.generelstart and motif.motifrelend < gene.generelend:
                motif_df.loc[i, 'overlap'] = 'Motif partially at gene start'
                motif_df.loc[i, 'RefBestHit'] = gene.RefBestHit
                print(motif_df.loc[i, 'overlap'])

            elif motif.motifrelend > gene.generelend and motif.motifrelstart > gene.generelstart and motif.motifrelstart <= gene.generelend:
                motif_df.loc[i ,'overlap'] = "Motif partially at gene end"
                motif_df.loc[i, 'RefBestHit'] = gene.RefBestHit
        print(motif_df.loc[i, 'overlap'])

    out = df.merge(motif_df, on=['SampleID','nodeID','RefBestHit','VariantName'], how='outer')
    out.loc[out['overlap'] == 'Intergenic', ['generelstart','generelend']] = out.loc[out['overlap'] == 'Intergenic', ['motifrelstart','motifrelend']]

    # ## correct motifori (1 if pointing to gene's 5'-end or 0 if reversed to 3'-end - in case gene orientation (on the genome) is - reverted the relative orientations of the motifs
    # for i, row in out.iterrows():
    #     if row.genestrand == '-':
    #         if row.motifori == 1:
    #             out.loc[i, 'motifori'] = 0
    #         if row.motifori == 0:
    #             out.loc[i, 'motifori'] = 1

    return out




motifs = {'olga':'CTTTGCAGA','michi1':"AAGATAG", 'michi2':'GCGATGG',"BT1042_rs":'CTCCCAA'}
mof = locate_motifs_in_genomic_region_of_interest(df, args.sampleDir, motifs=motifs)
# print(mof.sort_values('motif').loc[mof.SampleID == 'RefSeq',['motif', 'oldName','overlap','generelstart','generelend','motifrelstart','motifrelend', 'motifori']])
# mof.to_csv(f'{gene}_motifs_rel_coords.csv', index=False)
df = df.loc[df.SampleID.isin(onenodesamples)]
bl = compare_homologue_sequences(df, ref_fna=args.ref)

if args.mode == 'invivo':
    if not bl.empty:
        fin = bl.merge(mof, on = ['SampleID', 'Mouse', 'Cage', 'SampleType', 'Day', 'nodeID', 'N50',
               'nodeLength', 'OriScore', 'geneid', 'geneseq', 'generelstart',
               'genestart', 'geneend', 'generelend', 'genestrand', 'genelength',
               'yval', 'nodeN', 'RefBestHit', 'GeneIDwithBestHit', 'string',
           'VariantName', 'VariantFreq', 'oldName', 'info'], how='outer')
        fin.to_csv(f"{gene}_fragments_motifs_rel_coords.csv", index=False)

else:
    if not bl.empty:
        fin = bl.merge(mof, on = ['SampleID', 'Mix', 'Replicate', 'SampleType', 'Day', 'nodeID', 'N50',
               'nodeLength', 'OriScore', 'geneid', 'geneseq', 'generelstart',
               'genestart', 'geneend', 'generelend', 'genestrand', 'genelength',
               'yval', 'nodeN', 'RefBestHit', 'GeneIDwithBestHit', 'string',
           'VariantName', 'VariantFreq', 'oldName', 'info'], how='outer')
        fin.to_csv(f"{gene}_fragments_motifs_rel_coords.csv", index=False)

if not bl.empty:
    for name, group in bl.loc[bl.SampleID.isin(onenodesamples)].sort_values('genename').groupby('query'):
        print(f"\n\n\n\_______________{name}_______________\n")
        ### print mutational events
        for i, row in group.loc[group.description != 'Gene identical to RefSeq'].sort_values("genename").iterrows():
            print(f'{row.genename} - {row.oldName} | {row.description} | {row.fragseq} | mismatches: {row.mismatches}')
