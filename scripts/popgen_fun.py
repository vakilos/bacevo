import itertools
import logging
import os.path
import sys

import numpy as np
import pandas as pd

from bacevo_fun import *
def find_codon_haplotypes_from_isolates(poly_df):
    """A FinalPolyTable with all samples isolates derived from the same population. (needs to be updated according to metagenome function"""
    console = Console()
    poly_df = poly_df.loc[(poly_df.Type == "SNP") & (poly_df.cds != "Nan") & (poly_df.SampleProfile == 'isolate')]
    N_isolates = poly_df.drop_duplicates('SampleID').shape[0]
    N_polys = poly_df.drop_duplicates("PolyID").shape[0]
    console.print(f'A population has occured from {poly_df.drop_duplicates("SampleID").shape[0]} isolates . ')
    console.print(f'Will explore codons for {poly_df.shape[0]} polys!')
    haplo_df = pd.DataFrame()
    for i, (name, group) in zip(track(range(len(poly_df.groupby(["SampleID","GeneName"]).groups)), description='[b]Locating polys in codons...[/b]'),
                                poly_df.groupby(["SampleID","GeneName"])): ## find which polys from the same samples belong to the same codon
        for i,poly in group.iterrows():
            pos_idx = poly.Pos - 1
            gene_start = poly.start
            gene_length = len(poly.cds)
            #print(poly.GeneName)
            #print(poly.Pos, poly.start)
            #print(poly.strand)
            idx_on_gene = int(poly.Pos - float(poly.start)-1) ## genbank file uses as locators indexes you need to use to extract this sequence from the reference fasta. (start = location.start-1 end = end)
            try:
                if poly.strand == "1.0":
                    poly.Ref == poly.cds[idx_on_gene]
                else:
                    poly.Ref == str(Seq(poly.cds).reverse_complement()[idx_on_gene])
            except:
                print(f"WARNING:: REF from VCF does not match genbank annotation for {poly.GeneName}!!")

            idx_on_codon = idx_on_gene % 3  ## this is 0-index --> 0 = 1st , 1 = 2nd, 2 = 3rd
            idx_codon_on_cds = int(np.floor(idx_on_gene / 3)) ## the codon index on the cds
            if poly.strand == "1.0":
                ref_codon = chunks_h(poly.cds,3)[idx_codon_on_cds]
            else:
                ref_codon = chunks_h(str(Seq(poly.cds).reverse_complement()),3)[idx_codon_on_cds]
            #codon_start = 0
            #codon_end = 0
            haplo_df = pd.concat([haplo_df,pd.DataFrame({"SampleID":[name[0]],"PolyID":[poly.PolyID], "Alt":[poly.Alt], "SampleID":[poly.SampleID], "strand":[poly.strand],"cds":[poly.cds],
                                                     "GeneName":[poly.GeneName], "GeneLength":[gene_length], "GeneStart":[gene_start],
                                                     "RefCodon":[ref_codon],"Idx_on_codon":[idx_on_codon],"Idx_on_cds":[idx_codon_on_cds]})],ignore_index=True)

    haplo2_df = pd.DataFrame()
    for i, (name, group) in zip(track(range(len(haplo_df.groupby(["SampleID","GeneName","Idx_on_cds"]))), description='[b]Reconstructing sample codons...[/b]'),
                                haplo_df.groupby(["SampleID","GeneName","Idx_on_cds"])):## create SampleCodon - for multiple polys for the same codon in the isolate sample
        sample_codon = list(group["RefCodon"].values[0])
        #print(f"RefCodon={''.join(sample_codon)}, strand={group.strand.values[0]}")
        #print(f"Detected {group.shape[0]} polys in the codon.")
        for i, poly in group.iterrows():
            sample_codon[poly["Idx_on_codon"]] = poly.Alt
        sample_codon = "".join(sample_codon)
        #print(f'SampleCodon={"".join(sample_codon)}\n')
        haplo2_df = pd.concat([haplo2_df,pd.DataFrame({'SampleCodon':[sample_codon],'RefCodon':[group.RefCodon.values[0]],
                                                   'Index_on_cds':[name[2]],'GeneName':[name[1]], "GeneLength":[group.GeneLength.values[0]],"GeneStart":[group.GeneStart.values[0]]})], ignore_index=True)
        ## for each codon i need to count how many times different codons appear..
        ## counts in the SampleCodon DO NOT include Ref_codon..to estimate this ref_codon_counts = all_samples - counts for all alternative codons
    haplo2_df = haplo2_df.value_counts().to_frame().reset_index()
    return haplo2_df

def codon_haplotypes_from_bam_OLD(bam, sampleid, poly_df, mapping_quality=20, qc=20,codon_read_support=4):
    pd.set_option('display.max_rows', None)

    if not os.path.isfile(bam + '.bai'): os.system(f'samtools index -b {bam}')
    bamfile = pysam.AlignmentFile(bam, 'rb')
    haplotype_df = pd.DataFrame()
    breseq_only = []
    bam_only = []
    not_detected_genes = 0
    if not poly_df.empty:
        for indx, (name, group) in zip(
                track(range(len(poly_df.groupby("GeneName").groups)),
                      description=f'[b]Parsing codon haplotypes ({sampleid})...[/b]'),
                poly_df.groupby('GeneName')):
            cds = group.cds.values[0]
            gene_length = len(cds)

            if gene_length % 3 == 0:  ### work only with coding genes.
                codon_df = pd.DataFrame()
                gene_start = group.start.values[0]
                strand = group.strand.values[0]
                if strand == "1.0":
                    DNA = cds
                else:
                    DNA = str(Seq(cds).reverse_complement())

                # detected_polys = []
                ### split CDS to codons and their coordinates
                for idx_codon_on_cds, codon in enumerate(chunks_h(DNA, 3)):
                    codon_start = gene_start + (idx_codon_on_cds * 3)
                    codon_end = codon_start + 2

                    # assert DNA[0] == ref.loc[ref.Chrom == group.Chrom.values[0] + ".1", "seq"].item()[gene_start]
                    if not group.loc[group.Pos.between(codon_start + 1,
                                                       codon_end + 1)].empty:  ### select only codons with detected mutations
                        # print(f'\n\t\t\t\t\t {sampleid}:: {name} ({strand}) - [{gene_start} - {group.end.values[0]}] ')
                        # print(f'Codon ({idx_codon_on_cds}): {codon_start + 1}, {codon_end + 1}')
                        # print("Polys associated with codon:")
                        # print(group.loc[group.Pos.between(codon_start + 1,
                        #                                codon_end + 1)][['PolyID','TotalCov','AltCov']])
                        # print('\n')
                        # both_reads = False
                        # print(ref.loc[ref.Chrom == group.Chrom.values[0] + ".1", "seq"].item()[codon_start:codon_end + 1])
                        for read in bamfile.fetch(group.Chrom.values[0], codon_start, codon_end):
                            if read.mapping_quality >= mapping_quality and read.is_secondary is False and read.is_duplicate is False and read.is_qcfail is False:  # and read.is_proper_pair is True:
                                region = [x for x in read.get_aligned_pairs(matches_only=True, with_seq=True) if
                                          codon_start <= x[1] <= codon_end]
                                codon = {codon_start: "-", codon_start + 1: "-", codon_start + 2: "-"}

                                for position in region:
                                    if read.query_qualities[position[0]] >= qc:
                                        if position[2].islower():  ## mismatch with ref

                                            codon[position[1]] = read.query_sequence[position[0]]
                                            # detected_polys += [
                                            #   f'{group.Chrom.values[0]}:{position[1] + 1}:{DNA[position[1] - gene_start]}:{read.query_sequence[position[0]]}']
                                            # codon += read.query_sequence[position[0]]  ## extract alternative base from the read (position[0] is the read index)
                                        else:
                                            codon[position[1]] = position[2]
                                            # codon += position[2]  ## This is the reference base (third element in the tuple for a particular site

                                codon_df = pd.concat([codon_df, pd.DataFrame(
                                    {'SampleID': [sampleid], 'RefCodon': [chunks_h(DNA, 3)[idx_codon_on_cds]],
                                     'SampleCodon': [codon],
                                     'GeneName': [name],
                                     'Index_on_cds': [idx_codon_on_cds], "GeneLength": [gene_length],
                                     'GeneStart': [gene_start], 'GeneEnd':[group.end.values[0]], 'CodonStart': [codon_start]})], ignore_index=True)





                ## impute missing bases (not covered by read) with reference DNA base
                def impute_missing_with_ref(row):
                    samplecodon = list(row.SampleCodon.values())
                    for i, (x, y) in enumerate(zip(samplecodon, row.RefCodon)):
                        if x == '-': samplecodon[i] = y
                    return "".join(samplecodon)

                # print(codon_df)
                impute = True
                if not impute:
                    remove_incomplete = True
                if impute:
                    codon_df["SampleCodon"] = codon_df.apply(lambda x: impute_missing_with_ref(x), axis=1)
                else:
                    codon_df["SampleCodon"] = codon_df["SampleCodon"].apply(lambda x: ''.join(x.values()))
                    if remove_incomplete:
                        codon_df = codon_df.loc[~codon_df.SampleCodon.str.contains("-")]
                ## count different codons
                codon_df = codon_df.value_counts().to_frame().reset_index()  # .rename(columns={0:'Count'})
                #print(codon_df.sort_values('Index_on_cds'))
                # codon_df['PolyID'] = [poly.PolyID for x in codon_df.index]

                # detected_polys = pd.Series(detected_polys).value_counts().reset_index()
                # detected_polys = detected_polys.loc[detected_polys['count'] >= codon_read_support]['index'].tolist()
                ### maybe filter for Freq 0.01 as well?

                # breseq_only += [x for x in group.PolyID.values if x not in detected_polys]
                # bam_only += [x for x in detected_polys if x not in group.PolyID.values]
                not_detected = 0

                codon_df = codon_df.loc[
                    codon_df['count'].ge(codon_read_support)]  ## select codons support at least by x reads..
                codon_df['Freq'] = codon_df['count'].apply(lambda x: x / codon_df['count'].sum())

                for codoname, codon_group in codon_df.groupby('Index_on_cds'):
                    if codon_group.loc[codon_group.SampleCodon != codon_group.RefCodon].empty:
                        # print(f'{name} Could not retrieve mutation for this codon! (detected previously with breseq)')
                        not_detected = 1
                        # print(int(codon_group.GeneStart.values[0]) + (int(codoname) * 3) + 1);
                        # print(int(codon_group.GeneEnd.values[0]) - (int(codoname) *3) + 1)
                        # print(f'Codon_index_on_cds: {codoname}')
                        # print(f"{codon_group['SampleCodon'].values[0]}, {codon_group['RefCodon'].values[0]}");
                        # print(CDS(DNA).codons[codoname])
                        # print('\n\n\n\n\n');
                        # sys.exit()

                not_detected_genes += not_detected
                # print(f'Gene: {name}| {group.oldGeneName.values[0]} ({strand})')
                # if len(detected_polys) != len(group.PolyID.values):
                #     print(sorted(detected_polys))
                #     print(group.PolyID.sort_values().values)
                #     print(codon_df[['RefCodon','SampleCodon','Index_on_cds','CodonStart','GeneStart','count']].sort_values('CodonStart'))
                #     print('\n')

                haplotype_df = pd.concat([haplotype_df, codon_df], ignore_index=True)
    print(f'{not_detected_genes} have not been detected in the alignment! (despite called with breseq)')
    return haplotype_df

    # print(f"Polys only reported from breseq {len(breseq_only)} (out of {poly_df.loc[poly_df['SampleID'] == sampleid].shape[0]}), polys only detected in BAM {len(bam_only)}")
    #invalid = pd.concat([invalid, pd.DataFrame(
    #    {'SampleID': [sampleid for x in range(len(breseq_only) + len(bam_only))], 'PolyID': breseq_only + bam_only,
    #     'only': ['breseq' for x in breseq_only] + ['bam' for x in bam_only]})])


    # print(invalid.drop_duplicates("PolyID")['only'].value_counts())

    #invalid.to_csv('invalid.csv', index=False)

def codon_haplotypes_from_genomes_OLD(FastaFiles, poly_df):
    N_isolates = poly_df.drop_duplicates('SampleID').shape[0]
    N_polys = poly_df.drop_duplicates("PolyID").shape[0]
    console.print(f'Filtering out {poly_df.loc[poly_df["Coverage_mean"] < 20].SampleID.nunique()} genomes with < 20X average genome coverage...')
    remove_samples = list(poly_df.loc[poly_df["Coverage_mean"] < 20].SampleID.unique())
    cnt = 0
    for x in FastaFiles:
        if not os.path.dirname(os.path.basename(x)) in remove_samples:
            cnt += 1
            print(os.path.dirname(os.path.basename(x)))
    print(cnt)
    sys.exit()
    #FastaFiles = [x for x in FastaFiles if os.path.dirname(os.path.basename(x)) not in remove_samples]

    poly_df = poly_df.loc[poly_df['Coverage_mean'] >= 20]

    console.print(f'A population has occured from {poly_df.drop_duplicates("SampleID").shape[0]} isolates...')
    console.print(f'Will explore codons for {poly_df.shape[0]} polys!')
    posDF = pd.DataFrame()

    for i, (name, group) in zip(track(range(len(poly_df.drop_duplicates(['GeneName','PolyID']).groupby(["GeneName", "PolyID"]).groups)), description='[b]Locating polys in codons...[/b]'),
                                poly_df.drop_duplicates(["GeneName","PolyID"]).groupby(["GeneName","PolyID"])):


        gene_length = len(group.cds.values[0])
        if gene_length % 3 == 0:
            pos = int(group.Pos.values[0])
            genename = name[0]
            gene_start = group.start.values[0]
            strand = group.strand.values[0]
            if strand == "1.0":
                DNA = group.cds.values[0]
            else:
                DNA = str(Seq(group.cds.values[0]).reverse_complement())

            #pos_idx = pos - 1
            idx_on_gene = pos - gene_start - 1 # subtract 1 because gene start from genbank is 0-indexed
            idx_on_codon = idx_on_gene % 3
            # if genename == "BT_RS19990":
            #     print(f'{group.PolyID.values[0]} | Gene start = {gene_start} | Genelength={gene_length} | Index on gene{idx_on_gene}, ')
            codon_start = 0
            codon_end = 0
            idx_codon_on_cds = int(np.floor((idx_on_gene) / 3))
            if idx_on_codon == 0:
                codon_start = pos
            if idx_on_codon == 1:
                codon_start = pos - 1
            if idx_on_codon == 2:
                codon_start = pos - 2
            codon_end = codon_start + 2

            assert chunks_h(DNA,3)[idx_codon_on_cds] == ref.loc[ref.Chrom == group.Chrom.values[0]+'.1']['seq'].item()[codon_start-1:codon_end]

            #if genename == "BT_RS19990":
            #    print(f"{chunks_h(DNA,3)[idx_codon_on_cds]} , {genename} idx_on_cds={idx_codon_on_cds} idx_on_codon={idx_on_codon} | codon_start={codon_start} | {strand} | {DNA[:18]}, {ref.loc[ref.Chrom == group.Chrom.values[0]+'.1']['seq'].item()[codon_start-1:codon_end]}")
            #    print(ref.loc[ref.Chrom == group.Chrom.values[0]+".1"]['seq'].item()[gene_start: gene_start+gene_length])
            #    print(len(DNA),len(ref.loc[ref.Chrom == group.Chrom.values[0]+".1"]['seq'].item()[gene_start: gene_start+gene_length]))
            #    for x,y in zip(DNA, ref.loc[ref.Chrom == group.Chrom.values[0]+".1"]['seq'].item()[gene_start: gene_start+gene_length]):
            #        if x != y:
            #            print(f'Mismatch:{x} != {y}')
            #            sys.exit()
            posDF = pd.concat([posDF, pd.DataFrame({"Chrom":[group.Chrom.values[0]],'Pos':[pos], 'strand':[strand],
                                                    'GeneName':[genename], "GeneStart":[gene_start],
                                                    'GeneLength':[gene_length],'Index_on_cds':[idx_codon_on_cds], 'PolyID':[name[1]],
                                               "RefCodon":[chunks_h(DNA,3)[idx_codon_on_cds]],"codon_start":[codon_start],"codon_end":[codon_end],
                                                })])
    print(posDF.drop_duplicates(['GeneName','PolyID']).shape[0])
    print(f'{posDF.shape[0]} codons will be compared..')
    haplo = pd.DataFrame()
    ### this iterates over all samples in the directory - need to filter based on genome coverage which ones to include
    print(len(FastaFiles))
    for sample in FastaFiles:
        fasta = [x for x in FastaFiles if x == sample][0]
        fasta = fasta_to_df(fasta)
        for i, row in posDF.iterrows():

            codon = fasta.loc[fasta.Chrom == row.Chrom]['seq'].item()[row['codon_start']-1:row['codon_end']]
            haplo  = pd.concat([haplo, pd.DataFrame({"SampleID":[sample], "GeneName":[row.GeneName],
                                                     "GeneStart":[row.GeneStart],'GeneLength':[row.GeneLength],
                                                     'Index_on_cds':[row['Index_on_cds']],
                                                     'RefCodon':[row.RefCodon],
                                                     'SampleCodon':[codon]})])
    for name, group in haplo.groupby(['GeneName','Index_on_cds','GeneLength','GeneStart','RefCodon']):
        print(group["SampleCodon"].value_counts())
    #haplo = haplo.groupby(['GeneName','Index_on_cds','GeneLength','GeneStart','RefCodon'], as_index=False)['SampleCodon'].value_counts()#.rename(columns={'SampleCodon':'count'}).reset_index()
    #print(haplo)

    sys.exit()
def find_codon_haplotypes(bam_file, poly_df, ref_fasta, mapping_quality=20, codon_read_support=4, verbose=False):
    """bam_file should be the path to the BAM alignment file. poly should be a FinalPolyTable dataframe for
    a particular sample. verbose option prints parameter values. Use mapping quality >= 30 to get only uniquely mapped reads."""
    gene_check = "BT_RS02135"
    console = Console()
    if verbose: console.print(Panel(f"mapping_quality >= {mapping_quality}\ncodon_read_support >= {codon_read_support}", title="[magenta b]find_codon_haplotypes[/magenta b]", expand=False))
    console.print(f'This metagenome (SampleID={poly_df["SampleID"].values[0]}) is one population. ')

    console.print("[b]Parsing BAM file...[/b]")
    if not os.path.isfile(bam_file+'.bai'): os.system(f'samtools index -b {bam_file}')
    bamfile = pysam.AlignmentFile(bam_file, 'rb')
    ref = fasta_to_df(ref_fasta)
    poly_df = poly_df.loc[(poly_df.cds != 'Nan') & (poly_df.Type == "SNP")].sort_values('Pos')
    poly_df[['Pos','start','end']] = poly_df[['Pos','start','end']].astype(float)
    chroms = bamfile.references
    console.print(f'Will explore codons for {poly_df.shape[0]} polys!')
    haplotype_df = pd.DataFrame()
    fail_read_counter = 0
    pass_read_counter = 0
    duplicates = 0
    secondary = 0
    unique = 0
    paired = 0
    for indx, (i, poly) in zip(track(range(poly_df.shape[0]), description='[b]Parsing codon haplotypes...[/b]'),poly_df.iterrows()):
        codon_df = pd.DataFrame()
        pos = int(poly.Pos)
        gene_start = int(poly.start)
        gene_length = len(poly.cds)


        pos_idx = int(poly.Pos) - 1
        if gene_length % 3 == 0: ### work only with coding genes.
            idx_on_gene = int(poly.Pos - float(poly.start)) - 1  # ## ## genbank file uses as locators indexes you need to use to extract this sequence from the reference fasta. (start = location.start-1 end = end)

            if poly.strand == "1.0":
                DNA =  poly.cds
            else:
                DNA = str(Seq(poly["cds"]).reverse_complement())
                #print(f'Ref:{poly.Ref} --> cds_index: {poly.cds[idx_on_gene]} | {ref.loc[ref.Chrom == "NC_004663.1", "seq"].item()[int(poly.Pos)-1]}')
                #print(f'Ref:{poly.Ref} --> cds_index: {poly.cds[idx_on_gene-1]}')
                #print(poly.Ref, poly.cds[idx_on_gene +1])

            if poly.Ref != DNA[idx_on_gene]:
                if idx_on_gene >= 0:
                    logging.warning("Ref does not match cds! Index >= 0")
                    print(f"\tWARNING:: PolyID:{poly.PolyID} | Locus:{poly.GeneName} | strand: {poly.strand}|Start:{poly.start} POS:{poly.Pos} - index on gene {idx_on_gene}, Type: {poly['Type']}, Ancestry:{poly['Ancestry']}")
                    print(f'{ref.loc[ref.Chrom == poly.Chrom+".1","seq"].item()[int(poly.start):pos-1]}{color.RED}{ref.loc[ref.Chrom == poly.Chrom+".1","seq"].item()[pos-1]}{color.END}{ref.loc[ref.Chrom == poly.Chrom+".1","seq"].item()[pos:int(poly.end)]}')
                    print(f'length in fasta: {len(ref.loc[ref.Chrom == poly.Chrom+".1","seq"].item()[int(poly.start):int(poly.end)])}')
                    print(f'\n{DNA[:idx_on_gene]}{color.RED}{DNA[idx_on_gene]}{color.END}{DNA[idx_on_gene+1:]} ')
                    print(f'length of cds: {len(DNA)}')
                    print(f'\t\tRef --> {ref} != {DNA[idx_on_gene]} | GeneStart {ref.loc[ref.Chrom == poly.Chrom+".1","seq"].item()[gene_start]}={DNA[0]}, Alt = {poly.Alt}\n\n')
                else:
                    print(f"\tWARNING:: PolyID:{poly.PolyID} | Locus:{poly['GeneName']} | strand: {poly.strand}|Start:{poly.start} POS:{pos} - index on gene {idx_on_gene}, Type: {poly['Type']}, Ancestry:{poly['Ancestry']}")
                    print(f'{ref.loc[ref.Chrom == poly.Chrom+".1","seq"].item()[int(poly.start):pos-1]}{color.RED}{ref.loc[ref.Chrom == poly.Chrom+".1","seq"].item()[pos-1]}{color.END}{ref.loc[ref.Chrom == poly.Chrom+".1","seq"].item()[pos:int(poly.end)]}')
                    print(f'length in fasta: {len(ref.loc[ref.Chrom == poly.Chrom+".1","seq"].item()[int(poly.start):int(poly.end)])}')
                    print(f'\n{DNA[:idx_on_gene]}{color.RED}{DNA[idx_on_gene]}{color.END}{DNA[idx_on_gene+1:]} ')
                    print(f'length of cds: {len(DNA)}')
                    print(f'\t\tRef --> {ref} != {DNA[idx_on_gene]} | GeneStart {ref.loc[ref.Chrom == poly.Chrom+".1","seq"].item()[gene_start]}={DNA[0]}, Alt = {poly.Alt}\n\n')

                #print(f'{poly.PolyID} | {poly["oldGeneName"]} | start:{poly.start} |{poly.strand}|  {ref.loc[ref.Chrom == "NC_004663.1", "seq"].item()[int(poly.start)]} ?= {poly.cds[0]} ({poly.pseudo})')
            else:
                pass
                #print(f'Ref {poly.Ref} --> Reverse complement {str(Seq(poly.cds).reverse_complement()[idx_on_gene])} | {ref.loc[ref.Chrom == "NC_004663.1", "seq"].item()[int(poly.Pos)-1]}')
                #print(f'Ref {poly.Ref} --> Reverse complement {str(Seq(poly.cds).reverse_complement()[idx_on_gene-1])}')
                #print(f'{poly.Ref}, Reverse complement {str(Seq(poly.cds).reverse_complement()[idx_on_gene+1])})')
                #if poly.Ref != str(Seq(poly.cds).reverse_complement()[idx_on_gene]):
                    #print(f'{poly.PolyID} | {poly["oldGeneName"]} | end:{poly.end} |{poly.strand}|  {ref.loc[ref.Chrom == "NC_004663.1","seq"].item()[int(poly.end)]} ?= {poly.cds[0]}')

            #except:
            #    print(f"WARNING:: REF from VCF does not match genbank annotation for {poly.GeneName}!!\n")
            #    continue ## will skip those positions cause they are unexpected
            idx_on_codon = idx_on_gene % 3  ## this is 0-index --> 0 = 1st , 1 = 2nd, 2 = 3rd
            idx_codon_on_cds = int(np.floor(idx_on_gene / 3)) ## the codon index on the cds
            if verbose:
                print(f'{poly.GeneName}|| Index on Gene: {idx_on_gene}/{gene_length},Codon index {idx_codon_on_cds}/{gene_length /3}, Index on Codon:{idx_on_codon}')

            #print(f'{poly.PolyID} is in codon with index {idx_codon_on_cds} in the {idx_on_codon} position (index)')
            #print(chunks_h(poly.GeneSeq, 3))
            #print(f'Ref codon == {chunks_h(poly.GeneSeq, 3)[idx_codon_on_cds]}, {idx_codon_on_cds} ')
            #print(f'pos_idx={pos_idx}')

            if idx_on_codon == 0:
                codon_start, codon_end = pos_idx, pos_idx + 2
            if idx_on_codon == 1:
                codon_start, codon_end = pos_idx -1, pos_idx + 1
            if idx_on_codon == 2:
                codon_start, codon_end = pos_idx -2, pos_idx

            #print(f'CodonStart={codon_start} and CodonEnd={codon_end}')
            partial_reads = False
            for read in bamfile.fetch(poly.Chrom, codon_start, codon_end):
                if read.mapping_quality >= mapping_quality:# and read.is_secondary is False and read.is_proper_pair is True and read.is_duplicate is False:
                    pass_read_counter += 1
                    if read.mapping_quality == 255: unique+= 1
                    if read.is_secondary: secondary += 1
                    if read.is_proper_pair: paired += 1
                    if read.is_duplicate: duplicates += 1
                    ## get_aligned_pairs returns a tuple with indexes for read, reference + ref sequence (substitutions are lower-case)
                    #print(read.get_aligned_pairs(matches_only=True, with_seq=True))
                    #if len(re.findall('[acgt]', str([x[2] for x in read.get_aligned_pairs(matches_only=True, with_seq=True)]))) >= 2: ## only if there are at least 2 mutations on the read
                    region = [x for x in read.get_aligned_pairs(matches_only=True, with_seq=True) if
                              codon_start <= x[1] <= codon_end]

                    if len(region) == 3: ## makes sure read covers the entire codon
                        codon = ""
                        for position in region:
                            if position[2].islower(): ## mismatch with ref
                                codon += read.query_sequence[position[0]] ## extract alternative base from the read (position[0] is the read index)

                            else:
                                codon += position[2] ## This is the reference base (third element in the tuple for a particular site

                        codon_df = pd.concat([codon_df, pd.DataFrame(
                            {'RefCodon': [chunks_h(DNA, 3)[idx_codon_on_cds]], 'SampleCodon': [codon],
                             'GeneName': [poly.GeneName],
                             'Index_on_cds': [idx_codon_on_cds], "GeneLength": [gene_length],
                             'GeneStart': [gene_start]})], ignore_index=True)
                    else:
                        partial_reads = True

                else:
                    fail_read_counter += 1

        codon_df = codon_df.value_counts().to_frame().reset_index()#.rename(columns={0:'Count'})
        codon_df['PolyID'] = [poly.PolyID for x in codon_df.index]

        codon_df = codon_df.loc[codon_df['count'].ge(codon_read_support)] ## select codons support at least by x reads..
        if codon_df.shape[0] == 1:
            print(f'No polymorphism retrieved for {poly.PolyID}')
            if partial_reads:
                print("\tReads covering partial the codon detected. Maybe mutations are present in those reads.")

        haplotype_df = pd.concat([haplotype_df, codon_df], ignore_index=True)
    console.print(f'Unique={unique}, secondary={secondary}, duplicates={duplicates}, paired={paired}')
    console.print(f"[green b]{pass_read_counter}[/green b] reads made it through filters ([red b]{fail_read_counter}[/red b] failed).")
    if not haplotype_df.empty:
        return haplotype_df.sort_values(["GeneName",'Index_on_cds']).drop_duplicates(['GeneName','SampleCodon','Index_on_cds']) ## if there are two different polys in the same codon they will be represented twice (thus drop duplicates)
    else:
        console.print("[red b]FAILED. All reads have been filtered out!![/red b]")
        sys.exit()

def calculate_π_and_πN_πS(bam_file, poly_df):
    """bam_file should be the path to the BAM alignment file. poly should be a FinalPolyTable dataframe for
        a particular sample.
         For two sequences, the nucleotide diversity is the proportion of sites that contain non-identical nucleotides.
         For a population with two sequences,
         the nucleotide diversity is the nucleotide diversity of the pair times the frequence of each sequence in the population.
           (If there were more than two, the divrsity would be this product summed over all pairwise comparisons in the population)
         """
    console = Console()
    #os.chdir(os.path.dirname(bam_file))
    console.print("[b]Parsing BAM file...[/b]")
    poly_df = poly_df.loc[(poly_df.Type == "SNP") & (poly_df.cds != "Nan")].sort_values("Pos") ### extract polys within genes
    sampleID = poly_df.SampleID.values[0]
    if not os.path.isfile(bam_file + '.bai'): os.system(f'samtools index -b {bam_file}')
    bamfile = pysam.AlignmentFile(bam_file, 'rb')
    pi_site = pd.DataFrame()
    bases = list("ACGT")
    ## estimate π according to Nelson, Chase W.; Hughes, Austin L. (2015) [SNPGenie paper]
    """I have already charachterized SNPs as synonymous or nonsynonymous so π can be seperated into the respective πN and πS 
     where it is calculated for SNPs of either category."""
    for indx, (name, group) in zip(track(range(poly_df.drop_duplicates("Pos").shape[0]), description='[b]Calculating π for sites...[/b]'),
                                   poly_df.groupby("Pos")):
        pos_idx = name - 1
        gene_start = float(group.start.values[0]) ## I have only SNPs in cds in this dataset
        gene_length = len(group.cds.values[0])
        numerator = 0
        numerator_N = 0
        numerator_S = 0
        if gene_length % 3 == 0: ## runs only for coding genes...
            gene = group.cds.values[0]
            idx_on_gene = int(name - gene_start-1)
            idx_on_codon = idx_on_gene % 3  ## this is 0-index --> 0 = 1st , 1 = 2nd, 2 = 3rd
            idx_codon_on_cds = int(np.floor(idx_on_gene / 3)) ## the codon index on the cds

            pile = bamfile.count_coverage(contig=group.Chrom.values[0], start=pos_idx, stop=pos_idx+1, quality_threshold=20) ## count coverage returns an array of size 4 for coverages of A,C,G,T respectively
            denominator = (np.array(pile).sum()**2 - np.array(pile).sum()) / 2


            for (base1,base2), (count1,count2) in zip(itertools.combinations(bases,2), itertools.combinations(np.array(pile).flatten(),2)):
                if count1 > 0 and count2 > 0:
                    codon1 = list(chunks_h(gene,3)[idx_codon_on_cds])
                    codon1[idx_on_codon] = base1
                    codon1 = "".join(codon1)
                    codon2 = list(chunks_h(gene, 3)[idx_codon_on_cds])
                    codon2[idx_on_codon] = base2
                    codon2 = "".join(codon2)
                    Ndiff, Sdiff = Codon(codon1).Nd_Ns_for_codon_pair(Codon(codon2)) ## because they differ by 1 by default - either 1,0 or 0,1 is expected
                    numerator_N += Ndiff * count1 * count2
                    numerator_S += Sdiff * count1 * count2
                    numerator += numerator_N + numerator_S
            pi  = numerator / denominator
            piN = numerator_N / denominator
            piS = numerator_S / denominator
            pi_site = pd.concat([pi_site,pd.DataFrame({'SampleID':[sampleID], 'Chrom':[group.Chrom.values[0]],'Pos':[name],"GeneStart":[gene_start],"GeneLength":[gene_length],
                                               "Codon_index":[idx_codon_on_cds],'GeneName':[group.GeneName.values[0]],
                                               "pi":[pi], "piN":[piN], "piS":[piS]})], ignore_index=True)
    pi_codon = pd.DataFrame()
    for name, group in pi_site.groupby(['GeneName','Codon_index']): ### this needs to be fixed in the future
        region_df = (group[['pi','piN','piS']].agg('sum')/3)## 3 is the length of the codon
        pi_codon = pd.concat([pi_codon,pd.DataFrame({'SampleID':[sampleID],"GeneStart":[group.GeneStart.values[0]],"GeneLength":[group.GeneLength.values[0]],
                                                 "GeneName":[name[0]], "Codon_index":[name[1]],
                                 "codon_pi":[region_df.loc['pi']], 'codon_piN':[region_df.loc['piN']], 'codon_piS':[region_df.loc['piS']]})], ignore_index=True)
    pi_gene = pd.DataFrame()
    for name, group in pi_site.groupby(['GeneName']):
        region_df = (group[['pi', 'piN', 'piS']].agg('sum') / group.GeneLength.values[0])  ## divide by the length of the gene
        pi_gene = pd.concat([pi_gene,pd.DataFrame({'SampleID':[sampleID], "GeneStart":[group.GeneStart.values[0]],"GeneLength":[group.GeneLength.values[0]],"GeneName": [name],
                                                 "gene_pi": [region_df.loc['pi']], 'gene_piN': [region_df.loc['piN']],
                                                 'gene_piS': [region_df.loc['piS']]})], ignore_index=True)


    return pi_site, pi_codon, pi_gene

def calculate_π_and_πΝ_πS_isolates(poly_df):
    console = Console()
    base_count = pd.DataFrame()
    #DNA_bases = {}
    poly_df = poly_df.loc[(poly_df.Type == "SNP") & (poly_df.cds != "Nan")].sort_values("Pos")
    N_isolates = poly_df.drop_duplicates("SampleID").shape[0]
    #print(poly_df.groupby(['Pos'], as_index=True).agg({'Alt':'value_counts'}).rename(columns={"Alt":"Counts"}).reset_index() )
    pi_site = pd.DataFrame()
    for name, group in poly_df.groupby("Pos"):
        count_df = group["Alt"].value_counts().to_frame().reset_index().rename(columns={'Alt':'base'})
        print(group.GeneName)
        alt_count = count_df['count'].sum()
        gene_start = float(group.start.values[0]) ## I have only SNPs in cds in the dataset
        gene = group.cds.values[0]
        gene_length = len(gene)
        if gene_length % 3 == 0:
            ## add count for the reference base N_isolates - Count for all Alt
            count_df = pd.concat([count_df,pd.DataFrame({"base":[group['Ref'].values[0]], "count":[N_isolates - alt_count] })],ignore_index=True)
            print(count_df)
            ### find synonymous pairs and non synonymous..
            idx_on_gene = int(name - float(gene_start-1))
            idx_on_codon = idx_on_gene % 3  ## this is 0-index --> 0 = 1st , 1 = 2nd, 2 = 3rd
            idx_codon_on_cds = int(np.floor(idx_on_gene / 3)) ## the codon index on the cds
            print(idx_codon_on_cds, idx_on_gene, gene_length)
            count_df.set_index('base', inplace=True)
            numerator = 0
            numerator_N = 0
            numerator_S = 0
            denominator = (N_isolates**2 - N_isolates) / 2
            for base1,base2 in itertools.combinations(count_df.index,2):
                codon1 = list(chunks_h(gene,3)[idx_codon_on_cds])
                print(codon1)
                codon1[idx_on_codon] = base1
                print(codon1)
                codon1 = "".join(codon1)
                codon2 = list(chunks_h(gene, 3)[idx_codon_on_cds])
                codon2[idx_on_codon] = base2
                codon2 = "".join(codon2)
                Ndiff, Sdiff = Codon(codon1).Nd_Ns_for_codon_pair(Codon(codon2)) ## because they differ by 1 by default - either 1,0 or 0,1 is expected
                numerator_N += Ndiff * count_df.loc[base1, 'count'] * count_df.loc[base2, "count"]
                numerator_S += Sdiff * count_df.loc[base1, 'count'] * count_df.loc[base2, "count"]
                numerator += numerator_N + numerator_S
            pi  = numerator / denominator
            piN = numerator_N / denominator
            piS = numerator_S / denominator
            pi_site = pd.concat([pi_site,pd.DataFrame({"GeneStart":[gene_start], "GeneLength":[gene_length],'Pos':[name],"Codon_index":[idx_codon_on_cds],'GeneName':[group.GeneName.values[0]],
                                               "pi":[pi], "piN":[piN], "piS":[piS]})], ignore_index=True)
        ### to estimate over a region I need to average over the sites of the regions sum1-L(pi)/L , where L length of region.
    pi_codon = pd.DataFrame()
    for name, group in pi_site.groupby(['GeneName', 'Codon_index']):
        region_df = (group[['pi', 'piN', 'piS']].agg('sum') / 3)  ## 3 is the length of the codon
        pi_codon = pd.concat([pi_codon,
            pd.DataFrame({"GeneStart": [group.GeneStart.values[0]], "GeneLength": [group.GeneLength.values[0]],
                          "GeneName": [name[0]], "Codon_index": [name[1]],
                          "pi": [region_df.loc['pi']], 'piN': [region_df.loc['piN']], 'piS': [region_df.loc['piS']]})], ignore_index=True)
    pi_gene = pd.DataFrame()
    for name, group in pi_site.groupby(['GeneName']):
        region_df = (group[['pi', 'piN', 'piS']].agg('sum') / group.GeneLength.values[
            0])  ## divide by the length of the gene
        pi_gene = pd.concat([pi_gene, pd.DataFrame(
            {"GeneStart": [group.GeneStart.values[0]], "GeneLength": [group.GeneLength.values[0]],
             "GeneName": [name],
             "pi": [region_df.loc['pi']], 'piN': [region_df.loc['piN']],
             'piS': [region_df.loc['piS']]})], ignore_index=True)

    return pi_site, pi_codon, pi_gene


def pNps(Nd, N, Sd, S, printout=False):
    if printout: print(f'Nd={Nd}, N={N}, Sd={Sd}, S={S}')
    if N > 0:
        pN = Nd / N
    else:
        pN = np.nan ## cannot be defined
    if S > 0:
        pS = Sd / S
    else:
        #print('Expected Synonymous subsitutions is zero!!!')
        pS = np.nan ## cannot be defined
    if printout: print(f'pN={pN}, pS={pS}')
    return pN, pS

def dN_dS_from_pN_pS(pN, pS, printout=False):
    """This applies the correction from Jukes and Cantor (1969)"""
    if printout: print(f'pN={pN}, pS={pS}')
    dN, dS = np.nan, np.nan
    if isinstance(pN, float) and isinstance(pS,float):
        if pN < 0.75 and pS < 0.75:
            dN = -(3/4)*np.log(1 - (4*pN/3))
            dS = -(3/4)*np.log(1 - (4*pS/3))
        #elif pN < 0.75:
        #    logging.warning(f'pN >= 0.75, dN cannot be defined!')
        #elif pS < 0.75:
        #    logging.warning(f'pS >= 0.75, dS cannot be defined!')
    return dN,dS

class Nucleotide:
    "DNA nucleotide"
    def __init__(self, seq):
        self.seq = seq
        self.DNA_bases = {'A':'T','T':'A','C':'G','G':'C'}
        self.chemistry_group = {'A':'purine','G':'purine','C':'pyrimidine','T':'pyrimidine'}
    def __repr__(self):
        return f"Codon('{self.seq}')"
    def __str__(self):
        return self.seq

    @property
    def complementary(self):
        return self.DNA_bases[self.seq]
    @property
    def chemistry(self):
        return self.chemistry_group[self.seq]
class Codon:
    """Seq should be a codon sequence """

    def __init__(self, seq):
        self.seq = seq
        self.DNA_bases = list("ACGT")
        assert len(self.seq) == 3
    def __repr__(self):
        return f"Codon('{self.seq}')"
    def __str__(self):
        return self.seq
    @property
    def all_mutational_events(self):
        """Mutates a codon with all possible mutations (3*3=9) and returns a list of the mutated codons"""
        codon = list(self.seq)
        mutational_events = []
        for i in range(len(codon)):
            for base in self.DNA_bases:
                if base != codon[i]: ## mutate for the rest 3 DNA bases (Daa!)
                    mutant = codon[:]
                    mutant[i] = base
                    mutational_events.append("".join(mutant))
        return mutational_events

    @property
    def N_sites(self):
        counter = 0
        for mutant in self.all_mutational_events:
            if Seq(self.seq).translate(table='11') != Seq(mutant).translate(table='11'):
                counter += 1/3
        return counter
    @property
    def S_sites(self):
        return 3 - self.N_sites
    @property
    def aminoacid(self):
        return Seq(self.seq).translate(table='11')

    #@staticmethod
    def Nd_Ns_for_codon_pair(self, other):
        """Find evolutionary paths and calculate synonymous and non synonynmous sites between two codons.
        Returns total number of averaged over all paths observed nonsynonymous and synonymous mutations between two codons (sample - reference)"""
        N_paths = 0
        Nd = 0
        Sd = 0
        distance = np.sum([1 for x,y in zip(self.seq, other.seq) if x != y])
        if distance == 0:
            N_paths = 1 ## cannot divide with zero

        elif distance == 1:
            N_paths = 1
            if Seq(self.seq).translate(table='11') != Seq(other.seq).translate(table='11'):
                Nd += 1
            else:
                Sd += 1

        elif distance > 1:
            N_paths = 0
            Nd_values = []
            Sd_values = []
            for path in itertools.permutations([(i,y) for i,x,y in zip(range(3),self.seq, other.seq) if x != y]): ## loop over paths
                Nd_path = 0
                Sd_path = 0
                N_paths += 1
                mutant = list(self.seq) ## start of evolutionary path
                for mutation in path: ## loop over mutational events in the path
                    mutant_next = mutant[:]
                    mutant_next[mutation[0]] = mutation[1]

                    #print(f'{"".join(mutant)} --> {"".join(mutant_next)} || {Seq("".join(mutant)).translate(table="11")} --> {Seq("".join(mutant_next)).translate(table="11")}')

                    if Seq("".join(mutant)).translate(table='11') != Seq("".join(mutant_next)).translate(table='11'):
                        Nd_path += 1
                    else:
                        Sd_path += 1
                    mutant = mutant_next
                Nd_values.append(Nd_path)
                Sd_values.append(Sd_path)
                #print(f'Nd = {Nd_path} , Sd = {Sd_path}')
            Nd, Sd = sum(Nd_values), sum(Sd_values)
        Nd, Sd = Nd / N_paths, Sd / N_paths  ### to average across evolutionary paths.
        return Nd, Sd

class CDS:

    def __init__(self, seq):
        self.seq = seq
        self.console = Console()
    def __repr__(self):
        return self.codons
    def __str__(self):
        return self.seq

    @property
    def codons(self):
        try:
            len(self.seq) % 3 == 0
        except:
            self.console.log("CDS sequence should be a multiple of 3!")
        return [Codon(x) for x in chunks_h(self.seq, 3)] # return a list of Codon classes from a sequence

    @property
    def N_non_syn_sites(self):
        counter = 0
        for codon in self.codons:
            counter += codon.N_sites
        return counter

    @property
    def N_syn_sites(self):
        counter = 0
        for codon in self.codons:
            counter += codon.S_sites
        return counter

    def calculate_Nd_Sd_for_CDS(self, other):
        Nd = 0
        Sd = 0
        for codon1, codon2 in zip(self.codons, other.codons):
            codon_counter = codon1.Nd_Ns_for_codon_pair(codon2)
            Nd += codon_counter[0] ## returns Nd, Sd
            Sd += codon_counter[1]
        return Nd, Sd


def calculate_dNdS_from_HaploTable_test(haplo, codon_read_support=4):
    """HaploTable is the output of the find_codon_haplotypes function. Returns a table with pNpS and dNdS values for codon and gene level. """

    codon_df = pd.DataFrame()
    for i, (name, group) in zip(track(range(haplo.drop_duplicates(['GeneName', 'RefCodon', 'Index_on_cds']).shape[0]),
                                      description='[b]Calculating codon dN/dS[/b]'),
                                haplo.groupby(['GeneName', 'RefCodon', 'Index_on_cds'])):

        ### keep only the most abundant polymorphic state (similar to breseq behaviour), (get the first two including REF)
        keep_one_alt = False
        if keep_one_alt:
            group = group.sort_values('count', ascending=False).head(2)

        ref_codon = name[1]
        gene_length = group.GeneLength.values[0]
        gene_start = group.GeneStart.values[0]

        ### exclude reference states from the cal
        # total_coverage = group.loc[(group.SampleCodon != group.RefCodon) & (group['count'].astype(int) >= codon_read_support)][
        #    'count'].astype(int).sum()
        total_coverage = group.loc[group['count'].astype(int) >= codon_read_support]['count'].astype(int).sum()

        ### keep only the most abundant polymorphic state (similar to breseq behaviour), (get the first two including REF)
        keep_one_alt = False
        if keep_one_alt:
            group = group.sort_values('count', ascending=False).head(2)

        eN, eS = Codon(ref_codon).N_sites, Codon(
            ref_codon).S_sites  ## this is the calculation for the expected non synonymous and synonymous substitutions
        # if N > 0 and S > 0: ### those sites cannot be included

        # if S > 0: ## otherwise pS and dS cannot be defined

        eN = eN * total_coverage
        eS = eS * total_coverage

        codon_pN = 0
        codon_pS = 0
        codon_oN = 0  ### this is corrected for all alternative paths (can be float)
        codon_oS = 0  ### this is corrected for all alternative paths.(can be float)
        N_count = 0
        S_count = 0
        # if group.shape[0] > 1 and group.RefCodon.values[0]: ## if there is one entry it should be a fixed mutation
        for i, row in group.iterrows():  ## this iterates over all observed states of the codon in the data
            if row.SampleCodon != ref_codon and int(
                    row['count']) >= codon_read_support:  ### filter out states supported only by this number of reads

                alt_coverage = int(row['count'])  ## this is the coverage of each state
                ## this is the observed substitutions (for one state comparison)
                Nd, Sd = Codon(row.SampleCodon).Nd_Ns_for_codon_pair(Codon(ref_codon))
                codon_oS += Sd * alt_coverage#(alt_coverage / total_coverage)
                codon_oN += Nd * alt_coverage#(alt_coverage / total_coverage)


                ## this is for raw observations
                for i, x, y in zip(range(3), row.SampleCodon, ref_codon):
                    # print(i)
                    if x != y:
                        mutant = ref_codon
                        mutant = list(mutant)
                        mutant[i] = x
                        mutant = "".join(mutant)
                        if Codon(mutant).aminoacid != Codon(ref_codon).aminoacid:
                            N_count += 1
                        else:
                            S_count += 1

        if eN > 0:
            codon_pN = codon_oN / eN
        else:
            codon_pN = np.nan
        if eS > 0:
            codon_pS = codon_oS / eS
        else:
            codon_pS = np.nan

        ## apply the correction to convert to dNdS
        codon_dN, codon_dS = dN_dS_from_pN_pS(codon_pN, codon_pS, printout=False)

        codon_df = pd.concat([codon_df, pd.DataFrame(
            {"SampleID": [group.SampleID.values[0]], 'Codon_N_count': [N_count], 'Codon_S_count': [S_count],
             'Codon_oN': [round(codon_oN,5)], 'Codon_oS': [round(codon_oS,5)], 'Codon_eN': [round(eN,5)], 'Codon_eS': [round(eS,5)],
             'Codon_pN': [round(codon_pN,5)], "Codon_pS": [round(codon_pS,5)],
             "Codon_dN": [round(codon_dN,5)], 'Codon_dS': [round(codon_dS,5)],
             "Index_on_cds": [name[2]], "RefCodon": [ref_codon], "GeneName": [name[0]],
             "GeneStart": [gene_start], "GeneLength": [gene_length],
             })], ignore_index=True)

    """ Fractions can be defined only if pS > 0 and dS > 0 for pN/pS and dNdS respectively"""
    codon_df['Codon_pNpS'] = [np.nan for x in codon_df.index]
    for i, row in codon_df.iterrows():
        if row['Codon_pS'] > 0 and row['Codon_pS'] != np.nan:
            codon_df.loc[i, 'Codon_pNpS'] = row['Codon_pN'] / row['Codon_pS']
        else:
            codon_df.loc[i, 'Codon_pNpS'] = np.nan

    codon_df['Codon_dNdS'] = [np.nan for x in codon_df.index]
    for i, row in codon_df.iterrows():
        if row['Codon_dS'] > 0 and row['Codon_dS'] != np.nan:
            codon_df.loc[i, 'Codon_dNdS'] = row['Codon_dN'] / row['Codon_dS']
        else:
            codon_df.loc[i, 'Codon_dNdS'] = np.nan

    """ Statistics for the gene """
    gene_df = pd.DataFrame()
    # codN = 0
    if not codon_df.empty:
        for name, group in codon_df.groupby('GeneName'):
            # codN += 1

            gene_pN = np.nan
            gene_pS = np.nan
            gene_pNpS = np.nan
            gene_dNdS = np.nan
            gene_dN = np.nan
            gene_dS = np.nan
            gene_oN = group.Codon_oN.sum()
            gene_oS = group.Codon_oS.sum()
            gene_eN = group.Codon_eN.sum()
            gene_eS = group.Codon_eS.sum()
            gene_N_count = group.Codon_N_count.sum()
            gene_S_count = group.Codon_S_count.sum()

            # sum pN and pS over codons with pS > 0 and divide by the the number of valid codons (pS > 0)
            region_df = group[['Codon_pN', 'Codon_pS']].agg('sum') / group.shape[0]  ## each row is a codon

            ### try summing observed and expected values instead
            #region_df = group[['Codon_oN', 'Codon_oS','Codon_eN','Codon_eS']].agg('sum').to_frame().T

            #region_df.rename(columns ={x:y for x,y in zip(region_df.columns.values,[x.replace('Codon', "Gene") for x in region_df.columns])}, inplace=True)
            #region_df.reset_index(inplace=True)
            #region_df['Gene_pN'] = [np.nan for x in region_df.index]
            #region_df['Gene_pS'] = [np.nan for x in region_df.index]
            #region_df['Gene_pNpS'] = [np.nan for x in region_df.index]


            if gene_eN > 0:
                gene_pN = gene_oN / gene_eN
            if gene_eS > 0:
                gene_pS = gene_oS / gene_eS



            if gene_pS > 0:
                gene_pNpS = gene_pN / gene_pS
                gene_dN, gene_dS = dN_dS_from_pN_pS(gene_pN, gene_pS)
                gene_dNdS = gene_dN / gene_dS
            gene_df = pd.concat([gene_df, pd.DataFrame(
                {"SampleID": [group.SampleID.values[0]], "GeneStart": [group.GeneStart.values[0]],
                 "GeneLength": [group.GeneLength.values[0]],
                 "GeneName": [name],
                 "Gene_eN":[gene_eN], 'Gene_eS':[gene_eS],
                 "Gene_oN": [gene_oN], "Gene_oS": [gene_oS], 'Gene_N_count': [gene_N_count],
                 'Gene_S_count': [gene_S_count],
                 "Gene_pN": [gene_pN], 'Gene_pS': [gene_pS], 'Gene_pNpS': [gene_pNpS],
                 "Gene_dN": [gene_dN], 'Gene_dS': [gene_dS], "Gene_dNdS": [gene_dNdS]})], ignore_index=True)

            # if gene_eN > 0:
            #     gene_pN = gene_oN / gene_eN
            # if gene_eS > 0:
            #     gene_pS = gene_oS / gene_eS
            #
            # if gene_pS > 0:
            #     gene_pNpS = gene_pN / gene_pS
            #
            # if gene_pS > 0:
            #     gene_dN, gene_dS = dN_dS_from_pN_pS(gene_pN, gene_pS)
            #     gene_dNdS = gene_dN / gene_dS
            #
            #
            # gene_df = pd.concat([gene_df, pd.DataFrame(
            #     {"SampleID":[group.SampleID.values[0]],"GeneStart": [group.GeneStart.values[0]], "GeneLength": [group.GeneLength.values[0]],
            #      "GeneName": [name], "Gene_eN":[gene_eN], 'Gene_eS':[gene_eS],
            #      "Gene_oN":[gene_oN],"Gene_oS":[gene_oS], 'Gene_N_count':[gene_N_count],'Gene_S_count':[gene_S_count],
            #      "Gene_pN": [gene_pN], 'Gene_pS': [gene_pS],'Gene_pNpS':[gene_pNpS],
            #      "Gene_dN":[gene_dN], 'Gene_dS':[gene_dS], "Gene_dNdS":[gene_dNdS]})], ignore_index=True)
    else:
        logging.warning(f'No detected polymorphisms')
    return codon_df, gene_df


def calculate_dNdS_from_HaploTable(haplo, codon_read_support = 4):
    """HaploTable is the output of the find_codon_haplotypes function. Returns a table with pNpS and dNdS values for codon and gene level. """
    codon_df = pd.DataFrame()
    for i,(name, group) in zip(track(range(haplo.drop_duplicates(['GeneName','RefCodon','Index_on_cds']).shape[0]), description='[b]Calculating codon dN/dS[/b]'),
                               haplo.groupby(['GeneName','RefCodon','Index_on_cds'])):

        ### keep only the most abundant polymorphic state (similar to breseq behaviour), (get the first two including REF)
        keep_one_alt = True
        if keep_one_alt:
            group = group.sort_values('count', ascending=False).head(2)



        ref_codon = name[1]
        gene_length = group.GeneLength.values[0]
        gene_start = group.GeneStart.values[0]
        eN, eS = Codon(ref_codon).N_sites, Codon(ref_codon).S_sites ## this is the calculation for the expected non synonymous and synonymous substitutions
        #if N > 0 and S > 0: ### those sites cannot be included

        #if S > 0: ## otherwise pS and dS cannot be defined

        # exclude reference states from the weighted average
        total_coverage = group.loc[(group.SampleCodon != group.RefCodon) & (group['count'].astype(int) >= codon_read_support)]['count'].astype(int).sum()
        codon_pN = 0
        codon_pS = 0
        codon_oN = 0 ### this is corrected for all alternative paths (can be float)
        codon_oS = 0 ### this is corrected for all alternative paths.(can be float)
        N_count = 0
        S_count = 0
        #if group.shape[0] > 1 and group.RefCodon.values[0]: ## if there is one entry it should be a fixed mutation
        for i, row in group.iterrows(): ## this iterates over all observed states of the codon in the data
            if row.SampleCodon != ref_codon and int(row['count']) >= codon_read_support: ### filter out states supported only by this number of reads

                alt_coverage = int(row['count']) ## this is the coverage of each state
                ## this is the observed substitutions (for one state comparison)
                Nd, Sd = Codon(row.SampleCodon).Nd_Ns_for_codon_pair(Codon(ref_codon))
                codon_oS += Sd * (alt_coverage/total_coverage)
                codon_oN += Nd * (alt_coverage/total_coverage)

                # this is the raw observations
                for i,x,y in zip(range(3),row.SampleCodon, ref_codon):
                    #print(i)
                    if x != y:
                        mutant = ref_codon
                        mutant = list(mutant)
                        mutant[i] = x
                        mutant = "".join(mutant)
                        if Codon(mutant).aminoacid != Codon(ref_codon).aminoacid:
                            N_count += 1
                        else:
                            S_count += 1

        if eN > 0:
            codon_pN = codon_oN / eN
        else:
            codon_pN = np.nan
        if eS > 0:
            codon_pS = codon_oS / eS
        else:
            codon_pS = np.nan


        ## apply the correction to convert to dNdS
        codon_dN, codon_dS = dN_dS_from_pN_pS(codon_pN, codon_pS, printout=False)

        codon_df = pd.concat([codon_df,pd.DataFrame({"SampleID":[group.SampleID.values[0]],'Codon_N_count':[N_count],'Codon_S_count':[S_count],
                                                     'Codon_oN':[codon_oN],'Codon_oS':[codon_oS],'Codon_eN':[eN], 'Codon_eS':[eS],
                                                     'Codon_pN':[codon_pN], "Codon_pS":[codon_pS],
                                                 "Codon_dN":[codon_dN],'Codon_dS':[codon_dS],
                                                 "Index_on_cds":[name[2]],"RefCodon":[ref_codon],"GeneName":[name[0]],
                                                     "GeneStart":[gene_start], "GeneLength":[gene_length],
                                                     })], ignore_index=True)

    """ Fractions can be defined only if pS > 0 and dS > 0 for pN/pS and dNdS respectively"""
    codon_df['Codon_pNpS'] = [np.nan for x in codon_df.index]
    for i,row in codon_df.iterrows():
        if row['Codon_pS'] > 0 and row['Codon_pS'] != np.nan:
            codon_df.loc[i, 'Codon_pNpS'] = row['Codon_pN'] / row['Codon_pS']
        else:
            codon_df.loc[i, 'Codon_pNpS'] = np.nan

    codon_df['Codon_dNdS'] = [np.nan for x in codon_df.index]
    for i,row in codon_df.iterrows():
        if row['Codon_dS'] > 0 and row['Codon_dS'] != np.nan:
            codon_df.loc[i, 'Codon_dNdS'] = row['Codon_dN'] / row['Codon_dS']
        else:
            codon_df.loc[i, 'Codon_dNdS'] = np.nan



    """ Statistics for the gene """
    gene_df = pd.DataFrame()
    #codN = 0
    if not codon_df.empty:
        for name, group in codon_df.groupby('GeneName'):
            #codN += 1


            gene_pN = np.nan
            gene_pS = np.nan
            gene_pNpS = np.nan
            gene_dNdS = np.nan
            gene_dN = np.nan
            gene_dS = np.nan
            gene_oN = group.Codon_oN.sum()
            gene_oS = group.Codon_oS.sum()
            gene_eN = group.Codon_eN.sum()
            gene_eS = group.Codon_eS.sum()
            gene_N_count = group.Codon_N_count.sum()
            gene_S_count = group.Codon_S_count.sum()

            # sum pN and pS over codons with pS > 0 and divide by the the number of valid codons (pS > 0)
            region_df = group[['Codon_pN', 'Codon_pS']].agg('sum') / group.shape[0]  ## each row is a codon
            region_df.index = [x.replace('Codon', "Gene") for x in region_df.index]
            if region_df.loc['Gene_pS'] > 0:
                gene_pNpS = region_df.loc['Gene_pN'] / region_df.loc['Gene_pS']
                gene_dN, gene_dS = dN_dS_from_pN_pS(region_df.loc['Gene_pN'],region_df.loc['Gene_pS'] )
                gene_dNdS = gene_dN / gene_dS
            gene_df = pd.concat([gene_df, pd.DataFrame(
                {"SampleID": [group.SampleID.values[0]], "GeneStart": [group.GeneStart.values[0]],
                 "GeneLength": [group.GeneLength.values[0]],
                 "GeneName": [name],
                 "Gene_oN": [gene_oN], "Gene_oS": [gene_oS], 'Gene_N_count': [gene_N_count],
                 'Gene_S_count': [gene_S_count],
                 "Gene_pN": [gene_pN], 'Gene_pS': [gene_pS], 'Gene_pNpS': [gene_pNpS],
                 "Gene_dN": [gene_dN], 'Gene_dS': [gene_dS], "Gene_dNdS": [gene_dNdS]})], ignore_index=True)




            # if gene_eN > 0:
            #     gene_pN = gene_oN / gene_eN
            # if gene_eS > 0:
            #     gene_pS = gene_oS / gene_eS
            #
            # if gene_pS > 0:
            #     gene_pNpS = gene_pN / gene_pS
            #
            # if gene_pS > 0:
            #     gene_dN, gene_dS = dN_dS_from_pN_pS(gene_pN, gene_pS)
            #     gene_dNdS = gene_dN / gene_dS
            #
            #
            # gene_df = pd.concat([gene_df, pd.DataFrame(
            #     {"SampleID":[group.SampleID.values[0]],"GeneStart": [group.GeneStart.values[0]], "GeneLength": [group.GeneLength.values[0]],
            #      "GeneName": [name], "Gene_eN":[gene_eN], 'Gene_eS':[gene_eS],
            #      "Gene_oN":[gene_oN],"Gene_oS":[gene_oS], 'Gene_N_count':[gene_N_count],'Gene_S_count':[gene_S_count],
            #      "Gene_pN": [gene_pN], 'Gene_pS': [gene_pS],'Gene_pNpS':[gene_pNpS],
            #      "Gene_dN":[gene_dN], 'Gene_dS':[gene_dS], "Gene_dNdS":[gene_dNdS]})], ignore_index=True)
    else:
        logging.warning(f'No detected polymorphisms')
    return codon_df, gene_df

def codon_haplotypes_from_bam(bam, sampleid, poly_df, mapping_quality=20, qc=20, imputation = False):
    if not os.path.isfile(bam + '.bai'): os.system(f'samtools index -b {bam}')
    bamfile = pysam.AlignmentFile(bam, 'rb')
    haplotype_df = pd.DataFrame()
    breseq_only = []
    bam_only = []
    no_polymorphic_codon = 0
    missing_polys = []
    if not poly_df.empty:
        for indx, (name, group) in zip(
                track(range(len(poly_df.groupby("GeneName").groups)),
                      description=f'[b]Parsing codon haplotypes ({sampleid})...[/b]'),
                poly_df.groupby('GeneName')):
            cds = group.cds.values[0]
            gene_length = len(cds)

            if gene_length % 3 == 0:  ### work only with coding genes.
                codon_df = pd.DataFrame()
                gene_start = int(group.start.values[0])
                gene_end = int(group.end.values[0])
                strand = group.strand.values[0]
                if strand == "1.0":
                    strand = '+'
                    #DNA = cds
                else:
                    #DNA = str(Seq(cds).reverse_complement())
                    strand = "-"
                #detected_polys = []

                ### split CDS to codons and their coordinates
                codon_coords = find_codon_coords(cds, gene_start, return_all=True)
                for i,row in group.iterrows():
                    #print(f"{name} | {row.oldGeneName} ({strand}) [{gene_start} - {gene_end}] | {row.PolyID}\n{row.translation}")
                    codon_start, codon_end = codon_coords[int(float(row.CodonIndex))]

                    for read in bamfile.fetch(group.Chrom.values[0], codon_start, codon_end):
                        if read.mapping_quality >= mapping_quality and read.is_secondary is False and read.is_duplicate is False and read.is_qcfail is False:  # and read.is_proper_pair is True:
                            region = [x for x in read.get_aligned_pairs(matches_only=True, with_seq=True) if
                                      codon_start <= x[1] <= codon_end]
                            codon = {codon_start: "-", codon_start + 1: "-", codon_start + 2: "-"}

                            for position in region:
                                if read.query_qualities[position[0]] >= qc:
                                    if position[2].islower():  ## mismatch with ref
                                        codon[position[1]] = read.query_sequence[position[0]]
                                    else:
                                        codon[position[1]] = position[2]


                            codon_df = pd.concat([codon_df, pd.DataFrame(
                                {'SampleID': [sampleid], 'RefCodon': [str(CDS(cds).codons[int(float(row.CodonIndex))])],
                                 'SampleCodon': [codon],
                                 'GeneName': [name],
                                 'Index_on_cds': [int(float(row.CodonIndex))], "GeneLength": [gene_length],
                                 'GeneStart': [gene_start], 'CodonStart': [codon_start], 'Associated_poly':[row.PolyID]})], ignore_index=True)

                ## impute missing bases (not covered by read) with reference DNA base
                def impute_missing_with_ref(row):
                    samplecodon = list(row.SampleCodon.values())
                    for i, (x, y) in enumerate(zip(samplecodon, row.RefCodon)):
                        if x == '-': samplecodon[i] = y
                    return "".join(samplecodon)

                if not codon_df.empty:
                    impute = imputation
                    if not impute:
                        remove_incomplete = True
                    if impute:
                        codon_df["SampleCodon"] = codon_df.apply(lambda x: impute_missing_with_ref(x), axis=1)
                    else:
                        codon_df["SampleCodon"] = codon_df["SampleCodon"].apply(lambda x: ''.join(x.values()))
                        if remove_incomplete:
                            codon_df = codon_df.loc[~codon_df.SampleCodon.str.contains("-")]


                    ## count different codons
                    codon_df = codon_df.value_counts().to_frame().reset_index()  # .rename(columns={0:'Count'})
                    #print(codon_df);sys.exit()
                    ## select codons support at least by x reads..
                    #codon_df = codon_df.loc[codon_df['count'].ge(codon_read_support)]

                    # codon_df['PolyID'] = [poly.PolyID for x in codon_df.index]
                    #detected_polys = pd.Series(detected_polys).value_counts().reset_index()
                    #detected_polys = detected_polys.loc[detected_polys['count'] >= codon_read_support]['index'].tolist()
                    ### maybe filter for Freq 0.01 as well?
                    #breseq_only += [x for x in group.PolyID.values if x not in detected_polys]
                    #bam_only += [x for x in detected_polys if x not in group.PolyID.values]

                    for codoname,codon_group in codon_df.groupby('Index_on_cds'):
                        #print(codon_group.loc[codon_group.SampleCodon != codon_group.RefCodon])
                        if codon_group.loc[codon_group.SampleCodon != codon_group.RefCodon].empty:
                            #print(f'{name} Could not retrieve mutation for this codon! (detected previously with breseq)')
                            no_polymorphic_codon += 1
                            missing_polys.append(codon_group['Associated_poly'].values[0])

                    #codon_df['Freq'] = codon_df['count'].apply(lambda x: x / codon_df['count'].sum())

                    # print(f'Gene: {name}| {group.oldGeneName.values[0]} ({strand})')
                    # if len(detected_polys) != len(group.PolyID.values):
                    #     print(sorted(detected_polys))
                    #     print(group.PolyID.sort_values().values)
                    #     print(codon_df[['RefCodon','SampleCodon','Index_on_cds','CodonStart','GeneStart','count']].sort_values('CodonStart'))
                    #     print('\n')

                    haplotype_df = pd.concat([haplotype_df, codon_df], ignore_index=True)
    print(f'{no_polymorphic_codon} codons are not polymorphic after filtering reads in the alignment! (despite called with breseq)')


    return haplotype_df, missing_polys

def submit_dNdS_from_metagenome_job(sampleID, bamfile, poly_df, ref, outdir, mapping_quality=20):
    ram =8000
    cpus = 1
    with open(f'{sampleID}_dNdS_slurm_job.sh','w') as sl:
        sl.write("#!/bin/bash\n")
        sl.write("#\n")
        sl.write("#SBATCH --job-name=dNdS\n")
        sl.write(f"#SBATCH --mem={ram}\n")
        sl.write(f"#SBATCH --cpus-per-task={cpus}\n")
        #sl.write("#SBATCH --array=1-"+str(array_length)+'\n')
        sl.write('module load samtools\n')
        #sl.write(f'rsync -a {bamfile} {poly_df} {ref} $TMPDIR\n')
        sl.write(f'dNdS.py -i {bamfile}  -t $TMPDIR -r {ref} -p {poly_df} -o {outdir} --mq {mapping_quality}\n')
        sl.close()
    os.system(f'sbatch {sampleID}_dNdS_slurm_job.sh')

def submit_fst_from_metagenomes_job(bamDirs, poly_df, outdir, ancestralMode=False):
    ram =16000
    cpus = 1
    with open(f'fst_slurm_job.sh','w') as sl:
        sl.write("#!/bin/bash\n")
        sl.write("#\n")
        sl.write("#SBATCH --job-name=fst_metagenomes\n")
        sl.write(f"#SBATCH --mem={ram}\n")
        sl.write(f"#SBATCH --cpus-per-task={cpus}\n")
        #sl.write("#SBATCH --array=1-"+str(array_length)+'\n')
        sl.write('module load samtools\n')
        #sl.write(f'rsync -a {bamfile} {poly_df} {ref} $TMPDIR\n')
        if outdir:
            if not ancestralMode:
                sl.write(f'fst.py -d {" ".join(bamDirs)}  -t $TMPDIR -p {poly_df} -o {outdir} \n')
            else:
                sl.write(f'fst.py -d {" ".join(bamDirs)}  -t $TMPDIR -p {poly_df} -o {outdir} -a \n')
        else:
            if not ancestralMode:
                sl.write(f'fst.py -d {" ".join(bamDirs)}  -t $TMPDIR -p {poly_df}\n')
            else:
                sl.write(f'fst.py -d {" ".join(bamDirs)}  -t $TMPDIR -p {poly_df} -a \n')
        sl.close()
    os.system(f'sbatch fst_slurm_job.sh')

def submit_cnv_job(bamDirs, genbank, meta=False, outdir=False, name=False):
    ram = 16000
    cpus = 1
    with open(f'cnv_slurm_job.sh', 'w') as sl:
        sl.write("#!/bin/bash\n")
        sl.write("#\n")
        sl.write("#SBATCH --job-name=CNV\n")
        sl.write(f"#SBATCH --mem={ram}\n")
        sl.write(f"#SBATCH --cpus-per-task={cpus}\n")
        # sl.write("#SBATCH --array=1-"+str(array_length)+'\n')
        #sl.write('module load samtools\n')
        # sl.write(f'rsync -a {bamfile} {poly_df} {ref} $TMPDIR\n')
        if outdir:
            if meta:
                if name:
                    sl.write(f'CNV.py -d {",".join(bamDirs)}  -t $TMPDIR -g {genbank} -m {meta} -o {outdir} -n {name}\n')
                else:
                    sl.write(f'CNV.py -d {",".join(bamDirs)}  -t $TMPDIR -g {genbank} -m {meta} -o {outdir}\n')

            else:
                if name:
                    sl.write(f'CNV.py -d {",".join(bamDirs)}  -t $TMPDIR -g {genbank} -o {outdir}  -n {name}\n')
                else:
                    sl.write(f'CNV.py -d {",".join(bamDirs)}  -t $TMPDIR -g {genbank} -o {outdir}\n')

        else:
            if meta:
                if name:
                    sl.write(f'CNV.py -d {",".join(bamDirs)}  -t $TMPDIR -g {genbank} -m {meta} -n {name}\n')
                else:
                    sl.write(f'CNV.py -d {",".join(bamDirs)}  -t $TMPDIR -g {genbank} -m {meta} \n')
            else:
                if name:
                    sl.write(f'CNV.py -d {",".join(bamDirs)}  -t $TMPDIR -g {genbank} -n {name}\n')
                else:
                    sl.write(f'CNV.py -d {",".join(bamDirs)}  -t $TMPDIR -g {genbank}\n')
        sl.close()
    os.system(f'sbatch cnv_slurm_job.sh')


def show_table_for_dNdS(haplo):
    console = Console()
    table = Table(title="Sites for dN/dS calculation in")
    for x in ['# Genes', '# Codons']:
        table.add_column(x)
    table_row =  [str(len(haplo.GeneName.unique())), str(haplo.drop_duplicates(["GeneName",'RefCodon', 'Index_on_cds']).shape[0])]
    table.add_row(*table_row)
    console.print(table)

def Hs_Ht_for_metagenomePop_pair(bam1,bam2, pair_df):
    """The FinalPolyTable.tsv that contains polys for the two samples in comparison. bam1 and bam2 are paths to the respective bam files.
    Returns a dataframe with columns Fst, PopA, PopB"""
    if not os.path.isfile(bam1 + '.bai'):
        os.system(f'samtools index -b {bam1}')
    if not os.path.isfile(bam2 + '.bai'):
        os.system(f'samtools index -b {bam2}')
    bamfile1 = pysam.AlignmentFile(bam1, 'rb')
    bamfile2 = pysam.AlignmentFile(bam2, 'rb')

    pair_df = pair_df.drop_duplicates(['Chrom','Pos','Alt'])
    fst_df = pd.DataFrame()
    if len(pair_df.drop_duplicates("SampleID")['SampleID'].values) == 2: ## make sure no sample is missing from the df (after filtering):
        for i, poly in pair_df.iterrows():
            pos_idx = int(float(poly.Pos) - 1)
            popA = np.array(bamfile1.count_coverage(contig=poly.Chrom, start=pos_idx, stop=pos_idx+1, quality_threshold=20)).flatten() ## count coverage returns an array of size 4 for coverages of A,C,G,T respectively
            popB = np.array(bamfile2.count_coverage(contig=poly.Chrom, start=pos_idx, stop=pos_idx+1, quality_threshold=20)).flatten() ## count coverage returns an array of size 4 for coverages of A,C,G,T respectively

            popM = np.sum([popA,popB], axis=0)
            freqA = pd.Series(popA).value_counts(normalize=True).values
            freqB = pd.Series(popB).value_counts(normalize=True).values
            freqM = pd.Series(popM).value_counts(normalize=True).values
            Ha = 1 - np.sum(np.power(freqA, 2))
            Hb = 1 - np.sum(np.power(freqB, 2))
            Hs = np.mean([Ha, Hb])  ### average Heterozygosity in the subpopulation (average over the two populations)
            Ht = 1 - np.sum(np.power(freqM, 2))  ## average Heterozygosity in the metapopulation (two populations combined)
            if Ht > 0:
                fst = (Ht - Hs) / Ht
            else:
                fst = 0

            fst_df = pd.concat([fst_df,pd.DataFrame({'Fst':[fst], "Chrom":[poly.Chrom], "Pos":[poly.Pos],
                                                 'PopA':[pair_df.drop_duplicates("SampleID")["SampleID"].values[0]], 'PopB':[pair_df.drop_duplicates("SampleID")["SampleID"].values[1]]})],ignore_index=True)

        fst_df_mean = 0 #fst_df.agg(np.mean).to_frame().T.drop(columns="Pos")
        return fst_df_mean, fst_df[['Fst','Chrom','Pos']]

def locate_poly_on_cds(seq: str, poly: str, FivePrimeEnd: int, ThreePrimeEend: int, strand: str):
    """

    :param seq: cds sequence (str)
    :param poly: PolyID
    :param FivePrimeEnd: most 5 prime end coordinate of cds --> start (in Genbank)
    :param ThreePrimeEend: most 3 prime end coordinate of cds --> end (in Genbank)
    :param strand: Forward/plus strand (+) or reverse/minus strand (-) location of cds on the genome
    :return: A dataframe that associates the input polymorphim to the respective codon in the cds sequence
             and adds the mutated codon and the respective AA.

    """
    strand_param = ['+', '-']
    if strand not in strand_param:
        raise ValueError(f"Invalid value '{strand}'. Choose from {strand_param}.")
    if not len(poly.split(":")) == 4:
        raise ValueError(f"Invalid value '{poly}'. Mutations must be in following format: CHROM:POS:REF:ALT")

    chrom, pos, ref, alt = poly.split(":")
    pos = int(pos)
    """ locate position of polys on cds, position, gene start are 1-based index, must be converted to 0-based index
     five-prime-end - pos returns  the index of pos in the cds
     If cds is in - strand, then the placement of the mutation needs to be estimated from the three-prime-end.
     three-prime end - pos returns the index of pos in the cds 
     cds --> + strand
     mutation (pos) --> + strand
     placement --> depends on cds original strand. 

     Given cds='AAATTTCCCGGG' that starts at 10 and CHROM:20:A:C if cds on (+) then  
     mutant= 'AAATTTCCCGCG', if cds on - strand then mutant becomes 'ACATTTCCCGGG'

     First assign polys to codons then do the calculation for pnps for each codon. 
     """

    if strand == '-':
        placement = ThreePrimeEend - pos  ## index on cds
        alt = str(Seq(alt).complement())

        #seq = str(Seq(seq).complement())

    elif strand == "+":

        placement = pos - FivePrimeEnd - 1 ## index on cds
    else:
        print("SS"); sys.exit()
    codon_index_on_cds = int(np.floor(placement / 3))
    base_index_on_codon = int(placement % 3)
    #if len(alt) > 1:

    # print(f'Codon index in CDS: {codon_index_on_cds}')
    # print(f'Base index on codon: {base_index_on_codon}')
    # print(len(CDS(seq).codons))

    # print(f'({strand}) [{FivePrimeEnd} - {ThreePrimeEend}]')
    # print([(i,x) for i,x in enumerate(CDS(seq).codons)])
    # print(f'Length: {len(seq)}')
    # print(f"Placement: {placement}")
    # print(f'Codon index: {codon_index_on_cds}')
    # print(f'Base index: {base_index_on_codon}')
    # print(poly)
    # print(CDS(seq).codons[codon_index_on_cds])
    # print('\n')
    try: ## this exception was added because of two genes (CDS out of range)
        original_codon = CDS(seq).codons[codon_index_on_cds]
        mutant_codon = list(str(original_codon))
        mutant_codon[base_index_on_codon] = alt
        mutant_codon = Codon("".join(mutant_codon))
        #print(f'Original codon: {original_codon}\n'
        #f'Mutant codon: {mutant_codon}')


        return pd.DataFrame({'CodonID': [codon_index_on_cds], 'PolyID': [poly],
                         'RefCodon': [original_codon], 'ObsCodon': [mutant_codon],
                         'RefAA':[original_codon.aminoacid],'ObsAA':[mutant_codon.aminoacid]})
    except:
        print(f'Index out of range {strand} {FivePrimeEnd} {ThreePrimeEend} {FivePrimeEnd + len(seq) -1} {poly}')


def estimate_gene_pnps_stats(codon_df):
    """
    :param codon_df: output from locate_poly_on_cds(), is a dataframe where each row represents a polymorphism in a sample,
                     and associates it with the respective codon on the cds.
                     There might be multiple rows associated with one codon.
                     Columns: CodonID (index), PolyID, RefCodon, ObsCodon, RefAA, ObsAA
    :return: for each codon estimates the weight average from read support of observed Non-synonymous and synonymous polys,
    as well as expected number for each category and pN, pS estimates.
    """
    gene_df_out = pd.DataFrame()
    gene_eN = 0
    gene_eS = 0
    gene_oN = 0
    gene_oS = 0
    for codon, c_group in codon_df.groupby("CodonID"):
        gene_eN += c_group['RefCodon'].values[0].N_sites
        gene_eS += c_group['RefCodon'].values[0].S_sites
        #print(c_group)
        """Here I need a weighted average for all polys (mutated_codons) associated with this codon ID
           
        """
        #total_coverage = c_group['AltCov'].sum()
        #c_group['Freq'] = c_group['AltCov'] / total_coverage ## weigths --> fractions
        ### check if alwsays sum of AltCov < total_coverage
        #assert c_group['AltCov'].sum() <= total_coverage
        ### estimate weigthed average pN from all Non synonymous polys

        gene_oN += c_group.loc[c_group['RefAA'] != c_group['ObsAA']].shape[0]
        gene_oS += c_group.loc[c_group['RefAA'] == c_group['ObsAA']].shape[0]

    if gene_eN > 0:
        pN = gene_oN / gene_eN
    else:
        pN = np.nan
    if gene_eS > 0:
        pS = gene_oS / gene_eS
    else:
        pS = np.nan

    if pS > 0:
        gene_pnps = pN / pS
    else:
        gene_pnps = np.nan
    gene_df_out = pd.concat([gene_df_out, pd.DataFrame({'Gene_eN': [gene_eN], 'Gene_eS': [gene_eS],
                                            'Gene_oN': [gene_oN], 'Gene_oS': [gene_oS],
                                            "Gene_pN": [pN],'Gene_pS':[pS], 'Gene_pNpS':[gene_pnps],
                                            'RefCodon':[c_group['RefCodon'].values[0]]})])
    return gene_df_out


def find_codon_coords(seq, start, idx=0, return_all = False):
    """

    :param seq:
    :param start:
    :param idx:
    :param return_all: If this is true returns a dictionary with keys CodonIndex and values [start, end] for each codon. Else returns only coords for specific idx
    :return:
    """
    out_dict = {}
    the_codons  = CDS(seq).codons
    for i, codon in enumerate(the_codons):
        codon_start = start + (i*3)
        codon_end = codon_start + 2
        out_dict[i] = [codon_start, codon_end]
    if return_all:
        return out_dict
    else:
        return out_dict[idx]