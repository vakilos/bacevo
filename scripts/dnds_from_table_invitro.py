import glob
import itertools
import os.path
import sys
import pandas as pd

from popgen_fun import *
parser = argparse.ArgumentParser(description="By default exports files in the current working directory")
parser.add_argument('-d', dest='bamdirs', nargs='+',help='Directory(ies) containing sample folders with BAM files from the alignment. (or fasta files for genome mode)', default=False, required=False)
parser.add_argument('-p', dest='polytable', help='FinalPolyTable.tsv from bacevo pipeline.', required=True)
parser.add_argument('-n', dest='name', help='Name of the exported files for isolates population (genotypes for each individual) ',default=False)
parser.add_argument('-m', "--mq", dest='mq', help='Mapping quality to filter reads in the alignment when estimating dNdS, Default=20', default=20, type=int)
parser.add_argument('-s', "--support", dest='support', help='Reads supporting the alternative allele. Default=4', default=2)
parser.add_argument('-q', "--qc", dest='qc', help='Base quality filter. Default=30', default=20, type=int)
parser.add_argument('-o', dest='outdir', help="Directory to export files to. Default = creates dnds directory in CWD.", default=os.getcwd()+'/dnds')
parser.add_argument('-x', dest='mode',help='metagenome or genome mode. DEFAULT = metagenome', choices=['metagenome','genome'], default='metagenome')
parser.add_argument('-f', dest="BamFile", help='Path to one bam file instead of directory', default=None,)
parser.add_argument('-k', dest='maskFile', help='Path to a file with polys to be masked.', default=None)
parser.add_argument('-w', dest='over', help='Overwrite existing data. Default is FALSE', action='store_true', default=False)
parser.add_argument('-r', dest='ref', help='Reference fasta', required=True)

args = parser.parse_args()

os.chdir(os.getcwd())
cur = os.getcwd()
export_dir = args.outdir
if not os.path.isdir(export_dir):
    os.mkdir(export_dir)
console = Console(color_system='windows')
def find_bam_files(Dirs):
    Current = os.getcwd()
    bams =[]
    for Dir in Dirs:
        os.chdir(Dir)
        bams += [x for x in glob.glob(f"{Dir}/*/reference.bam")]
        os.chdir(Current)
    return bams
def find_fasta_files(Dirs):
    Current = os.getcwd()
    fastas = []
    for Dir in Dirs:
        os.chdir(Dir)
        fastas += [x for x in glob.glob(f'{Dir}/*/reference.fasta')]
        os.chdir(Current)
    return fastas

overwrite = args.over


use_cols = ['SampleID','Replicate','Mix','PolyID', 'Chrom','Pos','Type','Ref','Alt', 'AltCov','TotalCov',
                               'GeneName','oldGeneName','cds', 'start','end','strand', 'Ancestry', "Coverage_mean", "CodonIndex", 'translation']
if args.mode == 'genome':
    use_cols.remove("Dataset")
poly_df = pd.read_csv(args.polytable, sep='\t',
                      usecols= use_cols, dtype={'SampleID':str})


poly_df = poly_df.loc[(poly_df.Type == "SNP") & (poly_df.translation != "Nan") & (poly_df.cds != "Nan") & (poly_df.Ancestry == 'evolved')].sort_values("Pos")

if args.maskFile:
    masked = pd.read_csv(args.maskFile, usecols=["PolyID"]).drop_duplicates()
    poly_df = poly_df.loc[~poly_df.PolyID.isin(masked.PolyID.values)]

poly_df = poly_df.drop_duplicates(['SampleID','PolyID'])

poly_df[['start','end']] = poly_df[['start','end']].astype(float).astype(int)



#ref = fasta_to_df('/home/christos/bacevo/reference/Btheta-VPI-5482/GCF_000011065.1/GCF_000011065.1_ASM1106v1_genomic.fna')


ref = fasta_to_df(args.ref)
fasta = [x for x in SeqIO.parse(args.ref, "fasta")]
for x in fasta:
    x.id = x.id.strip(".1")
SeqIO.write(fasta, "ref.fa", 'fasta')

if args.mode == 'metagenome':

    if not args.BamFile:

        all_missing_polys = []

        bams = find_bam_files(args.bamdirs)

        for i,bam in enumerate(bams):
            if overwrite:
                os.system(f'rm {bam.replace(".bam","_MD.bam")} ')
                os.system(f'rm {bam.replace(".bam","_MD.bam.bai")} ')
                ## remove older index files

            if not os.path.isfile(f'{bam.replace(".bam","_MD.bam")}'):
                console.print("\nAdding MD tag...")
                os.system(f'samtools calmd -@ 4 -b {bam} ref.fa > {bam.replace(".bam","_MD.bam")}')

            sampleid = re.match(r'.+/(\d+)(?=/reference.bam)', bam).group(1)
            print(f'\t\t\t\t{sampleid} ({i+1}/{len(bams)})')
            #if not os.path.isfile(export_dir + "/" + sampleid + '_codon_dNdS.csv'):
            if not os.path.isfile(export_dir + "/" + sampleid + '_haplotable.csv'):
                haplotype_df, missing_polys = codon_haplotypes_from_bam(bam.replace(".bam","_MD.bam"), sampleid, poly_df.loc[poly_df.SampleID == sampleid], imputation=True)
                haplotype_df.to_csv(f'{export_dir}/{sampleid}_haplotable.csv', index=False)
                for poly in missing_polys:
                    if poly not in all_missing_polys:
                        all_missing_polys.append(poly)
            else:
                haplotype_df = pd.DataFrame()
                try:
                    haplotype_df = pd.read_csv(f'{export_dir}/{sampleid}_haplotable.csv', dtype={'SampleID':str})
                except:
                    print(f'Empty file: {sampleid}_haplotable.csv')
            codon_df, gene_df = pd.DataFrame(), pd.DataFrame()
            if not haplotype_df.empty:
                codon_df, gene_df = calculate_dNdS_from_HaploTable(haplotype_df, codon_read_support = args.support)
            if not codon_df.empty and not gene_df.empty:
                codon_df.merge(poly_df.drop_duplicates('SampleID')[['SampleID', 'Replicate', 'Mix']], on='SampleID',
                               how='left').to_csv(export_dir + "/" + sampleid + "_codon_dNdS.csv", index=None)
                assert gene_df.drop_duplicates(['GeneName','SampleID']).shape[0] == gene_df.shape[0]
                gene_df.merge(poly_df.drop_duplicates('SampleID')[['SampleID', 'Replicate', 'Mix']], on='SampleID',
                    how='left').to_csv(export_dir + "/" + sampleid + '_gene_dNdS.csv', index=None)
        print('\n\n')
        print(f"\n\n{len(all_missing_polys)} polys were not detected in the BAM files of this project (compared to Breseq calls)")
    else:
        bam = args.BamFile
        sampleid = re.match(r'.+/(\d+)(?=/reference_MD.bam)', bam).group(1)
        if not os.path.isfile(export_dir + "/" + sampleid + '_haplotable.csv'):
            haplotype_df, missing_polys = codon_haplotypes_from_bam(bam, sampleid,
                                                                    poly_df.loc[poly_df.SampleID == sampleid])
            haplotype_df.to_csv(f'{export_dir}/{sampleid}_haplotable.csv', index=False)
        else:
            haplotype_df = pd.DataFrame()
            try:
                haplotype_df = pd.read_csv(f'{export_dir}/{sampleid}_haplotable.csv', dtype={'SampleID': str})
            except:
                print(f'Empty file: {sampleid}_haplotable.csv')

        codon_df, gene_df = pd.DataFrame(), pd.DataFrame()
        if not haplotype_df.empty:
            codon_df, gene_df = calculate_dNdS_from_HaploTable(haplotype_df)
        if not codon_df.empty and not gene_df.empty:
            codon_df.merge(poly_df.drop_duplicates('SampleID')[['SampleID', 'Replicate', 'Mix']], on='SampleID',
                           how='left').to_csv(export_dir + "/" + sampleid + "_codon_dNdS.csv", index=None)
            assert gene_df.drop_duplicates(['GeneName','SampleID']).shape[0] == gene_df.shape[0]
            gene_df.merge(poly_df.drop_duplicates('SampleID')[['SampleID', 'Replicate', 'Mix']], on='SampleID',
                how='left').to_csv(export_dir + "/" + sampleid + '_gene_dNdS.csv', index=None)
        print(f"{sampleid}: Nothing to export!")








if args.mode == 'genome':
    fasta_files = find_fasta_files(args.bamdirs)
    haplotype_df = codon_haplotypes_from_genomes_OLD(fasta_files, poly_df)
    if not haplotype_df.empty:
        codon_df, gene_df = calculate_dNdS_from_HaploTable(haplotype_df)
        if not codon_df.empty and not gene_df.empty:
            codon_df.merge(poly_df.drop_duplicates('SampleID')[['SampleID', 'Mouse', 'Dataset','SampleType']], on='SampleID',
                               how='left').to_csv(export_dir + "/" + sampleid + "_codon_dNdS.csv", index=None)
            gene_df.merge(poly_df.drop_duplicates('SampleID')[['SampleID', 'Mouse', 'Dataset','SampleType']], on='SampleID',
                    how='left').to_csv(export_dir + "/" + sampleid + '_gene_dNdS.csv', index=None)