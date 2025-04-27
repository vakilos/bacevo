##!/apps/python3/3.9.0/bin/python3
import os
import sys
from bacevo_fun import *
from popgen_fun import *
console = Console()

parser = argparse.ArgumentParser(description="Creates a master table from breseq files combined with a meta data file. Run this script after ep_start.py.")
parser.add_argument('-d', dest='dir', help='Directory with breseq sample folders (ep_start.py directory structure)', default=os.getcwd())
#parser.add_argument('-r', dest='ref', help='directory to reference file. (fasta, genbank, gff3)', required=True)
parser.add_argument("-m", dest='meta', help='File with meta data for libraries. For invitro: Replicate, Day, Mix, SampleID, for invivo: Mouse, Day, ')
parser.add_argument("-n", dest='name', help="Add a tag to the exported master table.",default=None)
parser.add_argument("-a", dest='anno', help="Path to annotation file (tab-delimited) created by ep_annotable.py.", default='/home/christos/bacevo/reference/Btheta-VPI-5482/GCF_000011065.1/AnnotationTable.txt')
parser.add_argument("-c", dest='coverage', help ='File with genome coverage created with alignemnt investigation function. If not specified will expect files in sample directories to concatenate.', default=None)
parser.add_argument('-x', dest='mode', choices=['invivo','invitro'], default='invivo',help="Choose between invivo or invitro mode, this will affect meta data parsing.")
args = parser.parse_args()

if args.dir[-1] != "/": args.dir += "/"
DataDir = args.dir
ParentDir = Path(DataDir).parent.absolute()
MetaFile = args.meta
AlignmentFile = args.coverage
AnnotationFile = args.anno
NameTag = args.name+"_" if args.name else ""
BioUnit = 'Mouse' if args.mode == "invivo" else 'Replicate'
Condition = 'Cage' if args.mode == 'invivo' else 'Mix'
ancestralFreqThres = 0.9 #.99

os.chdir(ParentDir) ## change to parent directory
console.print(f'Working directory set to {ParentDir}')

if not AlignmentFile:
    concat_alignment_investigation(DataDir).to_csv(NameTag+"Alignment_investigation.csv",index=False) ## make file if not specified..
    AlignmentFile = NameTag+"Alignment_investigation.csv"
breseqFolders = [DataDir+folder for folder in os.listdir(DataDir) if os.path.isdir(DataDir+folder) is True]  ## this is full path

masterDf = pd.DataFrame()
print()
for index, sample in zip(track(range(len((breseqFolders))), description='Parsing VCFs...'), breseqFolders):
    vcfFile = sample+"/output.vcf"
    id = os.path.basename(sample)

    sampleDf = parse_vcf(id=str(id), vcfFile=vcfFile, meta_path=MetaFile, mode=args.mode)
    masterDf = pd.concat([masterDf, sampleDf], ignore_index=True)
    # parentDir = "/".join(args.dir.split("/")[:-2])
# os.chdir(parentDir)
masterDf.SampleID = masterDf.SampleID.astype(str)
masterDf.to_csv(NameTag+"Evopipe_master.txt", sep='\t', na_rep="NaN", index=False)
console.print(f'\nExported [u i]{NameTag+"Evopipe_master.txt"}[/u i] [bold green]SUCCESSFULLY![/bold green]')
masterDf = update_with_ancestral_polys(masterDf, ancestralFreqThres=ancestralFreqThres, bioUnit=BioUnit, condition=Condition, onlyAncestry=True, name=args.name)


##### import file from ep_align_stats.py
cov = pd.read_csv(AlignmentFile)
cov['SampleID'] = cov["SampleID"].astype(str)
masterDf = masterDf.merge(cov, on = ['SampleID',"Chrom"])

#print(masterDf.groupby("SampleID")["CoverageMean"].drop_duplicates())

annoDf = annotate_master_table(masterDf, AnnotationFile)

def assign_synonymy_new(Df):
    console.print('Assigning synonymy...')
    the_synonymy_dict = {}
    the_codonid_dict = {}
    for i, row in Df.drop_duplicates('PolyID').loc[(~Df.translation.isna()) & (Df.cds.str.len() % 3 == 0)].iterrows():
        strand = ""
        start = int(float(row.start))
        end = int(float(row.end))
        pos = int(row.Pos)
        if row.strand == -1.0:
            strand = '-'
        elif row.strand == 1.0:
            strand = '+'
        else:
            print('No DNA strand detected.');sys.exit()
        if row.Type == 'SNP':
            poly_context = locate_poly_on_cds(row.cds, row.PolyID, start, end, strand)
            the_codonid_dict[row.PolyID] = poly_context['CodonID'].item()
            if poly_context['RefAA'].item() == poly_context['ObsAA'].item():
                the_synonymy_dict[row.PolyID] = 'Synonymous'
            else:
                the_synonymy_dict[row.PolyID] = 'Nonsynonymous'


    Df['Synonymy'] = Df['PolyID'].map(the_synonymy_dict).fillna('undefined')
    Df['CodonIndex'] = Df['PolyID'].map(the_codonid_dict).fillna('undefined')
    Df['SampleCodon'] = Df.apply(lambda x: str(CDS(x.cds).codons[int(x.CodonIndex)]) if x.CodonIndex != 'undefined' else "undefined", axis=1)
    return Df


annoDf = assign_synonymy_new(annoDf)
annoDf = add_snp_type(annoDf)
annoDf.drop(columns="index", inplace=True)
annoDf.to_csv(NameTag+"AnnotatedPolyTable.txt", sep='\t', na_rep="", index=False)
console.print(f'\nExported [u i]{NameTag+"AnnotatedPolyTable.txt"}[/u i] [bold green]SUCCESSFULLY![/bold green]')
console.print('\nDONE.')


