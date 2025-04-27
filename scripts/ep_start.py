##!/apps/python3/3.9.0/bin/python3

from bacevo_fun import *

parser = argparse.ArgumentParser(description="""
Relies on a SLURM configuration to run breseq by submitting an array job. 
""")

parser.add_argument('-d', dest='bdir', nargs="+", help='directories with the BAM files, if multiple space seperate them',
                    required=True)
#parser.add_argument('-o', dest='odir', help='directory to copy files to')
parser.add_argument('-f', dest='File', help='file with sample names that need to be selected (True sample names). ONLY for -b modeRun ep_collect_bsf_samples.py to created this file.')
parser.add_argument('-c', dest='root', help='Root directory for the analysis.', default=os.getcwd())
parser.add_argument('-r', dest='ref', help='directory to reference file. (fasta, genbank, gff3)', default="/scratch/christos/bacteroides_evolution/reference/GCF_000011065.1_ASM1106v1_genomic.gbff")
parser.add_argument('-s', dest='subsample', help='Threshold for subsample reads from each fastq file.', default=False)
parser.add_argument('-m', dest='mode',  choices=['p','i'], default='p', help='If p is selected breseq runs in polymorphism mode (metagenomes). If i is selected runs in breseq in isolate mode. ')
parser.add_argument('-a', dest='array', default=None, help='submit different parts of the array job')
parser.add_argument('-e', '--ram',dest='ram', default="20000", help='Memory needed for the job. Default=20000.')
parser.add_argument('-x', dest='dataname', default='', help='Suffix to add to the data folder name.')
parser.add_argument("-l",dest="local", help='Run analysis locally (not with slurm) ', default=False)
args = parser.parse_args()

"""                                    Set up variables                  """

data_name = "data_"+args.dataname
BreseqMode = args.mode
RootDir = args.root
if RootDir[-1] != "/": RootDir += "/"
os.chdir(RootDir)  ### change to root directory...

BamDirList = args.bdir
BamDirList = [BamDir+'/' for BamDir in BamDirList if BamDir.endswith("/") is False]


os.mkdir('samples') if not os.path.isdir(RootDir+'samples') else os.system("rm  "+RootDir+'/samples/*')
OutDir = RootDir+'samples/'



for BamDir in BamDirList: ### more than one directories to look for BAM files
    os.chdir(BamDir)
    for File in glob.glob("*.bam"):
        print(File)
        if os.path.isfile(BamDir + File):
            os.system('ln -s '+ BamDir + File + " " + OutDir + os.path.basename(File)[:5]) ### the ID should be the first 5 charachters
    os.chdir(RootDir)
print(f"\n\nSymbolic links were created in '{OutDir}'.")

make_array_file(RootDir+"samples", outname=RootDir+'slurm_array_file.txt')

if not os.path.isdir(RootDir+data_name): os.mkdir(data_name) ## will create a data folder if it does not exist
if args.local is False:
    submit_array_evopipe_files_from_one_dir(RootDir+'slurm_array_file.txt', args.ref, LockFile=RootDir+data_name+'.lock', subsample=args.subsample, mode=BreseqMode, ram=args.ram)
    array_length = pd.read_csv(RootDir+'slurm_array_file.txt', header=None).shape[0]
    """ IMPORTANT: Submit from data directory as this will be the export directory.. """
    os.chdir(RootDir+data_name)
    os.system("sbatch -a "+args.array+" "+RootDir+"array_breseq_slurm.sh") if args.array else os.system("sbatch -a 1-"+str(array_length)+" "+RootDir+"array_breseq_slurm.sh")

else:
    print("YET NOT IMPLEMENTED")
