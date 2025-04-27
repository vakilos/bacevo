#!/apps/python3/3.9.0/bin/python3
import glob
import os
import sys
import pandas as pd
from bacevo_fun import *
import argparse
#from string import maketrans


parser=argparse.ArgumentParser(description="Estimates number of mapping statistics for a bam file, number of mapped reads etc.")
parser.add_argument('-i', dest='bamdir', help ='Full path to input file', default=None)
#parser.add_argument('-g', dest='genome', help = "Chromosome/ contig / genome name to look for the region of interest.", default = "NC_004663")
#parser.add_argument('-f',dest='form', help='Format of input file', choices=['sam','bam'], default="bam")
parser.add_argument('-f',dest='folder',help='Run in folder mode..Each sample is under a folder. Use parent directory to all folders.', default=None)
parser.add_argument('-n', dest='name', help='Project name to be added in the exported file name (for directories).', default=None)
parser.add_argument('-c', dest='concatenate',help='Use this argument to concatenate multiple AlignmentInvestigation files. Run from parent folder.', action='store_true')
parser.add_argument('-w', dest='window',help='length of window scanning the genome.', default=100, type=int)
parser.add_argument('-o', dest='overlap', help='Overalap of windows scanning the genome', default=0, type=int)
args=parser.parse_args()


console = Console(color_system='windows')

if args.name:
    proj_name = args.name+"_"
else:
    proj_name = ""

if args.concatenate:
    os.chdir(os.getcwd())
    filelist = [pd.read_csv(x) for x in glob.glob("*/Coverage_by_Position.csv", recursive=True)]

    pd.concat(filelist).to_csv('Coverage_By_Position_All_Files.csv', index=False)
    sys.exit()
currentDir=os.getcwd()
#bamdir=os.path.dirname(args.inDir)
#bamdir = "/media/christos/ssd/work/bacevo/data_invitro_isolates/69383/breseq/data/reference.bam"

if args.folder != None:
    os.chdir(args.folder)
    fdf = pd.DataFrame()
    for samplefolder in [x for x in os.listdir(args.folder) if os.path.isdir(x)]:
        console.print(f"\n\n\t\t\t[bold blue]{samplefolder}[/bold blue]")
        #print(color.BOLD+color.CYAN+"\n\t\tInvestigating {}\n".format(samplefolder)+color.END)
        os.chdir(samplefolder)
        #os.system("samtools view -@ 2 -bq 20 reference.bam > unique.bam")
        this_alignment = investigate_alignment("reference.bam", SampleID=samplefolder, window_size=args.window, overlap=args.overlap)
        this_alignment.to_csv(f"{samplefolder}_Alignment_investigation.csv")
        fdf = pd.concat([fdf, this_alignment])
        #os.system("rm unique.bam*")
        os.chdir(args.folder)
    fdf.to_csv("../"+proj_name+"Alignment_investigation.csv", index=False)
elif args.bamdir:
    os.chdir(os.path.dirname(args.bamdir))
    #os.system("samtools view -@ 2 -bq 20 reference.bam > unique.bam")
    odf = investigate_alignment("reference.bam", SampleID=os.path.basename(os.path.dirname(args.bamdir)), window_size=args.window,overlap=args.ovelap)
    odf.to_csv(args.bamdir.replace(".bam","_stats.csv"), index=False)
    #os.system("rm unique.bam*")