##!/apps/python3/3.9.0/bin/python3
import logging

from rich import print as rprint
from rich.console import Console
from rich.table import Table
from rich.progress import track
from rich.panel import  Panel
from rich.progress import Progress
from rich.tree import Tree
import sys,re,os, subprocess
import argparse
import pandas as pd
import numpy as np
import glob
#import pysam
import scipy.spatial.distance
from Bio.Seq import Seq
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from Bio import SeqIO
import itertools
import scipy.stats
import pysam
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

def chunks_h(l, n):
    """Yield successive n-sized chunks from l."""
    new_obj = []
    for i in range(0, len(l), n):
        new_obj.append(l[i:i + n])
    return new_obj
def chunks_v(l, n):
    """Yield successive n-sized chunks from l."""
    new_obj = []
    for i in range(int(len(l)/n)):
        new_obj.append(l[i::int(len(l)/n)])
    return new_obj
### this is needed for importing pysamstats on the cluster (it is not sudo installed)
#pysamstats_file_path = '/home/user/zioutis/pysamstats/'
#pysam_file_path = '/home/user/zioutis/pysam'
#sys.path.append(pysamstats_file_path)
#sys.path.append(pysam_file_path)
#import pysamstats



def reassign_names_vertically(df):
    """Df must have columns for Row (A,B,C) and Column (1,2,3,4) as well SampleID column for the name of each sample """

    Exp = re.compile(r'BSF_\d{4}_\w{9}_\d_([\w\d_-]{1,})_S(\d{5})\.bam$')  ## extract file label and BSF ID from file name
    well_numbers = np.arange(1,97)
    well_numbers_h = chunks_h(well_numbers, 12)
    well_numbers_v = chunks_v(well_numbers, 12)
    well_numbers_rows = [x for y in well_numbers_h for x in y]
    well_numbers_columns = [x for y in well_numbers_v for x in y]
    assign = pd.DataFrame()
    assignment_rows = {}
    for row,letter in zip(well_numbers_h, list("ABCDEFGH")):
        for index,column in enumerate(row):
            assign = assign.append({'Well':letter+str(index+1), 'r_order':int(column)}, ignore_index=True)
            assignment_rows[column] = letter+str(index+1)
    assign.set_index('r_order',inplace=True)
    wd = pd.DataFrame({'r_order':well_numbers_rows,"c_order":well_numbers_columns}).set_index('r_order')
    wd = wd.merge(assign, on='r_order')
    wd.reset_index(inplace=True)
    wd['c_order_well'] = wd.c_order.map(assignment_rows)
    df['Well'] = df['Row'] + df['Column'].astype(str)
    df = df.merge(wd, on='Well')
    df['TrueID'] = [df.loc[(df['Well'] == row['c_order_well']) & (df['Plate'] == row['Plate'])]['SampleID'].values[0] for index, row in df.iterrows()]
    #df['ResponsiblePerson'] = [df.loc[(df['Well'] == row['c_order_well']) & (df['Plate'] == row['Plate'])]['Responsible_person'].values[0] for index, row in df.iterrows()]
    filesToCopy = pd.read_csv('bsf_bam_filenames.txt', header=None)[0].tolist()  ### read file names from file not from directory
    sampleNamesToBam = {re.match(Exp, BAM).group(1): BAM for BAM in filesToCopy}
    df["bamFileName"] = df["SampleID"].map(sampleNamesToBam)
    df["TrueFileName"] = df['TrueID'].map(sampleNamesToBam)
    df['NewFileName'] = df['TrueFileName'].apply(lambda x: re.match(Exp, x).group(2) if x != np.nan else "empty")
    df.to_csv('all_samples_reassignments.csv')
    return  df

def investigate_alignment(bamfile, SampleID, window_size=100, overlap=0):
    console = Console(color_system='windows')
    os.system("samtools index -b "+bamfile)
    bamfile = pysam.AlignmentFile(bamfile, 'rb')
    unmapped = bamfile.unmapped  # WOULDN'T WORK !!
    mapped = bamfile.mapped
    all_reads = float(unmapped + mapped)
    chroms = list(bamfile.references)

    if all_reads == 0:
        print("No alignments found in the bam file.")
        sys.exit(0)

    ### sliding or not window function
    def window(lista, w, o):
        """Iterates over a list l and calculates a value for a window of width w and overlap o"""
        if o == 0:
            win_list = [sum(lista[i:i + w]) / float(w) for i in range(0, len(lista), w)]
        else:
            win_list = [sum(lista[i:i + w]) / float(w) for i in range(0, len(lista), w - o)]
        return win_list

    #### Estimate coverage all over the genome
    # samtools =subprocess.Popen("samtools depth -r "+args.genome+" -a "+args.inDir,shell=True, stdout=subprocess.PIPE)
    # shell_out=samtools.communicate()[0].split()
    ### Estimating reads stats and Nucleotide Identity Percentage distribution




    ### if Coverage_by_Position.csv exist skip this
    console.print(f'Estimating genome coverage...')
    if not os.path.isfile(f'{os.getcwd()}/Coverage_by_Position.csv'):
        cpDf = pd.DataFrame({"Chrom": [], "Pos": [], "Coverage": []})
        for chrom in chroms:
            ## this depends on pysamstats library
            #pos, cov = window([x[1] for x in pysamstats.load_coverage(bamfile, chrom=chrom)], 100,0), window([x[2] for x in pysamstats.load_coverage(bamfile, chrom=chrom)],100,0) ### this will subsample in a window fashion
            ## this does not
            coverage = np.sum(np.array(bamfile.count_coverage(contig=chrom, quality_threshold=20, read_callback="nofilter")),axis=0).tolist() ## count_coverage function returns an 4 length array count for each base for each position - need to sum over this
            pos, cov = window([i+1 for i in range(len(coverage))], window_size, overlap), \
                       window(coverage, window_size, overlap)  ### this will subsample in a window fashion

            cpDf = pd.concat([cpDf, pd.DataFrame({"SampleID":[SampleID for x in pos],"Chrom":[chrom for x in pos],'Pos':pos,'Coverage':cov})])
        cpDf.to_csv('Coverage_by_Position.csv',index=False)
    else:
        cpDf = pd.read_csv("Coverage_by_Position.csv")
    cDf = cpDf.groupby("Chrom").agg({"Coverage":["mean",'std']})
    cDf.columns = cDf.columns.map('_'.join).str.strip('_')
    cDf.reset_index(inplace=True)
    console.print('Parsing alignment statistics...')
    odf = pd.DataFrame()
    for chrom in chroms:
        chrom_mapped = 0
        chrom_reads = 0
        mq_reads = 0
        secondary = 0
        duplicate = 0
        proper_pair = 0
        paired = 0
        reverse = 0
        qc_fail = 0
        # nucl_ide=[]
        # read_pos=[]
        NIP_on_genome = {}
        pairs_distance = {}
        mapping_qualities = []
        for _,read in zip(track(range(sum(1 for _ in bamfile.fetch(chrom))), description=f'{chrom}'),bamfile.fetch(chrom)):
            mapping_qualities.append(float(read.mapping_quality))
            chrom_reads += 1
            if read.is_mapped:
                chrom_mapped += 1
            if read.mapping_quality > 20:
                mq_reads += 1
            if read.is_secondary:
                secondary += 1
            if read.is_duplicate:
                duplicate += 1
            if read.is_proper_pair:
                proper_pair += 1
            if read.is_reverse:
                reverse += 1
            if read.is_paired:
                paired += 1
            if read.is_qcfail:
                qc_fail += 1
            # if read.reference_start not in pairs_distance:
            #    pairs_distance[read.reference_start] = abs(read.reference_start - read.next_reference_start)
            read_len = 0
            nui = 0
            NPI = 0

            # if read.cigartuples is not None:
            #     for t in read.cigartuples:
            #         if t[0]== 0:
            #             nui+=float(t[1])
            #         read_len+=float(t[1])
            #
            #     NPI=nui/read_len
            #     #nucl_ide.append(NPI)
            #     mid_pos=0
            #     mapping=read.get_aligned_pairs()
            #     for i in range(len(mapping)):
            #         if mapping[i][1] is not None and mapping[len(mapping)-1-i][1] is not None:
            #             mid_pos+=(int(mapping[i][1])+int(mapping[len(mapping)-1-i][1]))/2
            #             break
            #         else:
            #             #print("Deletion or insertion")
            #             continue

            # NIP_on_genome[mid_pos] = NPI

        # df = pd.DataFrame(np.array(mapping_qualities), columns=["MQ"])
        # ax = df.plot(kind='density')
        # lines, labels = ax.get_legend_handles_labels()
        # ax.legend(lines,labels)
        # plt.show()

        # print("Number of reads: {}".format(int(all_reads)))
        # print("% of mapped reads: {}".format(round(100 * mapped / all_reads, 3)))
        # print("% of unmapped reads: {}".format(round(100 * unmapped / all_reads, 3)))
        # print("% with a pair: {}".format(round(100 * paired / all_reads, 3)))
        # print("% of proper pairs: {}".format(round(100 * proper_pair / all_reads, 3)))
        # print("% derive from reverse strand: {}".format(round(100 * reverse / all_reads, 3)))
        # print("% QC filter fail: {}".format(round(100 * qc_fail / all_reads, 3)))
        # print("% of reads with MQ > 20: {}".format(round(100 * mq_reads / all_reads, 3)))
        # print("% secondary alignments: {}".format(round(100 * secondary / all_reads, 3)))
        # print("% of duplicates: {}".format(round(100 * duplicate / all_reads, 3)))

        tdf = pd.DataFrame({"SampleID":[SampleID], "ReadCount":[mapped], "Chrom":[chrom],"ChromReads": [chrom_reads], "Mapped": [chrom_mapped], "Mapped_prc": [round(chrom_mapped / chrom_reads, 3)],
                            "Paired": [paired], "Paired_prc": [round(paired / chrom_reads, 3)],
                            "ReverseStrand": [reverse], "ReverseStrand_prc": [round(reverse / chrom_reads, 3)],
                            "MQpass": [mq_reads], "MQpass_prc": [round(mq_reads / chrom_reads, 3)],
                            "SecondaryAligned": [secondary], "SecondaryAligned_prc": [round(secondary / chrom_reads, 3)],
                            "Duplicates": [duplicate], "Duplicates_prc": [round(duplicate / chrom_reads, 3)]})
        odf = pd.concat([odf,tdf])
    odf = odf.merge(cDf, on='Chrom')
    return odf



def estimate_coverage(bam):
    os.system("samtools index "+bam)
    parsedBam = pysam.AlignmentFile(bam,"rb")
    chroms = list(parsedBam.references)
    df = pd.DataFrame({"SampleID":[ID for chrom in chroms],"Coverage":[np.mean(np.sum(np.array(parsedBam.count_coverage(contig=chrom, quality_threshold=20, read_callback="all")), axis=0)) for chrom in chroms], "Chrom":[x for x in chroms]})
    return df


def bam_to_fastq(dir):
    os.chdir(dir)
    for file in glob.glob("*.bam"):
        os.system("time samtols sort -n "+file+'-o '+file.replace('.bam','_sorted.bam'))
        os.system("time samtools fastq -@ 4 -1 "+file.replace(".bam","_1.fastq")+" -2 "+file.replace(".bam","_2.fastq")+" "+file.replace('.bam','_sorted.bam'))

def make_sample_folders_for_bam_files(dir, sublist):
    os.chdir(dir)
    FileList = [bam for bam in glob.glob("*.bam") if bam in sublist]
    #ID = re.compile(r'(\d+)_.+')
    #identifiers = [re.match(ID,x).group(1) for x in FileList]
    identifiers = [x.replace(".bam","") for x in FileList]
    for i in set(identifiers):
        os.mkdir(i)
        Files = [x for x in FileList if re.match(i,x)]
        for File in Files:
            os.system("mv "+File+" "+i)

def make_sample_folders_for_fastqgz_files(dir, regex=""):
    if regex != "":
        expr = regex
    else:
        expr = r'(\d+)_.+'
    os.chdir(dir)
    FileList = glob.glob("*.fastq.gz")
    ID = re.compile(expr)
    identifiers = [re.search(ID,x).group(0) for x in FileList]
    #identifiers = [x.replace(".fastq.gz","") for x in FileList]
    for i in set(identifiers):
        os.mkdir(i)
        Files = [x for x in FileList if re.search(i,x)]
        for File in Files:
            os.system("mv "+File+" "+i)

def fastp_folder_version_miseq(dir):
    """ Runs fastp for samples in organized in folders."""
    os.chdir(dir)
    root_dir = os.getcwd()
    ID = re.compile(r'([^\.]+)\.(\d).+')
    #FolderList= os.listdir()
    for folder in os.listdir():
        os.chdir(root_dir+"/"+folder)
        files = [file for file in glob.glob("*.[1,2].fastq.gz")]
        New_names =[re.match(ID, file).group(1)+"_"+re.match(ID, file).group(2)+"_paired.fastq.gz" for file in files]
        os.system("fastp -i "+files[0]+" -I "+files[1]+" -o "+New_names[0]+" -O "+New_names[1]+" -q 20 --detect_adapter_for_pe -l 80 --cut_tail --cut_tail_window_size 1  --trim_tail2 40")   ### automatically detects adapter sequences, filters reads with base quality > 20 and > 80 bp and clips the trailing 3' , windows size = 1 if below quality threshold (20).
def fastp_folder_version_hiseq(dir):
    """ Runs fastp for samples in organized in folders."""
    os.chdir(dir)
    root_dir = os.getcwd()
    ID = re.compile(r'([^\.]+)\.(\d).+')
    #FolderList= os.listdir()
    for folder in os.listdir():
        os.chdir(root_dir+"/"+folder)
        files = [file for file in glob.glob("*.[1,2].fastq.gz")]
        New_names =[re.match(ID, file).group(1)+"_"+re.match(ID, file).group(2)+"_paired.fastq.gz" for file in files]
        os.system("fastp -i "+files[0]+" -I "+files[1]+" -o "+New_names[0]+" -O "+New_names[1]+" -q 20 --detect_adapter_for_pe -l 80 --cut_tail --cut_tail_window_size 1")

def fastp_sample_version(dir):
    "Runs fastp for an individual sample"
    os.chdir(dir)
    ID = re.compile(r'([^\.]+)\.(\d).+')

def bowtie2_folder_version(dir, ref_dir):
    os.chdir(dir)
    root_dir = os.getcwd()

    for folder in os.listdir():
        os.chdir(root_dir+"/"+folder)
        print(os.getcwd())
        file1, file2 = [file for file in glob.glob("*[1,2]_paired.fastq.gz")]

        output_name =file1.replace("_1_paired.fastq.gz",'.bam')
        os.system("bowtie2 -x "+ref_dir+"/theta_ref -1 "+file1+" -2 "+file2+" | samtools view -bS -q 20 -hb -F 3844  -  | samtools sort  - > "+output_name) ## I have ran bowtie2 build for indexing reference genome.

def breseq_isolate_folder_version(dir, ref_gbk):
    os.chdir(dir)
    root_dir = os.getcwd()
    for folder in os.listdir():
        os.chdir(root_dir+"/"+folder)
        if not os.path.exists(root_dir+"/"+folder+'/breseq'):
            os.mkdir('breseq')
        else:
            os.system("rm -rf "+root_dir+"/"+folder+'/breseq/*') ## I need to remove breseq output as it does not overwrite files...
        file1, file2 = [file for file in glob.glob("*[1,2]_paired.fastq.gz")]
        output_name =file1.replace("_1_paired.fastq.gz",'')
        os.system("breseq -j 2  --base-quality-cutoff 20  -o breseq -n "+output_name+" -r "+ref_gbk+" "+file1+" "+file2)

def breseq_metagenome_version(dir, ref_gbk):
    os.chdir(dir)
    root_dir = os.getcwd()
    if not os.path.exists(root_dir+"/breseq"):
        os.mkdir('breseq')
    else:
        os.system("rm -rf "+root_dir+"/breseq/*") ## I need to remove breseq output as it does not overwrite files...
    file1, file2 = [file for file in glob.glob("*[1,2]_paired.fastq.gz")]
    output_name =file1.replace("_1_paired.fastq.gz",'')
    os.system("time breseq -j 4 -p --base-quality-cutoff 20  -o breseq -n "+output_name+" -r "+ref_gbk+" "+file1+" "+file2)
def bcftools_isolate_folder_version(dir,ref):
    os.chdir(dir)
    root_dir = os.getcwd()
    for folder in os.listdir():
        os.chdir(root_dir+"/"+folder)
        file1 = glob.glob("*.bam")[0]
        output_name = file1.replace(".bam",".vcf")
        os.system("bcftools mpileup -Ou -Q 20 -f "+ref+" "+file1+" | bcftools call - -o "+output_name+" --ploidy 1 -cv -Ob") ### concensus caller




def make_array_file(Directory, outname='slurm_array_file.txt'):
    currentDir = os.getcwd()
    os.chdir(Directory)
    if os.path.isfile(outname): os.system('rm '+outname) ## otherwise will append to the existing file
    os.system("ls -d $PWD/* | cat  > "+outname)
    os.chdir(currentDir)
    print(f"\n\nArray file has been created.")

def submit_array_breseq_assembly_isolates(array_file, ref, LockFile='foo.lock', subsample=False, ram="10000"):
    """Submits slurm scripts for running evopipe pipeline for each symbolic link for BAM file in the folders of the input directory.
    . Submit array jobs.Files will be exported to the input directory"""
    print("\n\nPreparing slurm array script...")
    print(f'::Breseq will run in isolate mode...')
    cpus = "1"
    with open("array_breseq_isolate_assebly.sh", 'w') as sl:
        sl.write("#!/bin/bash\n")
        sl.write("#\n")
        sl.write("#SBATCH --job-name=array_breseq\n")
        sl.write("#SBATCH --mem=" + ram + "\n")
        sl.write("#SBATCH --cpus-per-task=" + cpus + "\n")
        # sl.write("#SBATCH --array=1-"+str(array_length)+'\n')
        sl.write('module load samtools\n')
        sl.write('module load fastp\n')
        sl.write('module load seqtk\n')
        sl.write('module load breseq\n')
        sl.write('name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ' + array_file + ')\n')
        sl.write('BASE=${name##*/}\n')
        sl.write('\n')
        sl.write('link=`readlink -f $name`\n')
        sl.write('BASELINK=${link##*/}\n')
        sl.write("LOCKFILE=" + LockFile + "\n")
        sl.write("module load breseq\n")
        sl.write('set -e\n')
        sl.write('(\n')
        sl.write(' flock -w 3600 200\n')
        sl.write(" rsync -av " + ref+ " $link ${TMPDIR}\n")
        sl.write(") 200>$LOCKFILE\n")
        sl.write("samtools sort -@ " + cpus + " -n ${TMPDIR}/$BASELINK -o ${TMPDIR}/sorted.bam\n")
        sl.write("samtools fastq -@ " + cpus + " -1 ${TMPDIR}/foo_1.fq -2 ${TMPDIR}/foo_2.fq ${TMPDIR}/sorted.bam\n")
        sl.write(
            'fastp -i ${TMPDIR}/foo_1.fq -I ${TMPDIR}/foo_2.fq -o ${TMPDIR}/paired_1.fq -O ${TMPDIR}/paired_2.fq -q 20 --detect_adapter_for_pe -l 80 --cut_tail --cut_tail_window_size 1\n')
        if subsample:
            sl.write('seqtk sample -s 100 ${TMPDIR}/paired_1.fq ' + str(subsample) + ' > ${TMPDIR}/short_paired_1.fq\n')
            sl.write('seqtk sample -s 100 ${TMPDIR}/paired_2.fq ' + str(subsample) + ' > ${TMPDIR}/short_paired_2.fq\n')
            breIn1 = "short_paired_1.fq"
            breIn2 = "short_paired_2.fq"
        else:
            breIn1 = "paired_1.fq"
            breIn2 = "paired_2.fq"
        ### add rarefraction step here


        sl.write("breseq -j " + cpus + " --polymorphism-minimum-total-coverage-each-strand 2 "
                                                                       "--polymorphism-frequency-cutoff 0  --brief-html-output --no-javascript"
                                                                       "-o ${TMPDIR}"
                                                                       " -r ${TMPDIR}/" + os.path.basename(
                    ref) + "  ${TMPDIR}/" + breIn1 + " ${TMPDIR}/" + breIn2 + "\n")  ## breseq version installed on my server's home dir

        sl.write("mv ${TMPDIR}/output ${TMPDIR}/data\n")
        sl.write("mv ${TMPDIR}/data ${TMPDIR}/${BASE/.bam/}\n")
        sl.write(
            "ep_align_stats.py -i ${TMPDIR}/${BASE/.bam/}/reference.bam\n")  ## add investigate alignment here
        sl.write("cp -r ${TMPDIR}/${BASE/.bam/} $PWD\n")
        sl.write('rm -rf ${TMPDIR}\n')
    sl.close()
def submit_array_evopipe_files_from_one_dir(array_file,  ref_gbk, LockFile='foo.lock', subsample=False, mode='p', ram="20000"):
    """Submits slurm scripts for running evopipe pipeline for each symbolic link for BAM file in the input directory.
    The name is the TrueSampleID whereas the link points to the file that corresponds. Submit array job.
    Files will be exported to the input directory"""
    print("\n\nPreparing slurm array script...")
    print(f'::Breseq will run in "{mode}" mode...')
    cpus = "1"
    with open("array_breseq_slurm.sh",'w') as sl:
        sl.write("#!/bin/bash\n")
        sl.write("#\n")
        sl.write("#SBATCH --job-name=array_breseq\n")
        sl.write("#SBATCH --mem="+ram+"\n")
        sl.write("#SBATCH --cpus-per-task="+cpus+"\n")
        #sl.write("#SBATCH --array=1-"+str(array_length)+'\n')
        sl.write('module load samtools\n')
        sl.write('module load fastp\n')
        sl.write('module load seqtk\n')
        sl.write('module load breseq\n')
        sl.write('name=$(sed -n "$SLURM_ARRAY_TASK_ID"p '+array_file+')\n')
        sl.write('BASE=${name##*/}\n')
        sl.write('\n')
        sl.write('link=`readlink -f $name`\n')
        sl.write('BASELINK=${link##*/}\n')
        sl.write("LOCKFILE="+LockFile+"\n")
        sl.write("module load breseq\n")
        sl.write('set -e\n')
        sl.write('(\n')
        sl.write(' flock -w 3600 200\n')
        sl.write(" rsync -av "+ref_gbk+" $link ${TMPDIR}\n")
        sl.write(") 200>$LOCKFILE\n")
        sl.write("samtools sort -@ "+cpus+" -n ${TMPDIR}/$BASELINK -o ${TMPDIR}/sorted.bam\n")
        sl.write("samtools fastq -@ "+cpus+" -1 ${TMPDIR}/foo_1.fq -2 ${TMPDIR}/foo_2.fq ${TMPDIR}/sorted.bam\n")
        sl.write('fastp -i ${TMPDIR}/foo_1.fq -I ${TMPDIR}/foo_2.fq -o ${TMPDIR}/paired_1.fq -O ${TMPDIR}/paired_2.fq -q 20 --detect_adapter_for_pe -l 80 --cut_tail --cut_tail_window_size 1\n')
        if subsample:
            sl.write('seqtk sample -s 100 ${TMPDIR}/paired_1.fq '+str(subsample)+' > ${TMPDIR}/short_paired_1.fq\n')
            sl.write('seqtk sample -s 100 ${TMPDIR}/paired_2.fq '+str(subsample)+' > ${TMPDIR}/short_paired_2.fq\n')
            breIn1="short_paired_1.fq"
            breIn2="short_paired_2.fq"
        else:
            breIn1 = "paired_1.fq"
            breIn2 = "paired_2.fq"
        ### add rarefraction step here
        if mode == 'p':
            sl.write("~/breseq-0.35.5-Linux-x86_64/bin/breseq -j "+cpus+" --polymorphism-minimum-total-coverage-each-strand 2 "
                                                                        "--polymorphism-frequency-cutoff 0 -p --brief-html-output --no-javascript "
                                                                        "-o ${TMPDIR} "
                                                                        "-r ${TMPDIR}/"+os.path.basename(ref_gbk)+"  ${TMPDIR}/"+breIn1+" ${TMPDIR}/"+breIn2+"\n") ## breseq version installed on my server's home dir
        else:
            sl.write(
                "~/breseq-0.35.5-Linux-x86_64/bin/breseq -j " + cpus + " --polymorphism-minimum-total-coverage-each-strand 2 "
                                                                       "--polymorphism-frequency-cutoff 0  --brief-html-output --no-javascript"
                                                                       "-o ${TMPDIR}"
                                                                       " -r ${TMPDIR}/" + os.path.basename(ref_gbk) + "  ${TMPDIR}/"+breIn1+" ${TMPDIR}/"+breIn2+"\n")  ## breseq version installed on my server's home dir

        sl.write("mv ${TMPDIR}/output ${TMPDIR}/data\n")
        sl.write("mv ${TMPDIR}/data ${TMPDIR}/${BASE/.bam/}\n")
        sl.write("ep_align_stats.py -i ${TMPDIR}/${BASE/.bam/}/reference.bam\n")  ## add investigate alignment here
        sl.write("cp -r ${TMPDIR}/${BASE/.bam/} $PWD\n")
        #sl.write("mv ${TMPDIR}/08_mutation_identification/*coverage.tab ${TMPDIR}/08_mutation_identification/"+file.replace(".bam","")+"_coverage.tsv\n")
        #sl.write('rsync -av ${TMPDIR}/08_mutation_identification/'+file.replace(".bam","")+'_coverage.tsv $PWD\n')
        #sl.write("mv ${TMPDIR}/output/output.gd ${TMPDIR}/output/"+file.replace(".bam","")+"_output.gd\n")
        #sl.write('rsync -av ${TMPDIR}/output/'+file.replace(".bam",'')+'_output.gd $PWD\n')
        #sl.write("mv ${TMPDIR}/output/output.vcf ${TMPDIR}/output/" + file.replace(".bam", "") + "_output.vcf\n")
        #sl.write('rsync -av ${TMPDIR}/output/' + file.replace(".bam", '') + '_output.vcf $PWD\n')
        sl.write('rm -rf ${TMPDIR}\n')
        sl.close()

def submit_array_investigate_alignment_files_from_one_dir(array_file, LockFile='foo.lock', ram="20000"):

    print("\n\nPreparing slurm array script...")
    cpus = "1"
    with open("array_investigate_alignment_slurm.sh",'w') as sl:
        sl.write("#!/bin/bash\n")
        sl.write("#\n")
        sl.write("#SBATCH --job-name=array_inve_align\n")
        sl.write("#SBATCH --mem="+ram+"\n")
        sl.write("#SBATCH --cpus-per-task="+cpus+"\n")
        #sl.write("#SBATCH --array=1-"+str(array_length)+'\n')
        sl.write('module load samtools\n')
        sl.write('name=$(sed -n "$SLURM_ARRAY_TASK_ID"p '+array_file+')\n')
        sl.write('OUTDIR="$(dirname "$name")"\n')
        sl.write('BASE=${name##*/}\n')
        sl.write("LOCKFILE="+LockFile+"\n")
        sl.write('set -e\n')
        sl.write('(\n')
        sl.write(' flock -w 3600 200\n')
        sl.write(" rsync -av  $name ${TMPDIR}\n")
        sl.write(") 200>$LOCKFILE\n")
        sl.write("ep_align_stats.py -i ${TMPDIR}/reference.bam\n")  ## add investigate alignment here
        sl.write("cp -r ${TMPDIR}/reference_stats.csv $OUTDIR\n")
        #sl.write("mv ${TMPDIR}/08_mutation_identification/*coverage.tab ${TMPDIR}/08_mutation_identification/"+file.replace(".bam","")+"_coverage.tsv\n")
        #sl.write('rsync -av ${TMPDIR}/08_mutation_identification/'+file.replace(".bam","")+'_coverage.tsv $PWD\n')
        #sl.write("mv ${TMPDIR}/output/output.gd ${TMPDIR}/output/"+file.replace(".bam","")+"_output.gd\n")
        #sl.write('rsync -av ${TMPDIR}/output/'+file.replace(".bam",'')+'_output.gd $PWD\n')
        #sl.write("mv ${TMPDIR}/output/output.vcf ${TMPDIR}/output/" + file.replace(".bam", "") + "_output.vcf\n")
        #sl.write('rsync -av ${TMPDIR}/output/' + file.replace(".bam", '') + '_output.vcf $PWD\n')
        sl.write('rm -rf ${TMPDIR}\n')
    sl.close()


def make_another_array_file(Directory, outname='slurm_array_investigate_alignment.txt'):
    currentDir = os.getcwd()
    os.chdir(Directory)
    os.system("ls -d $PWD/*/reference.bam| cat > ../"+outname)
    os.chdir(currentDir)



############## POLY TABLE RELATED FUNCTIONS

def concat_alignment_investigation(Directory):
    """Concatenates output files from investigate_alignment function --> reference_stats.csv"""
    current_dir = os.getcwd()
    data_dir = Directory
    Files = [pd.read_csv(data_dir+folder+'/reference_stats.csv',skiprows=0) for folder in os.listdir(data_dir)]
    column_names = "SampleID,Chrom,ReadCount,Mapped,Mapped_prc," \
              "Paired,Paired_prc,ReverseStrand,ReverseStrand_prc," \
              "MQpass,MQpass_prc,SecondaryAligned,SecondaryAligned_prc," \
              "Duplicates,Duplicates_prc,Coverage_mean,Coverage_std".split(",")
    df = pd.concat(Files, ignore_index=True)
    df.columns = column_names
    return df

def parse_vcf(vcfFile, id, meta_path=False, mode='invivo'):
    assert mode in ['invivo','invitro']

    vcf = pd.read_csv(vcfFile, sep='\t', comment="#", header=None)
    vcf = vcf[[0,1,3,4,7]]
    vcf.rename(columns={0:"Chrom", 1:"Pos", 3:"Ref",4:"Alt",7:"Info"}, inplace=True)
    vcf = vcf.assign(SampleID=[id for x in vcf.index])
    #vcf=vcf.assign(Freq=[float(re.match(re.compile("^\w+=([^;]+);.+"),x).group(1)) for x in vcf["Info"]])
    vcf=vcf.assign(AltCov=[round(float(re.match(re.compile("^\w+=[^;]+;\w+=([^;]+).*"),x).group(1))) if re.match(re.compile("^\w+=[^;]+;\w+=([^;]+).*"),x) else 0 for x in vcf["Info"]])
    vcf=vcf.assign(TotalCov=[round(float(re.match(re.compile("^\w+=[^;]+;\w+=[^;]+;\w+=([\d\.]+).*"),x).group(1))) if re.match(re.compile("^\w+=[^;]+;\w+=[^;]+;\w+=([\d\.]+).*"),x) else 0 for x in vcf["Info"]])
    vcf=vcf.assign(RefCov = round(vcf["TotalCov"] - vcf["AltCov"]))
    vcf = vcf.assign(Freq=round((vcf["AltCov"] / vcf["TotalCov"]), 3))
    vcf = vcf.assign(Type=["" for x in vcf.index])
    vcf.drop(columns=["Info"], inplace=True)
    vcf = vcf.assign(Length = vcf["Alt"].apply(lambda x:len(x)) - vcf["Ref"].apply(lambda x:len(x)) )
    vcf["Pos"] = vcf["Pos"].astype(str)
    vcf = vcf.assign(PolyID=vcf[['Chrom', 'Pos', "Ref", "Alt"]].apply(lambda x: ':'.join(x), axis=1))

    if meta_path:
        metaDf = pd.read_csv(meta_path, sep='\t')
        metaDf.SampleID = metaDf.SampleID.astype(str)
        metaDf.set_index("SampleID", inplace=True)
        if mode == "invivo":
            vcf = vcf.assign(Day = [metaDf.loc[id, "Day"] for x in vcf.index])
            vcf = vcf.assign(Mouse = [metaDf.loc[id, "Mouse"] for x in vcf.index])
            vcf = vcf.assign(Cage = [metaDf.loc[id, "Cage"] for x in vcf.index])
        else:
            vcf = vcf.assign(Day=[metaDf.loc[id, "Day"] for x in vcf.index])
            vcf = vcf.assign(Replicate = [metaDf.loc[id, "Replicate"] for x in vcf.index])
            vcf = vcf.assign(Mix = [metaDf.loc[id, "Mix"] for x in vcf.index])
        if "MetagenomeID" in metaDf.columns:
            vcf = vcf.assign(MetagenomeID=[metaDf.loc[id, 'MetagenomeID'] for x in
                                           vcf.index])  ## isolates correspond to a metagenomic sample
        else:
            vcf = vcf.assign(
                MetagenomeID=[id for x in vcf.index])  ## metagenomes have the same MetagenomeID as SampleID
        if "SampleType" in metaDf.columns:
            vcf = vcf.assign(SampleType=[metaDf.loc[id, 'SampleType'] for x in vcf.index])

    additions = pd.DataFrame()
    drop_indexes = []
    for index, row in vcf.iterrows():
        if len(row["Alt"]) == len(row["Ref"]) and len(row["Alt"]) == 1:
            vcf.loc[index,"Type"] = "SNP"
        if len(row["Alt"]) > len(row["Ref"]):
            vcf.loc[index, "Type"] = "INS"
        if len(row["Alt"]) < len(row["Ref"]):
            vcf.loc[index, "Type"] = "DEL"
        if len(row["Alt"]) == len(row["Ref"]) and len(row["Alt"]) > 1: #### replace SUB with consequtive SNPs
            vcf.loc[index,"Type"] = "SUB"
            drop_indexes += [index]
            for base_idx, base in enumerate(row.Alt):
                mod_row = row.to_frame().T.reset_index().drop(columns='index')
                mod_row.loc[0,'Alt'] = base
                mod_row.loc[0,'Ref'] =  mod_row.loc[0.,'Ref'][base_idx]
                mod_row.loc[0,'Pos'] = str(int(mod_row.loc[0,'Pos'])+ base_idx)
                mod_row.loc[0,'Type'] = "SNP"
                mod_row.loc[0,"PolyID"] = f"{mod_row.loc[0,'Chrom']}:{mod_row.loc[0,'Pos']}:{mod_row.loc[0,'Ref']}:{mod_row.loc[0,'Alt']}"
                additions = pd.concat([additions,mod_row], ignore_index=True)
    vcf.drop(index=drop_indexes, inplace=True)
    vcf = pd.concat([vcf, additions],ignore_index=True)

    ### check if there is already the same polyID in the library and if there are duplicates and a digit to differentiate them.
    duplicates = vcf[vcf.duplicated(['PolyID'], keep=False)] ### find the duplicated entries and create a new dataframe with them
    vcf = vcf.sort_values("AltCov", ascending=False).drop_duplicates("PolyID") ## keep only the first entry from each duplicated entry (max AltCov)
    assert vcf.shape == vcf.drop_duplicates("PolyID").shape
    # for name, group in duplicates.groupby("PolyID"): ## now update the unique entries in the vcf dataframe combining AltCov from duplicates
    #     if printduplicates: print(f"\nWARNING:: Duplicated polymorphisms with PolyID --> {name}\n{group}")
        #vcf.loc[vcf["PolyID"] == name, "AltCov"] = group["AltCov"].sum()
        #vcf.loc[vcf["PolyID"] == name, "Freq"] = group["AltCov"].sum() / (group["AltCov"].sum() + group["RefCov"].mean())
    return vcf

def annotate_master_table(masterDf, anno):
    #masterDf = masterDf.assign(GeneID=[np.nan for x in masterDf.index])
    #masterDf.set_index("PolyID", inplace=True)
    annoDf = pd.read_csv(anno, sep='\t')
    annoDf["Chrom"] = annoDf["Chrom"].apply(lambda x: x.strip(".1")) ### output.vcf renames the chrom to NC_004663 now instead of NC_004663.1
    annoDf = annoDf.sort_values("start")
    annoDict = {}
    print()
    for i, (index,row) in zip(track(range(masterDf.drop_duplicates("PolyID").shape[0]), description='Annotating...'), masterDf.drop_duplicates("PolyID").iterrows()):
        pos = int(row["Pos"])
        chrom = row["Chrom"]

        annotation = annoDf.loc[(annoDf["Chrom"] == chrom) & (annoDf["start"] < pos) & (annoDf["end"] >= pos),"locus_tag"] ## if there is more than one gene it concatanates with comma
        if annotation.empty:
            #print(annotation)
            annoDict[row["PolyID"]] = "Intergenic"
            #masterDf.loc[masterDf["PolyID"] == row["PolyID"],"GeneID"] = annotation.values[0] ### extracts only the first hit (the variant could span more than one gene though)
        else:
            annoDict[row["PolyID"]] = annotation.values[0]

    annoDf.drop(columns="Chrom",inplace=True)
    masterDf["locus_tag"] = masterDf['PolyID'].map(annoDict)

    annotatedDf = masterDf.merge(annoDf.drop_duplicates("locus_tag"), on="locus_tag" ,how = "left")
    annotatedDf.reset_index(inplace=True)
    #if args.annoType == "contigs":
    #    annotatedDf.drop(columns="GeneID", inplace=True)
    #    annotatedDf.rename(columns={"RefGeneID":"GeneID","RefGeneName":"GeneName", "RefProteinID":"ProteinID", "RefProteinInfo":"ProteinInfo"}, inplace=True)
    return annotatedDf

def assign_synonymy(annotatedDf: object) -> object:
    """Annotated masterDf dataframe as input."""
    annotatedDf = annotatedDf.assign(Synonymy=["undefined" for x in annotatedDf.index])
    ref_path = '/home/christos/bacevo/reference/Btheta-VPI-5482/GCF_000011065.1/GCF_000011065.1_ASM1106v1_genomic.fna'
    ref_fasta = fasta_to_df(ref_path)
    funny_cds = []

    print("\n")
    for i, (index, row) in zip(track(range(annotatedDf.loc[annotatedDf["locus_tag"] != "Intergenic"].shape[0]),description='Assigning synomy...'),annotatedDf.loc[annotatedDf["locus_tag"] != "Intergenic"].iterrows()): ### select only intragenic polys
        if row["Type"] in ["SNP","INS","DEL", "SUB"] and isinstance(row['translation'], str) and len(row.cds) % 3 == 0: ## this will select only protein coding genes.
            geneStart = row["start"]
            alt = row["Alt"]
            pos = int(row["Pos"])
            ref = row["Ref"]
            orientation  = float(row["strand"])
            """If i have mutation on the first position of the cds --> index == 0 i need to change the pos_on_gene because gets -1 """
            idx_on_gene = int(pos - float(geneStart))-1 ## this is a subtraction of python indexes (relative to the genome sequence)
            #print(pos_on_gene, orientation, len(row.cds))
            #print(row['locus_tag'])


            #print(row['cds'])
            Type = row["Type"]
            if orientation == -1:
                DNA = str(Seq(row["cds"]).reverse_complement())
            else:
                DNA = row["cds"]
            protein = Seq(DNA).translate(table=11)
            mutatedProtein = ""
            mutatedDNA = list(str(DNA))
            if Type == "INS":
                proins = mutatedDNA[:idx_on_gene]
                preins = mutatedDNA[idx_on_gene:]
                #mutatedDNA = Seq("".join(proins + list(alt) + preins), generic_dna)
                plus = abs(len(alt) - len(ref))
                #print(f'Length of insertion {plus}')
                mutatedDNA = Seq("".join(proins + list(alt[plus:]) + preins))
                mutatedProtein = mutatedDNA.translate(table=11)
            if Type == "SNP":
                #print(f'Mutation occured in index {idx_on_gene} out of {len(row["cds"])} sites')
                if ref != DNA[idx_on_gene]:
                    if idx_on_gene >= 0:
                        logging.warning("Ref does not match cds! Index >= 0")
                        print(f"\tWARNING:: PolyID:{row.PolyID} | Locus:{row['locus_tag']} | strand: {orientation}|Start:{row.start} POS:{pos} - index on gene {idx_on_gene}, Type: {row['Type']}, Ancestry:{row['Ancestry']}")
                        print(f'{ref_fasta.loc[ref_fasta.Chrom == row.Chrom+".1","seq"].item()[int(row.start):pos-1]}{color.RED}{ref_fasta.loc[ref_fasta.Chrom == row.Chrom+".1","seq"].item()[pos-1]}{color.END}{ref_fasta.loc[ref_fasta.Chrom == row.Chrom+".1","seq"].item()[pos:int(row.end)]}')
                        print(f'length in fasta: {len(ref_fasta.loc[ref_fasta.Chrom == row.Chrom+".1","seq"].item()[int(row.start):int(row.end)])}')
                        print(f'\n{DNA[:idx_on_gene]}{color.RED}{DNA[idx_on_gene]}{color.END}{DNA[idx_on_gene+1:]} ')
                        print(f'length of cds: {len(DNA)}')
                        print(f'\t\tRef --> {ref} != {DNA[idx_on_gene]} | GeneStart {ref_fasta.loc[ref_fasta.Chrom == row.Chrom+".1","seq"].item()[int(float(geneStart))]}={DNA[0]}, Alt = {alt}\n\n')
                    else:
                        print(f"\tWARNING:: PolyID:{row.PolyID} | Locus:{row['locus_tag']} | strand: {orientation}|Start:{row.start} POS:{pos} - index on gene {idx_on_gene}, Type: {row['Type']}, Ancestry:{row['Ancestry']}")
                        print(f'{ref_fasta.loc[ref_fasta.Chrom == row.Chrom+".1","seq"].item()[int(row.start):pos-1]}{color.RED}{ref_fasta.loc[ref_fasta.Chrom == row.Chrom+".1","seq"].item()[pos-1]}{color.END}{ref_fasta.loc[ref_fasta.Chrom == row.Chrom+".1","seq"].item()[pos:int(row.end)]}')
                        print(f'length in fasta: {len(ref_fasta.loc[ref_fasta.Chrom == row.Chrom+".1","seq"].item()[int(row.start):int(row.end)])}')
                        print(f'\n{DNA[:idx_on_gene]}{color.RED}{DNA[idx_on_gene]}{color.END}{DNA[idx_on_gene+1:]} ')
                        print(f'length of cds: {len(DNA)}')
                        print(f'\t\tRef --> {ref} != {DNA[idx_on_gene]} | GeneStart {ref_fasta.loc[ref_fasta.Chrom == row.Chrom+".1","seq"].item()[int(float(geneStart))]}={DNA[0]}, Alt = {alt}\n\n')

                    funny_cds += [row['locus_tag']]

                mutatedDNA[idx_on_gene] = alt
                mutatedDNA = Seq("".join(mutatedDNA))
                mutatedProtein = mutatedDNA.translate(table=11)
            if Type == "DEL":
                minus = len(alt) - len(ref)
                #print(pos_on_gene)
                #print(ref)
                #print(alt)
                #print(minus)
                #print(DNA[:pos_on_gene-1]+color.RED+DNA[pos_on_gene:]+color.END)
                #for index, base in enumerate(list(alt)):
                    #mutatedDNA[pos_on_gene + index] = ""

                for index, base in enumerate(list(ref[minus:])):
                    if idx_on_gene + len(ref)+minus+index < len(DNA):
                        mutatedDNA[idx_on_gene+ len(ref)+minus+ index] = ""

                mutatedDNA = Seq("".join([base for base in mutatedDNA if base != ""]))
                mutatedProtein = mutatedDNA.translate(table=11)

            if Type == "SUB":
                minus = len(alt) - len(ref)
                for index, base in enumerate(list(alt)):
                    if idx_on_gene + len(ref)+minus+index < len(DNA):
                        mutatedDNA[idx_on_gene+index] = base
                mutatedDNA = Seq("".join(mutatedDNA))
                mutatedProtein = mutatedDNA.translate(table=11)

            if str(protein) == str(mutatedProtein):
                annotatedDf.loc[annotatedDf["PolyID"] == row["PolyID"], "Synonymy"] = "Synonymous"
            else:
                annotatedDf.loc[annotatedDf["PolyID"] == row["PolyID"], "Synonymy"] = "Nonsynonymous"
    print(pd.Series(funny_cds).value_counts())
    return annotatedDf

def add_snp_type(annotatedDf):
    print(color.BOLD+"\nAdding SNP types (transition - transversion)...\n"+color.END)
    transDict = {"A":"Purine", "G":"Purine", "T":"Pyrimidine","C":"Pyrimidine"}
    annotatedDf = annotatedDf.assign(MutationType=["undefined" for x in annotatedDf.index])
    for name, group in annotatedDf.loc[annotatedDf["Type"] == "SNP"].groupby("PolyID"):
        # alt = group["Alt"].values[0]
        # ref = group["Ref"].values[0]
        if transDict[group["Alt"].values[0]] == transDict[group["Ref"].values[0]]:
            annotatedDf.loc[annotatedDf["PolyID"] == name, "MutationType"] = "Transition"
        else:
            annotatedDf.loc[annotatedDf["PolyID"] == name, "MutationType"] = "Transversion"
    return annotatedDf

def update_with_ancestral_polys(df, ancestralFreqThres=0.9, bioUnit="Replicate", condition='Mix', name = 'bern',onlyAncestry=True):
    """Define ancestral polymorphisms..and update the dataframe with this information
    Set a threshold for fixed polymorphisms and choose for bioUnit between Replicate and Mouse and for condition Mix or Cage respectively."""
    print("\nAdding ancestry ('ancestral' or 'evolved').")
    # sdf["Ancestry"] = ["ancestral" if row["PolyID"] in sdf.loc[sdf[bioUnit] == 0, "PolyID"].tolist() else "evolved" for index,row in sdf.iterrows()]
    # sdf = sdf.loc[~sdf["PolyID"].isin(sdf.loc[sdf["Mouse"] ==0, "PolyID"].tolist())]
    ancestralDf = df.loc[df[bioUnit] == 0].sort_values("Freq")
    #ancestralDf.set_index("PolyID", inplace=True)
    #ancestralDf["Pos"] = ancestralDf["Pos"].astype(str)
    print(f"\t::Writing ancestral polymorphism to {name}_AncestralMutations.csv'.")
    ancestralDf.to_csv(f"{name}_AncestralMutations.csv")
    fixed  = ancestralDf.loc[ancestralDf["Freq"].ge(ancestralFreqThres)]
    df["Ancestry"] = ['' for x in df.index]
    for i, row in df.iterrows():
        if row.PolyID in fixed.PolyID.values:
            df.loc[i, "Ancestry"] = 'fixed'
        elif row.PolyID in ancestralDf.PolyID.values:
            df.loc[i, "Ancestry"] = 'ancestral'
        else:
            df.loc[i, 'Ancestry'] = 'evolved'
    #df["Ancestry"] = ["ancestral" if row["PolyID"] in ancestralDf.PolyID.values else "evolved" for index, row in df.iterrows()]
    #df['Ancestry'] = df['Ancestry'].apply(lambda x:"fixed" if )
    """Remove biounit 0 (day 0 )"""
    print("\t::Removing the ancestral polymorphisms from the dataset.")
    df = df.loc[df[bioUnit] != 0]

    ### remove fixed polymorpisms in the ancestral
    assert ancestralDf.shape[0] == len(set(ancestralDf.index.values))

    if onlyAncestry is False:
        print("Adding Reversed column for different evolutionary scenarios.")
        """If a polymorphism is present in the ancestral population (fixed) and then it is not present in any of the evolved libraries, it would not be reported as polymorphism. Thus 
        I have to recreate these reversed to the reference state events. If it is present then the frequency should be 1 - current freq"""
        print(
            "Update the dataset using fixed polymorphisms in the ancestral population. Fixed polymorphisms should have frequency > {}".format(
                ancestralFreqThres))
        print(
            "Will create a Reversed column with entries: 0: for de novo mutations, 1:for Reversion (adds entry), 2: mutation to sth else than ref with freq higher than the reference,"
            "3: mutaion to sth else than ref with freq lower than ref state, 4: Ancestral polymorphims that remains fixed., 5: Reversion to ancestral state but not fixed. ")
        df = df.assign(Reversed=[999 for x in df.index])
        df.loc[~df.PolyID.isin(ancestralDf.index.values), "Reversed"] = 0

        revdf = pd.DataFrame(columns=df.columns)
        from datetime import datetime
        startTime = datetime.now()
        entries = ancestralDf.shape[0]
        for name, group in df.groupby("SampleID"):
            assert len(group["PolyID"].values) == len(group["PolyID"].drop_duplicates().values)
            # group.set_index("PolyID", inplace=True)
            print("\nSample: {}".format(name))

            entryIndex = 0
            for polyID in ancestralDf.index.values:  ### for only fixed mutations in the ancestral population (above the threshold)
                entryIndex += 1
                sys.stdout.write('\r' + color.RED + color.BOLD + "Running time: " + str((datetime.now() - startTime)) + " (" + str(
                        round((entryIndex * 100) / entries, 1)) + "%)" + color.END)
                ancState = ancestralDf.loc[polyID, "Alt"]
                refState = ancestralDf.loc[polyID, "Ref"]
                if polyID in group["PolyID"].values:  ### scenarios 4,5 in the table (google slides: Genetic diversification in the murine gut.. - slide: Ancestral Population)
                    # print(group.loc[group["PolyID"] == polyID, ["PolyID","Freq"]])
                    idx = (df['PolyID'] == polyID) & (df["SampleID"] == name)
                    df.loc[idx, ["Alt", "Ref"]] = df.loc[idx, ["Ref", "Alt"]].values
                    df.loc[idx, ["RefCov", "AltCov"]] = df.loc[idx, ["AltCov", "RefCov"]].values
                    # print(df.loc[idx, "Freq"].apply(lambda x: 1.0 - float(x)).values)
                    df.loc[idx, "Freq"] = df.loc[idx, "Freq"].apply(lambda x: 1.0 - float(x)).values
                    # print(df.loc[idx, "Freq"])
                    df.loc[idx, "PolyID"] = ancestralDf.loc[polyID, ["Chrom", "Pos", "Alt", "Ref"]].str.cat(sep=":")
                    if group.loc[group["PolyID"] == polyID, "Freq"].values >= ancestralFreqThres:  ## scenario 4 in the table
                        df.loc[idx, "Reversed"] = 4
                        df.loc[idx, "Ancestry"] = 'ancestral'
                        # print("4")
                    else:
                        df.loc[idx, "Reversed"] = 5
                        df.loc[idx, "Ancestry"] = 'evolved'
                        # print("5")
                else:  ## scenarios 1,2,3
                    if ":".join(polyID.split(":")[:-1]) not in [":".join(x.split(":")[:-1]) for x in group["PolyID"].values]:  ## scenario 1 that there is no detected polymorphism
                        newRow = ancestralDf.loc[polyID,].copy()
                        newRow["Ref"] = ancState
                        newRow["Alt"] = refState
                        newRow["RefCov"] = 0
                        newRow["AltCov"] = 100
                        newRow["TotalCov"] = 100
                        newRow["PolyID"] = newRow["Chrom"] + ":" + newRow["Pos"] + ":" + ancState + ":" + refState
                        newRow["SampleID"] = group["SampleID"].values[0]
                        newRow[bioUnit] = group[bioUnit].values[0]
                        newRow[condition] = group[condition].values[0]
                        newRow["Day"] = group["Day"].values[0]
                        newRow["Freq"] = 1.0
                        newRow["Reversed"] = 1
                        newRow["Ancestry"] = "evolved"
                        revdf = revdf.append(newRow, ignore_index=True)
                        # print("1")
                    else:  ##senaria 2,3
                        for index, row in group.iterrows():
                            if ancestralDf.loc[polyID, ["Chrom", "Pos", "Ref"]].equals(row[["Chrom", "Pos", "Ref"]]) and ancestralDf.loc[polyID, "Alt"] != row["Alt"]:
                                # if ancestralDf[polyID, "Chrom"] == row["Chrom"] and str(ancestralDf[polyID, "Pos"]) == str(row["Pos"]) and
                                # if ":".join(polyID.split(":")[-1]) == ":".join(entry.split(":")[-1]) and group[entry,"Alt"] == refState
                                idx = (df['PolyID'] == polyID) & (df["SampleID"] == name)
                                df.loc[
                                    idx, "Ref"] = ancState  ## in both cases the new Reference is the ancestral state (fixed alternative in the ancestral population)
                                df.loc[idx, "Ancestry"] = "evolved"
                                df.loc[idx, "PolyID"] = ancestralDf.loc[polyID, ["Chrom", "Pos", "Alt", "Ref"]].str.cat(
                                    sep=":")
                                if float(row["Freq"]) > 0.5:  # scenario 2
                                    df.loc[idx, "Reversed"] = 2
                                    # print("2")
                                else:  ## scenario 3
                                    df.loc[idx, "Alt"] = row["Ref"]

                                    df.loc[idx, "Freq"] = df.loc[idx, "Freq"].apply(lambda x: 1.0 - float(x)).values
                                    df.loc[idx, "Reversed"] = 3
                                    # print("3")
                                df.loc[idx, "RefCov"] = (df.loc[idx, "AltCov"] / df.loc[idx, "Freq"]).astype(int).values  ### estimate RefCov from frequency and altCov.
                                df.loc[idx, "TotalCov"] = df.loc[idx, "RefCov"] + df.loc[idx, "AltCov"].values

        originalpols = df["PolyID"].unique()
        reversepols = revdf["PolyID"].unique()
        revdf["Day"] = revdf["Day"].astype(int)
        revdf[bioUnit] = revdf[bioUnit].astype(int)
        revdf["Pos"] = revdf["Pos"].astype(int)
        revdf["Length"] = revdf["Length"].astype(int)
        revdf["SampleID"] = revdf["SampleID"].astype(str)

        df = df.append(revdf, ignore_index=True)

    return df

#################################################
def tsne_from_dataframe(CountsDf , metric='euclidean', find_perplexity=True):
    import scipy
    from sklearn.manifold import TSNE
    from scipy.spatial.distance import squareform
    from scipy.spatial.distance import pdist
    from sklearn.preprocessing import StandardScaler
    import plotly.express as px
    X = squareform(pdist(CountsDf.values, metric=metric))
    scaler = StandardScaler()
    #X_norm = scaler.fit_transform(X)
    if find_perplexity:
        perplexity = np.arange(5,100, 5)
        divergence = []

        for i in perplexity:
            model = TSNE(n_components=2, init="pca", perplexity=i)
            reduced = model.fit_transform(X)
            divergence.append(model.kl_divergence_)

        fig = px.line(x=perplexity, y=divergence, markers=True)
        fig.update_layout(xaxis_title="Perplexity Values", yaxis_title="Divergence")
        fig.update_traces(line_color="red", line_width=1)
        fig.show()

    tsne = TSNE(n_components=2, random_state=42, perplexity=80)
    X_tsne = tsne.fit_transform(squareform(pdist(CountsDf.values, metric=metric)))
    print(f' kl divergence: {tsne.kl_divergence_}')
    print(X_tsne)
    fig = px.scatter(x=X_tsne[:, 0], y=X_tsne[:, 1], color=CountsDf.index.get_level_values(1))
    fig.update_layout(
    title="t-SNE visualization of Customer Churn dataset",
    xaxis_title="First t-SNE",
    yaxis_title="Second t-SNE",)
    fig.show()

    return pd.DataFrame(data=X_tsne, columns=['tsne1','tsne2'],index=CountsDf.index.values)

    from skbio.stats.distance import DistanceMatrix
    #PcoaDf = My_pcoa.samples[["PC1","PC2"]].set_index(CountsDf.index)
    #pc1_prop_explained, pc2_prop_explained = My_pcoa.proportion_explained.loc["PC1"], My_pcoa.proportion_explained.loc["PC2"]
    #return PcoaDf, pc1_prop_explained, pc2_prop_explained


def pcoa_from_dataframe(CountsDf , metric='braycurtis'):
    import scipy
    from skbio.stats.distance import DistanceMatrix
    from skbio.stats.ordination import pcoa
    My_pcoa = pcoa(scipy.spatial.distance.pdist(CountsDf.values, metric=metric))
    PcoaDf = My_pcoa.samples[["PC1","PC2"]].set_index(CountsDf.index)
    pc1_prop_explained, pc2_prop_explained = My_pcoa.proportion_explained.loc["PC1"], My_pcoa.proportion_explained.loc["PC2"]
    return PcoaDf, pc1_prop_explained, pc2_prop_explained

def benchmark_coverage_threshold(Df):
    Console.log(f'Evaluate relationship between genome coverage and number of polys.')
    """Filters polys for a range of TotalCov values and counts the number of polys - Compares this count to the mean library coverage"""
    nf = pd.DataFrame()
    cov_range = [x for x in range(0,1000,50)]
    for cov_thres in track(cov_range, description='Evaluating Polys for coverage threshold...'):
        tf = Df.loc[Df["TotalCov"].astype(float) >= cov_thres]
        tf["PolyCount"] = tf.groupby("SampleID")["PolyID"].transform(len)
        tf['TotalCovThres'] = [cov_thres for x in tf.index]
        tf = tf.drop_duplicates(['SampleID',"Chrom"])[["SampleID","PolyCount", "Chrom", "Coverage_mean","Coverage_std", "TotalCovThres"]]
        tf['cor'] = tf["Coverage_mean"].corr(tf.PolyCount)
        nf = nf.append(tf)
    return nf



    ### sliding or not window function
    def window(lista, w, o):
        """Iterates over a list l and calculates a value for a window of width w and overlap o"""
        if o == 0:
            win_list = [sum(lista[i:i + w]) / float(w) for i in range(0, len(lista), w)]
        else:
            win_list = [sum(lista[i:i + w]) / float(w) for i in range(0, len(lista), w - o)]
        return win_list

    #### Estimate coverage all over the genome
    # samtools =subprocess.Popen("samtools depth -r "+args.genome+" -a "+args.inDir,shell=True, stdout=subprocess.PIPE)
    # shell_out=samtools.communicate()[0].split()
    ### Estimating reads stats and Nucleotide Identity Percentage distribution


    cpDf = pd.DataFrame({"Chrom":[],"Pos":[],"Coverage":[]})
    for chrom in ["NC_004663","NC_004703"]:
        pos, cov = window([x[1] for x in pysamstats.load_coverage(bamfile, chrom=chrom)], 100,0), window([x[2] for x in pysamstats.load_coverage(bamfile, chrom=chrom)],100,0) ### this will subsample in a window fashion
        cpDf = cpDf.append(pd.DataFrame({"SampleID":[SampleID for x in pos],"Chrom":[chrom for x in pos],'Pos':pos,'Coverage':cov}))
    cpDf.to_csv('Coverage_by_Position.csv',index=None)
    cDf = cpDf.groupby("Chrom").agg({"Coverage":["mean",'std']})
    cDf.columns = cDf.columns.map('_'.join).str.strip('_')
    cDf.reset_index(inplace=True)
    odf = pd.DataFrame()
    for chrom in chroms:
        mq_reads = 0
        secondary = 0
        duplicate = 0
        proper_pair = 0
        paired = 0
        reverse = 0
        qc_fail = 0
        # nucl_ide=[]
        # read_pos=[]
        NIP_on_genome = {}
        pairs_distance = {}
        mapping_qualities = []
        for read in bamfile.fetch(chrom):
            mapping_qualities.append(float(read.mapping_quality))
            if read.mapping_quality > 20:
                mq_reads += 1
            if read.is_secondary:
                secondary += 1
            if read.is_duplicate:
                duplicate += 1
            if read.is_proper_pair:
                proper_pair += 1
            if read.is_reverse:
                reverse += 1
            if read.is_paired:
                paired += 1
            if read.is_qcfail:
                qc_fail += 1
            # if read.reference_start not in pairs_distance:
            #    pairs_distance[read.reference_start] = abs(read.reference_start - read.next_reference_start)
            read_len = 0
            nui = 0
            NPI = 0

            # if read.cigartuples is not None:
            #     for t in read.cigartuples:
            #         if t[0]== 0:
            #             nui+=float(t[1])
            #         read_len+=float(t[1])
            #
            #     NPI=nui/read_len
            #     #nucl_ide.append(NPI)
            #     mid_pos=0
            #     mapping=read.get_aligned_pairs()
            #     for i in range(len(mapping)):
            #         if mapping[i][1] is not None and mapping[len(mapping)-1-i][1] is not None:
            #             mid_pos+=(int(mapping[i][1])+int(mapping[len(mapping)-1-i][1]))/2
            #             break
            #         else:
            #             #print("Deletion or insertion")
            #             continue

            # NIP_on_genome[mid_pos] = NPI

        # df = pd.DataFrame(np.array(mapping_qualities), columns=["MQ"])
        # ax = df.plot(kind='density')
        # lines, labels = ax.get_legend_handles_labels()
        # ax.legend(lines,labels)
        # plt.show()

        print("Number of reads: {}".format(int(all_reads)))
        print("% of mapped reads: {}".format(round(100 * mapped / all_reads, 3)))
        print("% of unmapped reads: {}".format(round(100 * unmapped / all_reads, 3)))
        print("% with a pair: {}".format(round(100 * paired / all_reads, 3)))
        print("% of proper pairs: {}".format(round(100 * proper_pair / all_reads, 3)))
        print("% derive from reverse strand: {}".format(round(100 * reverse / all_reads, 3)))
        print("% QC filter fail: {}".format(round(100 * qc_fail / all_reads, 3)))
        print("% of reads with MQ > 20: {}".format(round(100 * mq_reads / all_reads, 3)))
        print("% secondary alignments: {}".format(round(100 * secondary / all_reads, 3)))
        print("% of duplicates: {}".format(round(100 * duplicate / all_reads, 3)))

        tdf = pd.DataFrame({"SampleID":[SampleID],"Chrom":chrom,"ReadCount": [all_reads], "Mapped": [mapped], "Mapped_prc": [round(unmapped / all_reads, 3)],
                            "Paired": [proper_pair], "Paired_prc": [round(paired / all_reads, 3)],
                            "ReverseStrand": [reverse], "ReverseStrand_prc": [round(reverse / all_reads, 3)],
                            "MQpass": [mq_reads], "MQpass_prc": [round(mq_reads / all_reads, 3)],
                            "SecondaryAligned": [secondary], "SecondaryAligned_prc": [round(secondary / all_reads, 3)],
                            "Duplicates": [duplicate], "Duplicates_prc": [round(duplicate / all_reads, 3)]})
        odf = odf.append(tdf)
    odf = odf.merge(cDf, on='Chrom')
    return odf

def distance_between_runs(df, binary=True):
    """ PcoA of samples based on variant pattern both as presence/absence and with frequencies"""
    outDf = pd.DataFrame()
    for x in np.arange(0, 20, 5):
        for y in np.arange(0.2,1, 0.1):
            for z in [["shared"],["all",'shared','VBCF','BSF']]: ## added all label just to be able to use it as label for the filter later
                filtering = (df['RunPrevalence'] < y) & (df.AltCov > x) & (df.SequencingPattern.isin(z))
                totalPolys = df.loc[filtering].drop_duplicates("PolyID").shape[0]
                table = pd.pivot(data=df.loc[filtering], index='SampleID', columns='PolyID', values="Freq").fillna(0)
                binary_table = table.copy()
                if binary is True:
                    binary_table[binary_table > 0] = 1
                import seaborn as sns
                from skbio.stats.distance import DistanceMatrix
                Distance_matrix = DistanceMatrix(scipy.spatial.distance.pdist(binary_table.values)).data
                dm = pd.DataFrame(index=binary_table.index.values, columns=binary_table.index.values, data=Distance_matrix)
                dm = dm.stack().reset_index().rename(columns={'level_0': 's1', 'level_1': 's2', 0: 'Distance'})
                RunDistance = []
                for i, row in dm.iterrows():
                    if (str(row['s1']).startswith("5") and str(row['s2']).startswith("6")) or (str(row['s1']).startswith("6") and str(row['s2']).startswith("5")):  ### select only library pairs from different runs
                        RunDistance.append(row['Distance'])
                outDf = outDf.append(pd.DataFrame({"RunDistanceMean":[np.mean(RunDistance)], 'RunDistanceStd':[np.std(RunDistance)], "AltCovFilter":[x], "RunPrevalenceFilter":[y], "SequencingPatternFilter":z[0], "PolyCount":[totalPolys]}))
    return outDf

def compare_isolates_to_metagenomes(df):

    metagenomes_with_isolates = [ row.MetagenomeID for i,row in
                                  df.drop_duplicates("SampleID").groupby('MetagenomeID').size().reset_index().iterrows()
                                  if row[0] > 1]
    for name, group in df.loc[df.MetagenomeID.isin(metagenomes_with_isolates)].groupby("MetagenomeID"):
        N_isolates = group.drop_duplicates("SampleID").shape[0] - 1 ## subtract the metagenomic sample
        N_polys = group.drop_duplicates("PolyID").shape[0]
        print(f"{N_isolates} isolates and {N_polys} different polys for metagenome {name}.")
        #group.SampleID = group.SampleID.astype(str)
        ## make matrix with poly catalog based on the union.
        cat = pd.pivot(data=group, index="SampleID", columns="PolyID", values="Freq").fillna(0)
        ##find those present in metagenomes BUT absent in isolates
        isolates = [x for x in cat.index.values if x != name]
        print(name)
        print(len(isolates))
        metagenome_polys = cat.loc[:,cat.drop(index=isolates).sum(axis=0) > 0].shape[1]
        not_in_isolates = cat.loc[:,cat.drop(index=name).sum(axis=0) == 0].shape[1] ## only metagenome has the same sampleID with metagenomeID ;) quite handy
        isolate_polys = cat.loc[:, cat.drop(index=name).sum(axis=0) > 0].shape[1]
        not_in_metagenome = cat.loc[:,cat.drop(index=isolates).sum(axis=0) == 0].shape[1]

        fraction_in_metagenome = 1 - (not_in_metagenome / isolate_polys)
        fraction_in_isolates = 1 - (not_in_isolates / metagenome_polys)
        ##find those present in isolates BUT absent from metagenome
        #len([poly for poly in cat.columns if cat.loc[name,poly] == 0 and cat[poly].sum() > 0])
        print(f'{fraction_in_isolates:.2%} of metagenome polys present in the isolates and \n'
              f'{fraction_in_metagenome:.2%} of polys from isolates present in the metagenome\n\n' )
def add_group_poly_prevalence_and_mean(df, column='', withSampleType=True):
    if withSampleType:

        group_sizes = df.drop_duplicates("SampleID")[[column,'SampleType']].value_counts().to_dict() ## size of groups
        group_sizes = {"".join([str(x) for x in k]):v for k,v in group_sizes.items()}
        #print(group_sizes)
        df['tmp'] = df[column].astype(str) + df['SampleType']
        df['tmpsize'] = df["tmp"].map(group_sizes)
        df[column+"Prevalence"] = df.groupby(['PolyID',column,"SampleType"]).Freq.transform(len)
        df[column+"PrevalencePercent"] = df[column+'Prevalence'] / df['tmpsize']
        df[column+"MeanFreq"] = df.groupby(["PolyID",column, 'SampleType']).Freq.transform("mean")
        #print(df[[column+'Prevalence',column+'PrevalencePercent', "SampleType",'PolyID']].head(10))
        df.drop(columns=['tmp','tmpsize'], inplace=True)



def plot_polycount_by_day_type(df, save=False,plot=True):
    """Polymorphisms count per Day"""
    CountByDayDf = df.groupby(["Cage", "Mouse", "Day", "Type"]).size().unstack(fill_value=0).reset_index()
    CountByDayDf = CountByDayDf.assign(Total=CountByDayDf[df["Type"].unique().tolist()].sum(axis=1))
    m = pd.melt(CountByDayDf, value_vars=df["Type"].unique().tolist().append("Total"),
                value_name="Count", id_vars=["Cage", "Mouse", "Day"], var_name="Type")
    m.sort_values(by=["Mouse", "Day"], inplace=True)
    m["Day"] = m["Day"].astype(int)
    m["Mouse"] = m["Mouse"].astype(int)
    # m.to_csv("PolyCountsR.csv", sep='\t', index=False)
    g = sns.FacetGrid(data=m, col="Type", hue="Mouse", sharey=False, legend_out=True)
    g.map(sns.lineplot, "Day", "Count", estimator=None)
    g.set(xticks=sorted(m["Day"].unique()))
    g.add_legend()
    if save: g.savefig("Figures/1.PolymorphismsCount.svg")
    g.add_legend()
    if plot: plt.show()
    plt.close("all")
def plot_polycount_by_synonymy(df, save=False,plot=True):
    CountSynNonSynDf = df.groupby(["Mouse", "Day", "Synonymy"]).size().unstack(fill_value=0).reset_index()
    m = pd.melt(CountSynNonSynDf, value_vars=df["Synonymy"].unique().tolist(), value_name="Count",
                id_vars=["Mouse", "Day"], var_name="Synonymy")
    m.sort_values(by=["Mouse", "Day"], inplace=True)
    m["Day"] = m["Day"].astype(int)
    m["Mouse"] = m["Mouse"].astype(int)
    h = sns.FacetGrid(m, col="Mouse", hue="Synonymy")
    h.map(sns.lineplot, "Day", "Count", estimator=None)
    h.add_legend()
    h.set(xticks=sorted(df["Day"].unique()))
    if save: h.savefig("Figures/2.SynNonSynCount.svg")
    if plot: plt.show()
    plt.close("all")
def plot_tstv_by_day_mouse(df, save=False, plot=True):
    CountTransDf = df.loc[df["Type"] == "SNP"].groupby(["Mouse", "Day", "MutationType"]).size().unstack(
        fill_value=0).reset_index()
    CountTransDf = CountTransDf.assign(Ratio=CountTransDf["Transition"] / CountTransDf["Transversion"])
    h = sns.FacetGrid(CountTransDf, col="Mouse", col_wrap=2, sharex=True)
    h.map(sns.lineplot, "Day", "Ratio")
    h.set(xticks=sorted(CountTransDf["Day"].unique()))
    h.set_ylabels("Ts / Tv")
    for ax in h.axes.flat:
        title = ax.get_title()
        thisMouse = int(re.match(re.compile("[^\d]*(\d+).*"), title.split("|")[0]).group(1))
        # sns.lineplot(x="Day", y="Freq", data=ancestralDf.loc[ancestralDf[bioUnit] == thisMouse], color="magenta", legend=False )
        # thisDay = int(re.match(re.compile("[^\d*]*(\d+).*"),title.split("|")[1]).group(1))
        axdf = df.loc[df["Mouse"] == thisMouse, ["Day", "PolyCount"]].drop_duplicates("Day").sort_values("Day")
        ax2 = ax.twinx()
        # ax2.plot(axdf["Day"], axdf["PolyCount"], color="red",linestyle='--', linewidth=2, markersize=12)
        sns.lineplot(x="Day", y="PolyCount", data=axdf, color="r", linestyle="--", ax=ax2)
        ax2.set_ylabel("Polymorphisms")
        ax2.grid(False)
        # ax2.set_ylabel("PolyCount", fontsize=10, color='red')
    h.add_legend()
    h.set_axis_labels("Days post inoculation")
    if save: h.savefig("Figures/2b.TsTvRatio.svg")
    if plot: plt.show()
    plt.close("all")
    return CountTransDf
def plot_counts_by_day(df, save=False,plot=True):
    CountByDayDf = df.groupby(["Cage", "Mouse", "Day", "Type"]).size().unstack(fill_value=0).reset_index()
    CountByDayDf = CountByDayDf.assign(Total=CountByDayDf[df["Type"].unique().tolist()].sum(axis=1))
    f = plt.figure()
    ax = sns.boxplot(x="Day", y="Total", data=CountByDayDf, palette="Set3")
    ax = sns.swarmplot(x="Day", y="Total", data=CountByDayDf, dodge=True, color=".2")
    ax.set_xlabel("Days post inculation", fontsize=16)
    ax.set_ylabel("Polymorphisms", fontsize=16)
    if save: f.savefig("Figures/4.PolyByDay.svg")
    if plot: plt.show()
    from scipy import stats
    np.random.seed(1)
    import itertools
    for x in itertools.combinations(CountByDayDf["Day"].drop_duplicates(), 2):
        t, p = stats.ttest_ind(CountByDayDf.loc[CountByDayDf["Day"] == x[0], "Total"],
                               CountByDayDf.loc[CountByDayDf["Day"] == x[1], "Total"], equal_var=False)
        print(str(x[0]) + ":" + str(x[1]) + " Ttest statistic: {}, p-value: {}".format(t, p))
    plt.close("all")
def plot_polycount_coverage(df,save=False, plot=True):
    from scipy.stats import linregress
    slope, intercept, Rvalue, pvalue, stderr, = linregress(df.drop_duplicates("SampleID")["Coverage_mean"],
                                                           df.drop_duplicates("SampleID")["PolyCount"])
    print("{}, {}, {}, {}".format(slope, Rvalue, pvalue, stderr))
    g = sns.jointplot(x="Coverage_mean", y="PolyCount", data=df.drop_duplicates("SampleID"), kind="reg", color="m",
                      height=7)
    g.set_axis_labels(" Mean Genome Coverage", "Polymorphisms")
    if plot: plt.show()
    if save:g.savefig("Figures/5b.PolyCountCoverage.svg")
    plt.close("all")
def plot_polys_on_genome(df, save=False, plot=True):
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    # Initialize the FacetGrid object
    df["LMouse"] = df['Mouse'].apply(lambda x: "Mouse " + str(x))
    df["LDay"] = df['Day'].apply(lambda x: "Day " + str(x))
    df["Label"] = df[["LMouse", "LDay"]].astype(str).apply(lambda x: " - ".join(x), axis=1)
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    g = sns.FacetGrid(df.sort_values(["Mouse", "Day"]), row="Label", hue="Label", aspect=25, height=.5, palette=pal)
    # Draw the densities in a few steps
    g.map(sns.kdeplot, x="Pos", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
    g.map(sns.kdeplot, x="Pos", clip_on=False, color="w", lw=1.8, bw=.2)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)

    g.map(label, "Pos")
    # Set the subplots to overlap
    # g.fig.subplots_adjust(hspace=-.1)
    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[])
    plt.show()

def poly_turnover(df, save=False,plot=True, NameTag="", outdir=''):
    if NameTag != "": NameTag += "_"
    turn = pd.DataFrame(
        columns=["Dataset", "Cage","Mouse", "Day", "PolyID", "Synonymy", "new", "seen", "justseen", "lost", "justlost", "cumseen"])
    df["Cage"] = df.Cage.astype(str)
    # df["seen"] = [0 for x in df.index]
    # df["justseen"] = [0 for x in df.index]
    # df["new"] = [0 for x in df.index]
    # df["justlost"] = [0 for x in df.index]
    # df["lost"] = [0 for x in df.index]
    days = df["Day"].drop_duplicates().sort_values().tolist()
    from pathlib import Path
    # my_file = Path(currentDir+"turn.obj")
    # if os.path.isfile(my_file):
    #     turn = pickle.load(open("tmp/turn.obj",'rb'))
    # else:
    # turnDf = pd.melt(turnDf,id_vars=[bioUnit,"Day"], value_vars=["new","seen","justseen","lost","justlost"], value_name="Count")
    console = Console()
    console.print("Estimating poly turnover... ")
    for idx, (name,group) in zip(track(range(len(list(df.groupby(['Mouse', 'PolyID','Dataset','Cage']).groups))),description='Turnover estimation...'),df.groupby(
            ["Mouse", "PolyID", "Dataset",'Cage'])):  ### index in the group object is the same as in the original dataframe
        presAbsDict = {day: (1 if day in group["Day"].values else 0) for day in days}
        presAbsList = [1 if day in group["Day"].values else 0 for day in days]
        timeIndex = 0
        for day in days:
            if day in group["Day"].tolist():
                if timeIndex == 0:  ## this would be the first time point that the poly appears, as days are ordered.
                    # df.loc[ group.loc[group['Day'] == day].index, "new"] = 1
                    ## distinction between ancestral and evolved for the first time point - if ancestral do NOT count it as new but only in cumseen
                    if group.Ancestry.values[0] == 'evolved':
                        turn = pd.concat([turn,pd.DataFrame({"Cage":[name[3]],"Dataset":[name[2]],"Mouse": [name[0]], "Day": [day], "PolyID": [name[1]],
                                                     "Synonymy": [group["Synonymy"].values[0]],
                                                     "new": [1], "seen": [0], "justseen": [0], "lost": [0],
                                                     "justlost": [0], "cumseen": [1]})])
                    else:
                        turn = pd.concat([turn,pd.DataFrame({"Cage":[name[3]],"Dataset":[name[2]],"Mouse": [name[0]], "Day": [day], "PolyID": [name[1]],
                                                     "Synonymy": [group["Synonymy"].values[0]],
                                                     "new": [0], "seen": [0], "justseen": [0], "lost": [0],
                                                     "justlost": [0], "cumseen": [1]})])
                    timeIndex += 1
                else:
                    if sum([presAbsDict[x] for x in presAbsDict if x < day]) >= 1 and presAbsDict[
                        days[days.index(day) - 1]] != 1:  ## has been seen at least one time before
                        # df.loc[(df["Day"] == day) & (df[bioUnit] == name[0]) & (df["PolyID"] == name[1]), "seen"] = 1
                        turn = pd.concat([turn,pd.DataFrame({"Cage":[name[3]],"Dataset":[name[2]],"Mouse": [name[0]], "Day": [day], "PolyID": [name[1]],
                                                         "Synonymy": [group["Synonymy"].values[0]],
                                                         "new": [0], "seen": [1], "justseen": [0], "lost": [0],
                                                         "justlost": [0], "cumseen": [1]})])
                    else:  ## the previous timepoint
                        # df.loc[ (df["Day"] == day) & (df[bioUnit] == name[0]) & (df["PolyID"] == name[1]), "justseen"] = 1
                        turn = pd.concat([turn, pd.DataFrame({"Cage":[name[3]],"Dataset":[name[2]], "Mouse": [name[0]], "Day": [day], "PolyID": [name[1]],
                                                         "Synonymy": [group["Synonymy"].values[0]],
                                                         "new": [0], "seen": [1], "justseen": [1], "lost": [0],
                                                         "justlost": [0], "cumseen": [1]})])
            else:  ## if day not in group["Day"] - lost or not seen yet. this means there is no entry in the dataframe
                if sum([presAbsDict[x] for x in presAbsDict if x < day]) >= 1 and presAbsDict[
                    days[days.index(day) - 1]] != 1:  ## it has been seen, but lost
                    # df.loc[(df["Day"] == day) & (df[bioUnit] == name[0]) & (df["PolyID"] == name[1]) , "lost"] = 1
                    turn = pd.concat([turn, pd.DataFrame({"Cage":[name[3]],"Dataset":[name[2]], "Mouse": [name[0]], "Day": [day], "PolyID": [name[1]],
                                                     "Synonymy": [group["Synonymy"].values[0]],
                                                     "new": [0], "seen": [0], "justseen": [0], "lost": [1],
                                                     "justlost": [0], "cumseen": [1]})])
                if sum([presAbsDict[x] for x in presAbsDict if x < day]) >= 1 and presAbsDict[
                    days[days.index(day) - 1]] == 1:  ## lost but seen in the previous timepoint

                    # df.loc[ (df["Day"] == day) & (df[bioUnit] == name[0]) & (df["PolyID"] == name[1]), "justlost"] = 1
                    turn = pd.concat([turn, pd.DataFrame({"Cage":[name[3]],"Dataset":[name[2]], "Mouse": [name[0]], "Day": [day], "PolyID": [name[1]],
                                                     "Synonymy": [group["Synonymy"].values[0]],
                                                     "new": [0], "seen": [0], "justseen": [0], "lost": [1],
                                                     "justlost": [1], "cumseen": [1]})])
                else:  ## has not been seen yet
                    turn = pd.concat([turn, pd.DataFrame(
                        {"Cage":[name[3]],"Dataset":[name[2]], "Mouse": [name[0]], "Day": [day], "PolyID": [name[1]],
                         "Synonymy": [group["Synonymy"].values[0]],
                         "new": [0], "seen": [0], "justseen": [0], "lost": [0], "justlost": [0], "cumseen": [0]})])
    # pickle.dump(turn, open(currentDir+"turn.obj",'wb'))

    # new, seen, justseen, lost, just lost
    aggregation = {"new": "sum", "seen": "sum", "justseen": "sum", "lost": "sum", 'justlost': "sum", "cumseen": "sum"}
    turnDf = turn.groupby(["Dataset","Mouse", "Day", "Cage"]).agg(aggregation).reset_index()
    #turnDf["Mouse"] = turnDf["Mouse"].astype(int)
    turnDf["Day"] = turnDf["Day"].astype(int)
    turnDf = pd.merge(df.drop_duplicates("SampleID")[["Dataset","Mouse", "Day", "PolyCount","Cage"]], turnDf, on=["Mouse", "Day","Cage"],
                      how="right", suffixes=("_y", ""))
    turnDf.drop(list(turnDf.filter(regex="_y")), axis=1, inplace=True)
    # g = sns.catplot(x="variable", y="Count",col=bioUnit, hue="Count", data=turnDf, saturation=.5,kind="bar", ci=None, aspect=.6)
    # (g.set_axis_labels("Day", "Count").set_titles("{col_name} {col_var}"))
    # rows = 2
    # columns = 3
    # fig = plt.figure()
    # # axis = plt.sublot(rows, columns index)
    #
    # barwidth = 0.3
    # plotIndex = 1
    # colors = ["red", "pink", "black", "blue", "cyan", "green", "orange"]
    categories = ["justseen", "seen", "new", "justlost", "lost", 'cumseen', "PolyCount"]
    # print(turnDf.columns)

    # for name, group in turnDf.groupby("Mouse"):
    #     y_pos = []
    #     plt.subplot(rows, columns, plotIndex)
    #     barIndex = 0
    #     for cat in categories:
    #         bars = group.sort_values("Day")[cat]
    #         y_pos = [x + (barIndex * barwidth) for x in np.arange(3 * len(days), step=3)]
    #         plt.bar(y_pos, bars, width=barwidth, color=colors[barIndex], label=cat)
    #         barIndex += 1
    #     plt.xticks([x + (len(categories) * barwidth / 3) for x in np.arange(3 * len(days), step=3)], days)
    #     plotIndex += 1
    # if plot: plt.show()
    meltedTurn = pd.melt(turnDf, id_vars=["Dataset","Day", "Mouse","Cage"], var_name="Category", value_vars=categories,
                         value_name="value")
    meltedTurn.replace({"Category": {"PolyCount": "Count for this timepoint",
                                     "new": "De novo",
                                     "justseen": "Maintained from previous timepoint",
                                     "seen": "Detected in this timepoint and previously",
                                     "justlost": "Lost from previous timepoint",
                                     "lost": "Detected previously but not present",
                                     "cumseen": "All polymorphisms detected"}}, inplace=True)

    meltedTurn.to_csv(outdir+""+NameTag+"TurnMeltedToR.csv")

    # pal = ["#ea0064", "#606fd0", "#31874d", "#505d54", "#f86b39", "#cab04c"]
    # g = sns.FacetGrid(data=meltedTurn, sharex=True, row="Category", hue="Mouse", palette=pal)
    # g = (g.map(plt.plot, "Day", "value").add_legend())
    # if plot: plt.show()

    # colorDict = {k: v for k, v in zip([x for x in range(1, 7)], ["green", "red", "black", "blue", "purple", "orange"])}
    # fig = plt.figure()
    # plt.style.use("fivethirtyeight")
    # import matplotlib.gridspec as gridspec
    # gs = gridspec.GridSpec(6, 5)
    # for index, cat in enumerate(categories):
    #     ax = plt.subplot(gs[index, :])
    #     ax.set_xticks(group.sort_values("Day")["Day"])
    #     ax.set_title(cat, fontdict={'fontsize': 8})
    #     for name, group in turnDf.groupby("Mouse"):
    #         ax.plot(group.sort_values("Day")["Day"], group.sort_values("Day")[cat], c=colorDict[int(name)], linewidth=1)
    #
    # from matplotlib.lines import Line2D
    # legend_elements = [Line2D([0], [0], color=colorDict[k], label=str(k)) for k in colorDict]
    # fig.legend(handles=legend_elements, loc='lower center', ncol=3)
    # if plot: plt.show()

    """ Estimate the empirical cumulative distribution.
     it is the number of occurences of data <= x  divided by len(data) 
     def ecdf(x):
        xs = np.sort(x)
        ys = np.arange(1, len(xs)+1) / float(len(xs))
        return xs, ys"""

    #plt.close("all")

    # from statsmodels.distributions.empirical_distribution import ECDF
    # ecdf = ECDF(turnDf["new"])
    # plt.plot(ecdf.x, ecdf.y)
    # plt.show()

    from scipy.stats import poisson
    mean, variance = np.mean(turnDf["new"]), np.var(turnDf["new"])
    print("Mean: {}, Variance:{}".format(mean, variance))
    Max, Min = np.max(turnDf["new"]), np.min(turnDf["new"])
    print("Min:{}, Max: {}".format(Min, Max))
    bins = np.linspace(Min, Max, 50)
    # poisson_values = poisson.rvs(mean, size=1000)
    # negBinom_values = scipy.stats.nbinom(m)
    # x = np.arange(len(turnDf["new"].unique()))
    # params = poisson.fit(turnDf["new"])
    # pdf_fitted = poisson.pdf(x, *params[:-2], loc=params[-2], scale=params[-1])
    # original = poisson.pdf(x)
    # x = np.arange(poisson.ppf(0.01, mean),poisson.ppf(0.99, mean))
    # plt.plot(x, poisson.pmf(x, mean), 'bo', ms=8, label='poisson pmf')
    # plt.hist(turnDf["new"], bins = bins, density=True )
    # plt.plot(x, pdf_fitted, '-r', x,original, 'b-')
    # plt.show()
def testNameTag(NameTag=""):
    if NameTag != "": NameTag +="_"
    print(NameTag)
def plot_maintainance(df, plot=True):
    days = df["Day"].drop_duplicates().sort_values().tolist()
    maintainDf = pd.DataFrame()
    for name, group in df.groupby(["Mouse", "PolyID"]):
        presAbsList = [1 if day in group["Day"].values else 0 for day in days]
        cnt = 0
        for x, y in zip(presAbsList, days):
            if x == 0:
                if cnt > 0:
                    maintainDf = maintainDf.append(pd.DataFrame({"Day": [y], "Mouse": [name[0]], "MainScore": [cnt]}),
                                                   ignore_index=True)
                cnt = 0
            else:
                cnt += 1
    plt.close('all')
    mainValues = np.arange(1, 6)
    counts = df["Maintenance"].value_counts()
    print(counts)
    print(df["Maintenance"].unique())
    # plt.hist(maintainDf["MainScore"], bins=mainValues)
    plt.title("Number of consequetive appearances of polymorphisms.")
    plt.bar(mainValues, counts.values, edgecolor="black")
    plt.xticks(mainValues)
    plt.show()



def most_mutated_genes(df, variable=False):
    df = df.loc[df['GeneName'] != 'Intergenic']
    descr_dict = df.drop_duplicates("GeneName").set_index('GeneName')["cddDescription"].to_dict()
    if not variable:
       df = df.groupby("GeneName").size().reset_index().rename(columns={0:'counts'}).sort_values('counts',ascending=False)
    else:
        df = df.groupby([variable, "GeneName"]).size().reset_index().rename(columns={0:'counts'}).sort_values('counts',ascending=False)
    df["cddDescription"] = df['GeneName'].map(descr_dict)
    return df

def cog_enrichment_by_group(df,variable="Day", NameTag=""):
    if NameTag != "": NameTag += "_"
    ## COG family enrichment in polys
    cogDf = pd.DataFrame()
    for name, group in df.groupby(["SampleID","CogFamily"]):
        Counts = group.shape[0]
        CogFamilySize = group["CogFamilySize"].astype(float).values[0]
        polymorphic_percent = (Counts*100) / CogFamilySize
        cogDf = cogDf.append(pd.DataFrame({variable:[group[variable].values[0]],
                                           "SampleID":[name[0]],
                                           "CogFamily":[name[1]],
                                           "PolymorphicPercent":[round(polymorphic_percent,6)]}))
    cogDf = cogDf.groupby([variable,'CogFamily']).agg({'PolymorphicPercent':'mean'}).reset_index()
    cogDf.to_csv("tmp/"+NameTag+"CogFamilyPolymorpicPercent.csv", index=False)

def compartment_polys(df):
    print("All polys detected")
    print(df.loc[df.Day == 28].groupby('SampleType').size().reset_index().rename(columns={0:'PolyCount'}))
    print("Compartment poly catalog")
    print(df.loc[df.Day == 28].drop_duplicates(["PolyID",'SampleType']).groupby('SampleType').size().reset_index().rename(columns={0:'PolyCount'}))

def compartment_enriched_polys(df, howmany=30):
    df = df.loc[df.Day == 28]
    compartment_enriched = []
    ## get samples that have high sampletype prevalence but low prevalence in general..
    for name, group in df.sort_values(by=["SampleTypePrevalencePercent","Prevalence","SampleTypeMeanFreq"], ascending=[False,True,False]).groupby("SampleType"):
        #print('\n'+name)
        #print(group.head(howmany)[['GeneName', 'cddDescription', 'CogFamily']])
        #group['SampleType'] = [name for x in group.index]
        #compartment_enriched = compartment_enriched.append(group[['Prevale''GeneName', 'GeneDescription', 'CogName', "SampleType"]])
        this_compartment = group.head(howmany)['PolyID'].tolist()
        compartment_enriched += this_compartment

    return list(set(compartment_enriched))
def prepare_files_for_heatmap(df, NameTag = "", mode='invivo'):
    if NameTag != "": NameTag += "_"
    pf = pd.pivot(df, index="PolyID", columns='SampleID', values="Freq").fillna(0).reset_index()
    pf.to_csv("tmp/"+NameTag+"PolyFrequencyHeatmap.tsv",sep='\t',index=False)
    pf.merge(df.drop_duplicates('PolyID'), on='PolyID', how='left').to_csv('tmp/'+NameTag+'PolyFrequencyHeatmapMetaRows.tsv',sep='\t',index=False)
    if mode == 'invivo':
        df.drop_duplicates('SampleID')[["SampleID","SampleType",'Day','Run','Cage','Mouse']].to_csv('tmp/'+NameTag+'PolyFrequencyHeatmapMetaColumns.tsv',sep='\t',index=False)
    else:
        df.drop_duplicates('SampleID')[["SampleID","SampleType",'Day','Run','Mix','Replicate']].to_csv('tmp/'+NameTag+'PolyFrequencyHeatmapMetaColumns.tsv',sep='\t',index=False)

def LDA_filtering_distance_between_runs(Df, binary=True):
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA

    X = pd.get_dummies(Df.Run.values) ### create dummy arrays for the categorial variable
    y = Df.PolyID
    model = LDA().fit(X, y)
    LDA_dict = {k: v for k, v in zip(model.classes_, [x[0] for x in model.coef_])} ## create dict with PolyID and its respective LDA coefficient.
    Df["LDAcoef"] = Df.PolyID.map(LDA_dict)

    """ PcoA of samples based on variant pattern both as presence/absence and with frequencies"""
    outDf = pd.DataFrame()
    for threshold in np.arange(0.1, 3, 0.1):
        LDA_filter = (Df.LDAcoef.between(-threshold, threshold))
        totalPolys = Df.loc[LDA_filter].drop_duplicates("PolyID").shape[0]
        table = pd.pivot(data=Df.loc[LDA_filter], index='SampleID', columns='PolyID', values="Freq").fillna(0)
        binary_table = table.copy()
        if binary is True:
            binary_table[binary_table > 0] = 1
        import seaborn as sns
        from skbio.stats.distance import DistanceMatrix
        Distance_matrix = DistanceMatrix(scipy.spatial.distance.pdist(binary_table.values)).data
        if Distance_matrix.shape[0] > 1:
            dm = pd.DataFrame(index=binary_table.index.values, columns=binary_table.index.values,data=Distance_matrix)
            dm = dm.stack().reset_index().rename(columns={'level_0': 's1', 'level_1': 's2', 0: 'Distance'})
            RunDistance = []
            for i, row in dm.iterrows():
                if (str(row['s1']).startswith("5") and str(row['s2']).startswith("6")) or (str(row['s1']).startswith("6") and str(row['s2']).startswith("5")):  ### select only library pairs from different runs
                    RunDistance.append(row['Distance'])
            outDf = outDf.append(pd.DataFrame(
                {"RunDistanceMean": [np.mean(RunDistance)], 'RunDistanceStd': [np.std(RunDistance)],
                 "LDACoefFilter": [threshold],"PolyCount": [totalPolys]}))

    return Df, outDf

def run_LDA(Df, dependent_variable="Mix", independent_variable="PolyID"):
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
    X = pd.get_dummies(Df[dependent_variable].values) ### create dummy arrays for the categorial variable
    y = Df[independent_variable]
    model = LDA().fit(X, y)
    LDA_dict = {k: v for k, v in zip(model.classes_, [x[0] for x in model.coef_])} ## create dict with PolyID and its respective LDA coefficient.
    Df["LDAcoef"] = Df[independent_variable].map(LDA_dict)
    return Df

def make_pcoa(df, binary=True, color='Dataset', style=None, plot=True, keep_columns=None, NameTag="", metric='euclidean'):

    """Will export binary_pcoa.tsv file in the tmp of the project dir. keep columns should be a list of columns to keep in the export file"""
    # table = pd.pivot(data=df.loc[filtering], index='SampleID', columns='PolyID', values="Freq").fillna(0)
    if NameTag != "": NameTag += "_"
    table = pd.pivot(data=df, index='SampleID', columns='PolyID', values="Freq").fillna(0)
    table = table.loc[(table != 0).any(axis=1), (table != 0).any(axis=0)] ## make sure there are no zero rows and columns
    binary_table = table.copy()
    if binary is True:
        binary_table[binary_table > 0] = 1
    binary_pcoa, var1, var2 = pcoa_from_dataframe(binary_table, metric=metric)
    binary_pcoa.reset_index(inplace=True)
    if keep_columns:
        binary_pcoa = binary_pcoa.merge(df.drop_duplicates("SampleID")[keep_columns], on='SampleID', how='left')
    else:
        binary_pcoa = binary_pcoa.merge(df.drop_duplicates("SampleID"), on='SampleID', how='left')

    binary_pcoa = binary_pcoa.merge(df.loc[df.Chrom == 'NC_004663'].drop_duplicates("SampleID")[['SampleID','Coverage_mean']], on='SampleID', how='left')
    binary_pcoa.to_csv('tmp/'+NameTag+'binary_pcoa.csv', index=None)

    if plot is True:
        sns.scatterplot(data=binary_pcoa, x='PC1', y='PC2', hue=color, style=style, alpha=0.7,)
        #for i,row in binary_pcoa.loc[binary_pcoa.SampleID.isin([50986, 69384, 69397,50982])].iterrows():
        #    plt.text(row.PC1, row.PC2, row.SampleID, color='red')
        plt.show()
    print(f"{NameTag} | PC1 {round(var1*100,2)} %, PC2 {round(var2*100,2)} %")

def run_permanova(Df, distance_metric='euclidean', variable='Mix', binary=True):
    """Permutational Multivariate Analysis of Variance (PERMANOVA) is a non-parametric method that tests whether two or more groups of objects (e.g., samples) are significantly different
    based on a categorical factor. It is conceptually similar to ANOVA except that it operates on a distance matrix, which allows for multivariate analysis. PERMANOVA computes a pseudo-F statistic.
    Statistical significance is assessed via a permutation test. The assignment of objects to groups (grouping) is randomly permuted a number of times (controlled via permutations).
    A pseudo-F statistic is computed for each permutation and the p-value is the proportion of permuted pseudo-F statisics that are equal to or greater than the original (unpermuted) pseudo-F statistic.
    """
    from skbio.stats.distance import permanova
    import scipy
    from skbio.stats.distance import DistanceMatrix
 #### do i need to use some transformation ?? root, clr?
    table = pd.pivot(data=Df, index='SampleID', columns='PolyID', values="Freq").fillna(0)
    if binary is True:
        table[table > 0] = 1
    Distance_matrix = DistanceMatrix(scipy.spatial.distance.pdist(table.values, metric=distance_metric), ids=table.index)
    #print("Running PERMANOVA...")
    #factors =["infantID","Age","SampleType" ]
    ## create new combinational factors
    #PatientplusAge = MasterDf['DM']+MasterDf['Age'].astype(str)
    results = pd.DataFrame()
    meta = Df.drop_duplicates("SampleID").merge(table.reset_index(), on='SampleID', how='right').set_index("SampleID")
    perm = permanova(Distance_matrix,meta[variable], permutations=1000)
    print(perm)
    thisDf  = pd.DataFrame(dict(zip(["F-statistic", "p-value", 'factor'],[[perm['test statistic']], [perm['p-value']], [variable]])))
    results = results.append(thisDf)
    return results.sort_values(['p-value','F-statistic'])


def make_random_forest_model(df, variable="Mix", binary=False, oneHotEncode=False):
    """Random forest analysis"""
    """We need to separate data into the features: data we use for the prediction and target: the data we want to predict."""

    table = pd.pivot(data=df, index='SampleID', columns='PolyID', values="Freq").fillna(0)
    if binary is True:
        table[table > 0] = 1
    features = table
    featureColumns = table.columns.values ### Polys
    target_col = table.merge(df.drop_duplicates("SampleID")[["SampleID", variable]].set_index("SampleID"), on='SampleID',how='left')[variable]
    if oneHotEncode is True:
        target_onehotEncode = pd.get_dummies(target_col)
        target = np.array(target_onehotEncode)
    else:
        target = np.array(target_col)
    print(target)
    from sklearn.model_selection import train_test_split
    #train_features, test_features, train_labels, test_labels =train_test_split(features, target, test_size=10, train_size=0.9, random_state=95)
    #print('Training Features Shape:', train_features.shape)
    #print('Training Labels Shape:', train_labels.shape)
    #print('Testing Features Shape:', test_features.shape)
    #print('Testing Labels Shape:', test_labels.shape)

    """ Train the model"""
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.ensemble import RandomForestClassifier
    # instantiate model with 1000 decision trees
    #rf = RandomForestRegressor(n_estimators=10000, random_state=95)
    rf = RandomForestClassifier(n_estimators=1000, random_state=95)
    # Train the model on training data

    accuracies = []
    # for feature in featureColumns:
    #     newfeatures= features.drop([feature], axis =1)
    #     newfeatures = newfeatures.astype(float).values
    #     rf.fit(newfeatures,target)
    #     accuracy = rf.score(newfeatures, target)
    #     accuracies.append(accuracy)
    # ##plot accuracy values
    # ind = np.arange(len(accuracies))
    # plt.bar(ind, accuracies, orientation='vertical')
    # plt.xticks(ind, featureColumns, rotation='vertical')
    # plt.ylabel("Mean Accuracy");plt.xlabel("SNP")
    # plt.title("Mean Accuracy removing SNP")
    # plt.show()

    rf.fit(features.astype(float), target)


    ## Make predictions
    predictions = rf.predict(features)
    print(predictions)
    #errors = sum([1.0 for x,y in zip(target, predictions) if x!=y])
    errors = np.sum(np.not_equal(target,predictions))
    accuracy = 100 - (100* (errors/len(predictions)))
    print ("This model's accuracy {}%".format(accuracy))


    importances = list(rf.feature_importances_)
    print ("This is the score (mean accurary): {}".format(rf.score(features, target)))
    feature_importances = sorted([(feature,round(importance,4)) for feature,importance in zip(featureColumns, importances) if importance >= 0] , key= lambda x:x[1],reverse=True)
    #[print("Polys: {:20} Importance: {}".format(*pair)) for pair in feature_importances if pair[1] > 0]
    imp_df = pd.DataFrame({'feature':[x[0] for x in feature_importances], 'importance':[x[1] for x in feature_importances]})

    ##plot importance values
    # ind = np.arange(len(feature_importances))
    # plt.bar(ind, [x[1] for x in feature_importances], orientation='vertical')
    # plt.xticks(ind,[x[0] for x in feature_importances], rotation='vertical')
    # plt.ylabel("Importance");plt.xlabel("Variable")
    # plt.title("Variable Importance")
    # plt.show()
    return imp_df

def run_logistic_regression(df, binary=False):
    import statsmodels.api as sm
    from sklearn.linear_model import LogisticRegression
    table = pd.pivot(data=df, index='SampleID', columns='PolyID', values="Freq").fillna(0)
    binary_table = table.copy()
    if binary is True:
        binary_table[binary_table > 0] = 1
    meta = df.drop_duplicates("SampleID").merge(table.reset_index(), on='SampleID', how='right').set_index("SampleID")
    dep = 'Mix'
    for x in table.columns:
        X = pd.get_dummies(table[x].values)
        X = sm.add_constant(X)
        y = meta[dep].values
        model = sm.Logit(y,X, missing='drop')
        res = model.fit_regularized()
        print(res.summary())

####### Annotation related functions
def gbff_to_df(file_path, protein_file_path=None):
    """
    For each feauture:
    type: rRNA - tRNA - CDS, tmRNA, regulatory, ncRNA
    location: [XX:XX](+) (- if transcribed from the reverse strand) - UPDATE - now location is XX...XX
    qualifiers: locus_tag, codon_start, transl_table, product, protein_id, db_xref, translation (they have values lists)
    """
    if protein_file_path:
        prot_df = fasta_to_df(protein_file_path)
    genbank_df = pd.DataFrame()
    prot_expr = re.compile(r'(?<=RefSeq:).+')
    for seq_record in SeqIO.parse(file_path, "genbank"):
        for feature in seq_record.features:
            if feature.type in ["CDS",'tRNA','rRNA', 'ncRNA','tmRNA', 'regulatory']:
                protein_id = "" ### need to search for protein id  in inference key cause for some entries the protein_id key is empty
                translation = ""
                pseudo = 'no'
                if feature.qualifiers:
                    if feature.type == "CDS":
                        if 'protein_id' not in feature.qualifiers:
                            if re.search(prot_expr, feature.qualifiers['inference'][0]):
                                protein_id = re.search(prot_expr, feature.qualifiers['inference'][0]).group(0)

                            if protein_file_path and protein_id in prot_df.Chrom.values:
                                translation = prot_df.loc[prot_df.Chrom == protein_id]['seq'].values[0]

                        else:
                            protein_id = feature.qualifiers['protein_id'][0]
                        if not protein_file_path and 'translation' in feature.qualifiers:
                            translation = feature.qualifiers['translation'][0] ## some entries have protein id but not translation, so use the protein.faa instead
                        elif protein_file_path and protein_id in prot_df.Chrom.values:
                            translation = prot_df.loc[prot_df.Chrom == protein_id]['seq'].values[0]

                    if 'pseudo' in feature.qualifiers:
                        pseudo = 'yes'

                #start, end, orientation =list(re.match(re.compile(r'\[<*(\d+):>*(\d+)\]\(([\+\-])\).*'), str(feature.location)).groups())
                if protein_id != '':
                    this_df = pd.DataFrame({"Chrom":seq_record.id,'type':[feature.type], 'start':[feature.location.start], 'end':[feature.location.end],'strand':[feature.location.strand], 'codon_start':feature.qualifiers['codon_start'],
                                                             'old_locus_tag':[""],'locus_tag':feature.qualifiers['locus_tag'],'product':feature.qualifiers['product'],
                                                             "protein_id":[protein_id], "translation":[translation], 'pseudo':[pseudo], 'DNA':[feature.location.extract(seq_record).seq]})
                    if 'old_locus_tag' in feature.qualifiers:
                        this_df['old_locus_tag'] = feature.qualifiers['old_locus_tag'][0]

                    genbank_df = pd.concat([genbank_df,this_df], ignore_index=True)

                else:
                    if feature.type == 'regulatory':
                        this_df = pd.DataFrame(
                            {'Chrom': [seq_record.id], 'type': [feature.type], 'start': [feature.location.start],
                             'end': [feature.location.end], 'strand': [feature.location.strand], 'codon_start': [""],
                             'old_locus_tag': [""], 'locus_tag': feature.qualifiers['db_xref'],
                             'product': feature.qualifiers['note'],
                             "protein_id": [""], "translation": [""], 'pseudo':[pseudo],'DNA':[feature.location.extract(seq_record).seq]})
                        if 'old_locus_tag' in feature.qualifiers:
                            this_df['old_locus_tag'] = feature.qualifiers['old_locus_tag'][0]
                        genbank_df = pd.concat([genbank_df, this_df], ignore_index=True)

                    else:

                        this_df = pd.DataFrame({'Chrom':[seq_record.id],'type':[feature.type], 'start':[feature.location.start], 'end':[feature.location.end],'strand':[feature.location.strand], 'codon_start':[""],
                                                                 'old_locus_tag':[""],'locus_tag':feature.qualifiers['locus_tag'],'product':feature.qualifiers['product'],
                                                                 "protein_id":[""], "translation":[""], 'pseudo':[pseudo],'DNA':[feature.location.extract(seq_record).seq]})
                        if 'old_locus_tag' in feature.qualifiers:
                            this_df['old_locus_tag'] = feature.qualifiers['old_locus_tag'][0]
                        genbank_df = pd.concat([genbank_df,this_df], ignore_index=True)
    return genbank_df

def dna_fasta_from_genbank(genbank_file):
    f = gbff_to_df(genbank_file)
    #for i, row in f.iterrows():
    #    print(row['locus_tag'])
    #    print(row['product'])
    sequences = [SeqIO.SeqRecord(id=row['locus_tag'], name='',description='', seq=row['DNA']) for i,row in f.iterrows()]

    with open(genbank_file.replace(".gbff",'_genes.fa'),'w') as output:
        SeqIO.write(sequences, output, 'fasta')
def cds_from_fna_to_df(file_path):
    cds_df = pd.DataFrame()
    for rec in SeqIO.parse(file_path,'fasta'):
        regex = re.compile(r'.*\[locus_tag=(.+)\]\s+\[protein=(.+)\]\s+\[[^=]+=(.+)\]\s+\[location=.+\]')
        locus_tag, protein, protein_id = re.match(regex, rec.description).groups()
        cds_df = pd.concat([cds_df,pd.DataFrame({'locus_tag':[locus_tag],'protein':[protein], "protein_id":[protein_id],'cds':["".join(rec.seq)]})])
    return cds_df

def fasta_to_df(file_path):
    fasta_df = pd.DataFrame()
    for rec in SeqIO.parse(file_path,'fasta'):
        fasta_df = pd.concat([fasta_df,pd.DataFrame({'Chrom':[rec.id],'seq':["".join(rec.seq)], 'Description':[rec.description]})])
    return fasta_df

########

def copy_number_variation_from_metagenome(bamfile, df):
    """
    :param bamfile: path to bamfile
    :param df: a parsed genbank file with gbff_to_df()
    :return: returns a pandas dataframe with coverage per locus_tag
    """
    #console = Console()
    #console.print("[b]Parsing BAM file...[/b]")
    sampleID = os.path.basename(os.path.dirname(bamfile)) ## expects SampleName/reference.bam
    if not os.path.isfile(bamfile+'.bai'): os.system(f'samtools index -b {bamfile}')
    bamfile = pysam.AlignmentFile(bamfile, 'rb')
    out_df = pd.DataFrame()
    for i, row in df.iterrows(): ## dataframe with all loci from genbank file
        locus = np.array(bamfile.count_coverage(contig=row.Chrom.rstrip(".1"), start=row.start, stop=row.end, quality_threshold=30))
        locus_coverage = np.sum(locus) / locus.shape[1]
        out_df = out_df.append(pd.DataFrame({"locus_tag":row['locus_tag'], "start":row.start, "end":row.end, 'coverage':[locus_coverage], "SampleID":[sampleID]}))
    return out_df
def find_best_number_of_clusters(CountsDf, metric='euclidean', exportFolder="",maxCluster=100, title=""):
    CountsDf = CountsDf[CountsDf.sum(axis=0) > 0] ##remove zero rows
    pcoaDf, v1, v2 = pcoa_from_dataframe(CountsDf, metric=metric)
    if maxCluster >= pcoaDf.shape[0]:
        maxCluster = pcoaDf.shape[0]  / 2
    from sklearn.metrics import calinski_harabasz_score
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score, silhouette_samples
    import matplotlib.cm as cm
    ch_scores = []
    sil_scores = []
    for n_clusters in range(2,maxCluster):
        #fig, (ax1, ax2) = plt.subplots(1,2)
        #fig.set_size_inches(18,7)
        #ax1.set_ylabel([0, len(pcoaDf.values) + (n_clusters+1)*10])
        kmeans = KMeans(n_clusters=n_clusters, random_state=9).fit(pcoaDf.values)
        kmean_clusters = kmeans.labels_
        ch_score = calinski_harabasz_score(pcoaDf.values, kmean_clusters)
        sample_silhoutte_values = silhouette_samples(pcoaDf.values, kmean_clusters)
        sil_score = silhouette_score(pcoaDf.values, kmean_clusters)
        sil_scores.append(sil_score)
        ch_scores.append(ch_score)
        y_lower = 10
        # for i in range(n_clusters):
        #     ith_cluster_silhouette_values =  sample_silhoutte_values[kmean_clusters == i]
        #     ith_cluster_silhouette_values.sort()
        #     size_cluster_i = ith_cluster_silhouette_values.shape[0]
        #     y_upper = y_lower + size_cluster_i
        #     color = cm.nipy_spectral(float(i) / n_clusters)
        #     ax1.fill_betweenx(np.arange(y_lower, y_upper),
        #                   0, ith_cluster_silhouette_values,
        #                   facecolor=color, edgecolor=color, alpha=0.7)
        #
        #     ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
        #     y_lower = y_upper + 10
        # # The vertical line for average silhouette score of all the values
        # ax1.axvline(x=sil_score, color="red", linestyle="--")
        # ax1.set_yticks([])  # Clear the yaxis labels / ticks
        # ax1.set_title("The silhouette plot for the various clusters.")
        # ax1.set_xlabel("The silhouette coefficient values")
        # ax1.set_ylabel("Cluster label")
        # colors = cm.nipy_spectral(kmean_clusters.astype(float) / n_clusters)
        # ax2.scatter(pcoaDf.values[:, 0], pcoaDf.values[:, 1], marker='.', s=30, lw=0, alpha=0.7,
        #         c=colors, edgecolor='k')
        # # Labeling the clusters
        # centers = kmeans.cluster_centers_
        # # Draw white circles at cluster centers
        # ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
        # c="white", alpha=1, s=200, edgecolor='k')
        #
        # for i, c in enumerate(centers):
        #     ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
        #             s=50, edgecolor='k')
        # ax2.set_title("The visualization of the clustered data.")
        # ax2.set_xlabel("PC1")
        # ax2.set_ylabel("PC2")
        # plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
        #           "with n_clusters = %d" % n_clusters),
        #          fontsize=14, fontweight='bold')
        # if exportFolder == "":
        #     fig.savefig('/media/christos/ssd/work/Infants/figures/Silhouette_figures'+str(n_clusters)+'.svg')
        # else:
        #     fig.savefig(exportFolder+'Silhouette_figures' + str(n_clusters) + '.svg')
    #plt.show()
    plt.close()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot([str(i) for i in range(2,maxCluster)], np.log(ch_scores)/10.0, label="Calinski-Harabz score (log /10)")
    ax.plot([str(i) for i in range(2,maxCluster)], np.array(sil_scores), label="Sihlouette score")
    ax.legend()
    plt.title(title)
    plt.xlabel("Number of Clusters")
    plt.tight_layout()
    if exportFolder == "":
        fig.savefig('/run/media/christos/ssd/work/bacevo/invivo/tmp/Evaluate_clusters'+title+'.png')

    else:
        fig.savefig(exportFolder+'Evaluate_clusters'+title+'.png')
    #     fig.savefig(exportFolder +'number-of-clusters.png')
    plt.close()
    #plt.show()

#### clustering polys
def clustering_general(Df, metric='braycurtis', nclusters=4, eps=0.1, outfolder=''):
    MasterDf = pd.DataFrame()
    sns.set(style="whitegrid")
    plt.rcParams.update({'axes.labelsize': 'x-large',
                         'xtick.labelsize': 'large',
                         'ytick.labelsize': 'large'})

    pcoaDf, v1, v2 = pcoa_from_dataframe(Df, metric=metric)
    import sklearn.cluster as skc
    from sklearn.cluster import DBSCAN
    from sklearn import metrics
    kmeans = skc.KMeans(n_clusters=nclusters).fit(pcoaDf.values)
    # import sklearn.preprocessing.StandardScaler as StandardScaler
    from sklearn import preprocessing
    # import sklearn.preprocessing.StandardScaler as StandardScaler

    X = preprocessing.StandardScaler().fit_transform(pcoaDf.values)
    db = DBSCAN(eps=eps).fit(X)
    dbscan_clusters = DBSCAN().fit_predict(X)
    # print (dbscan_clusters)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    ##number of clusters in labels ignoring noise if present (-1 is the 'cluster' label for noisy samples (not core)
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    clusters = kmeans.labels_
    db_clusters = labels
    MasterDf = MasterDf.assign(PC1=pcoaDf.PC1)
    MasterDf = MasterDf.assign(PC2=pcoaDf.PC2)
    MasterDf = MasterDf.assign(kmCluster=clusters)  ### these are the clusters from kmeans
    MasterDf = MasterDf.assign(dbCluster=db_clusters)

    centroids = kmeans.cluster_centers_
    #pcoaDf = pcoaDf.assign(Cluster=clusters)  ### these are the clusters from kmeans

    ### DBSCAN plotting
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    unique_labels = set(labels)  ### dbscan labels
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    for k, col in zip(unique_labels, colors):
        if k == -1:
            ## Black used for noise
            col = [0, 0, 0, 1]
        class_member_mask = (labels == k)
        xy = X[class_member_mask & core_samples_mask]
        db_centroid_x, db_centroid_y = np.mean(xy[:, 0]), np.mean(xy[:, 1])

        ax.scatter(xy[:, 0], xy[:, 1], marker='o', color=tuple(col), edgecolor='k', label=k, s=90)
        xy = X[class_member_mask & ~core_samples_mask]

        ax.scatter(xy[:, 0], xy[:, 1], marker='o', color=tuple(col), edgecolor='k', s=30)
        #if k != -1:
        #    ax.annotate(MostAbundantDBscan.loc[k, 'label'], (db_centroid_x+0.05, db_centroid_y+0.05))
    handles, labels = ax.get_legend_handles_labels()
    legend1 = ax.legend(handles, labels, title=r'dbscan Cluster', fontsize=14, frameon=True, bbox_to_anchor=(1.301, 1))
    plt.title('PCoA-BrayCurtis: Estimated number of clusters: %d' % n_clusters_)
    fig.text(0.03, 0.5, "PC2 [{0:.0%}]".format(v2), rotation='vertical', va='center', fontsize=16)
    fig.text(0.4, 0.03, "PC1 [{0:.0%} ]".format(v1), va='center', fontsize=16)
    # fig.tight_layout()
    plt.subplots_adjust(right=0.8)
    plt.show()
    #fig.savefig(outfolder + "/Dbscan_clusters.png")
    #fig.savefig(outfolder + "/Dbscan_clusters.svg")

    # plt.close()

    fig = plt.figure(figsize=(10, 8))
    sizefactor = 2
    # pal = dict(Stool=sns.xkcd_rgb["brown"], Skin=sns.xkcd_rgb["peach"], Oral=sns.xkcd_rgb["salmon"])
    ax1 = fig.add_subplot(111)
    grouped = pcoaDf.groupby("Cluster")
    pal = sns.color_palette("Set2", grouped.ngroups)
    Index = 0
    for name, group in grouped:
        for x, y in zip(group["PC1"], group["PC2"]):
            ax1.plot([centroids[int(name)][0], x], [centroids[int(name)][1], y], color='k', alpha=0.2)
        ax1.scatter(group["PC1"], group["PC2"], color=pal[Index], s=90, alpha=0.8, marker='o', edgecolor='k',
                    label=name)
        ax1.scatter(centroids[int(name)][0], centroids[int(name)][1], s=100, color="white", marker="^", edgecolor='k',
                    linewidths=1)
        # ax1.annotate(keySpeciesDict[name], (centroids[int(name)][0],centroids[int(name)][1]))
        #ax1.annotate(MostAbundant.loc[name, 'label'], (centroids[int(name)][0]+0.01, centroids[int(name)][1]+0.01))
        Index += 1
    handles, labels = ax1.get_legend_handles_labels()
    legend1 = ax1.legend(handles, labels, title=r'kmeans Cluster', fontsize=14, frameon=True, bbox_to_anchor=(1.305, 1))
    legend1.get_title().set_fontsize(14)
    fig.text(0.03, 0.5, "PC2 [{0:.0%}]".format(v2), rotation='vertical', va='center', fontsize=16)
    fig.text(0.4, 0.03, "PC1 [{0:.0%}]".format(v1), va='center', fontsize=16)
    # fig.suptitle(title)
    # fig.tight_layout()
    plt.subplots_adjust(right=0.8)
    plt.show()
    #fig.savefig(outfolder + "/Kmeans_clusters.png")
    #fig.savefig(outfolder + '/kmeans_clusters.svg')


def find_late_prevalent_polys(df, cutoff=2):
    """

    :param df: dataframe of FinalPolyTable.tsv
    :return: dataframe with PolyID colymn with all late (appearing at day 21 onwards) polys in at least three trajectories
    """
    out = pd.DataFrame()
    for poly, group in df.loc[df.SampleType == 'Feces'].drop_duplicates(['SampleID', 'PolyID']).groupby('PolyID'):
        cnt = 0

        for mouse, sub in group.groupby("Mouse"):


            if (sub.Day.sort_values().tolist() == [21,28] or sub.Day.sort_values().tolist() == [28]):
                cnt += 1
        if cnt > cutoff:
            out = pd.concat([out, pd.DataFrame({'PolyID': [poly]})])

    return out

def find_parallel_polys(df, cutoff=2):
    """

    :param df: dataframe of FinalPolyTable.tsv
    :return: dataframe with PolyID colymn with all late (appearing at day 21 onwards) polys in at least three trajectories
    """
    out = pd.DataFrame()
    for poly, group in df.loc[df.SampleType == 'Feces'].drop_duplicates(['SampleID', 'PolyID']).groupby('PolyID'):
        parallel = len(group.Mouse.unique())
        if parallel > cutoff:
            out = pd.concat([out, pd.DataFrame({'PolyID': [poly]})])

    return out




def find_early_persistent_polys(df, cutoff=2):
    """

    :param df: dataframe of FinalPolyTable.tsv
    :return: dataframe with PolyID column with all early and persistent poyls in at least three trajectories
    """
    out = pd.DataFrame()
    for poly, group in df.loc[df.SampleType == 'Feces'].drop_duplicates(['SampleID','PolyID']).groupby("PolyID"):
        cnt = 0
        for mouse, sub in group.groupby("Mouse"):
            if sub.shape[0] >= 4: ### present in at least 4 time points
                cnt += 1
        if cnt > cutoff: ## present in at least 3 trajectories
            out = pd.concat([out, pd.DataFrame({"PolyID":[poly]})])
    return out

def find_individual_peristent_polys(df):
    """

    :param df: dataframe of FinalPolyTable.tsv
    :return: dataframe with PolyID column with all early and persistent polys that appear only in one trajectory
    """
    out = pd.DataFrame()
    for poly, group in df.loc[df.SampleType == 'Feces'].drop_duplicates(['SampleID','PolyID']).groupby("PolyID"):
        cnt = 0
        for mouse, sub in group.groupby("Mouse"):
            if sub.shape[0] >= 4: ### present in at least 4 time points
                cnt += 1
        if cnt == 1: ## present in only one trajectory
            out = pd.concat([out, pd.DataFrame({"PolyID":[poly]})])
    return out


def add_trajectories_appearing(df):
    """

    :param df: dataframe of FinalPolyTable.tsv
    :return: dataframe with PolyID colymn and number of trajectories detected (at least once)
    """
    out = pd.DataFrame()
    df['TrajectoriesDetected'] = [0 for x in df.index]
    for poly, group in df.loc[df.SampleType == 'Feces'].drop_duplicates(['SampleID', 'PolyID']).groupby('PolyID'):
        df.loc[df.PolyID == poly, "TrajectoriesDetected"] = group.Mouse.unique().size

    return df