#!/apps/python3/3.9.0/bin/python3
import os

import pandas as pd

from bacevo_fun import *
from itertools import groupby
parser = argparse.ArgumentParser(description="Poly turnover estimation ")
parser.add_argument('-p', dest='poly', help='FinalPolyTable', required=True)
parser.add_argument('-o', dest='outdir', help='Output directory. Default = CWD.', default=os.getcwd())
parser.add_argument('-m', dest='mask', help='file with polys to mask. Expects a csv file with PolyID columnname.', required=True)

args = parser.parse_args()



if args.outdir == "":
    outdir = os.path.dirname(args.poly)
else:
    outdir = args.outdir
if not outdir.endswith("/"):
    outdir += "/"
df = pd.read_csv(args.poly,sep='\t', usecols =['Chrom','Mouse','Day', 'Freq','PolyID','Dataset','SampleID','PolyID', "TotalCov",'Coverage_mean','Ancestry', 'SampleType','AltCov'])
df = df.loc[df.SampleType == 'Feces']
keep_samples = df.loc[(df.Chrom == 'NC_004663') & (df['Coverage_mean'] >= 100), 'SampleID'].unique()
df = df.loc[(df.TotalCov.ge(100)) & (df.SampleID.isin(keep_samples))].drop_duplicates(['SampleID','PolyID'])

if args.mask:
    mask = pd.read_csv(args.mask)
    df = df[~df.PolyID.isin(mask.PolyID.unique())]


pr = pd.DataFrame()
for i, ((mouse,poly), group) in zip(track(range(len(df.drop_duplicates(['SampleID','PolyID']).groupby(["Mouse",'PolyID']).groups)), description='Analyzing trajectory..'),df.drop_duplicates(['SampleID','PolyID']).groupby(["Mouse",'PolyID'])):
    days = df.Day.sort_values().unique()
    presAbsDict = {day: (1 if day in group["Day"].values else 0) for day in days}
    day_dict = {0:3, 1:7, 2:14, 3:21, 4:28}
    traj = [1 if day in group["Day"].values else 0 for day in days]
    mouse_preva = sum(traj)
    g = groupby(traj, key=lambda x:x>0.0)
    mouse_persi= sum(max([list(s) for v, s in g if v > 0.0], key=len, default=0))
    pr = pd.concat([pr, pd.DataFrame(
        {'Mouse': [mouse], 'PolyID': [poly], 'TrajectoryPersistence': [mouse_persi], 'TrajectoryPrevalence':[mouse_preva]})])
pr.to_csv('~/bacevo/invivo/tmp/TrajPer.csv',index=False)

tr = pd.DataFrame()
### poly turnover 1. SNP count, 2. new mutations 3. cumulative number of mutations for each trajectory
for i, ((mouse,poly), group) in zip(track(range(len(df.drop_duplicates(['SampleID','PolyID']).groupby(["Mouse",'PolyID']).groups)), description='Analyzing trajectory..'),df.drop_duplicates(['SampleID','PolyID']).groupby(["Mouse",'PolyID'])):
    days = df.Day.sort_values().unique()
    presAbsDict = {day: (1 if day in group["Day"].values else 0) for day in days}
    day_dict = {0:3, 1:7, 2:14, 3:21, 4:28}
    traj = [1 if day in group["Day"].values else 0 for day in days]

    ## first time the
    first_time = 0
    #em_check = 0 ## this switches to 1 once poly is present in the trajectory (for cumulative number of different polys)
    c_index = 0
    for i,day in enumerate(traj):
        bf_check = 0
        ls_check = 0

        if day == 1 and i == 0 and group.Ancestry.values[0] == 'evolved': ## first timepoint
            #em_check = 1
            #c_index += 1
            tr = pd.concat([tr, pd.DataFrame({"PolyID":[poly],'Mouse':[mouse],'Day':[day_dict[i]], 'Dataset':[group.Dataset.values[0]],
                                              "Total":[1],'FirstEmerge':[1], 'Denovo':[1]})])
        elif day == 1 and i == 0 and group.Ancestry.values[0] == 'ancestral':
            tr = pd.concat([tr, pd.DataFrame({"PolyID":[poly],'Mouse':[mouse],'Day':[day_dict[i]], 'Dataset':[group.Dataset.values[0]],
                                              "Total":[1],'FirstEmerge':[0], 'Denovo':[0]})])
        # elif day == 0 and i == 0:
        #     tr = pd.concat([tr, pd.DataFrame({"PolyID":[poly],'Mouse':[mouse],'Day':[day_dict[day]], 'Dataset':[group.Dataset.values[0]],
        #                                       "Total":[0],'NewSpots':[0], 'Denovo':[0]})])

        elif day == 1 and i > 0 :
            #em_check = 1
            ## check if poly has been seen before
            if sum(traj[:i]) == 0:  ## not present before if this is zero
                bf_check = 1
            ## check if poly has been present in the previous timepoint, else is a de novo poly
            if traj[i-1] == 0:
                ls_check = 1
            tr = pd.concat([tr, pd.DataFrame({"PolyID":[poly],'Mouse':[mouse],'Day':[day_dict[i]], 'Dataset':[group.Dataset.values[0]],
                                              "Total":[1],'FirstEmerge':[bf_check], 'Denovo':[ls_check]})])

        else: ## here poly is absent

            tr = pd.concat([tr, pd.DataFrame({"PolyID":[poly],'Mouse':[mouse],'Day':[day_dict[i]], 'Dataset':[group.Dataset.values[0]],
                                              "Total":[0],'FirstEmerge':[0], 'Denovo':[0]})])

tr = tr.groupby(['Mouse', 'Day', 'Dataset']).agg({'Total':'sum', 'FirstEmerge':'sum', 'Denovo':'sum'}).reset_index()
tr.sort_values(['Mouse','Day'], inplace=True)
tr['CumDenovo'] = tr.groupby(['Mouse'])['Denovo'].transform(lambda x: x.cumsum())
tr['CumFirstEmerge'] = tr.groupby(['Mouse'])['FirstEmerge'].transform(lambda x: x.cumsum())
tr['CumTotal'] = tr.groupby(['Mouse'])['Total'].transform(lambda x: x.cumsum())

tr.to_csv("~/bacevo/invivo/tmp/Turnover.csv", index=False)


