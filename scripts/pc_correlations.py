import time
import multiprocessing as mp
import scipy.stats
from bacevo_fun import *

"""      Detect correlated polymorphisms            """

def make_poly_correlation_comparisons(name,group):
    """
    Will only calculate correlation of samples with at max 1 missing sample in their trajectory
    :param name:
    :param group:
    :return:
    """
    overwrite = args.overwrite
    if overwrite is True and os.path.isfile(f"tmp/{name}_Poly_correlation.csv"):
        os.system(f'rm tmp/{name}_Poly_correlation.csv')

    if not os.path.isfile(f"tmp/{name}_Poly_correlation.csv"):  ## check if file already available and skip

        cov_df = pd.DataFrame()
        N_polys = group.drop_duplicates("PolyID").shape[0]

        #print(f'Working on {N_polys} polys ({N_polys*(N_polys-1)/2} combinations) for Mouse {name} ')


        for idx, (x, y) in zip(track(range(len(list(itertools.combinations(group.PolyID.unique(), 2)))),
                                     description=f'Poly comparisons for Mouse {name}', total=N_polys*(N_polys-1)/2),
                               itertools.combinations(group.PolyID.unique(), 2)):
            days = [3, 7, 14, 21, 28]
            freq_x = []
            freq_y = []
            missing_samples_x = 0
            missing_samples_y = 0
            for day in days:
                if group.loc[(group.PolyID == x) & (group.Day == day)].empty is False:
                    freq_x.append(group.loc[(group.PolyID == x) & (group.Day == day)]["Freq"].values[0])
                else:
                    freq_x.append(0)
                    missing_samples_x += 1
                if group.loc[(group.PolyID == y) & (group.Day == day)].empty is False:
                    freq_y.append(group.loc[(group.PolyID == y) & (group.Day == day)]["Freq"].values[0])
                else:
                    freq_y.append(0)
                    missing_samples_y += 1
            if not missing_samples_x > 1 and not missing_samples_y > 1:
                corr, p_value = scipy.stats.pearsonr(freq_x, freq_y)
            else:
                corr, p_value = 0,1
            cov_df = pd.concat([cov_df, pd.DataFrame({'polyA': [x], 'polyB': [y], 'cor': [corr], 'p-value':[p_value],'Mouse': [name]})],
                               ignore_index=True)
        cov_df.to_csv(f'tmp/{name}_Poly_correlation.csv', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest='threads', help='How many threads to be used. Default=2', type=int, default=2)
    parser.add_argument('-d', dest='workdir',help='Working directory.., Default is current working directory', default=os.getcwd())
    parser.add_argument('-p', dest='polytable', help='Path to a FinalPolyTable.tsv.', required=True)
    parser.add_argument('-m', dest='mode', help='Invitro or invivo', choices=['invitro','invivo'], default='invivo')
    parser.add_argument('-w', dest='overwrite', help='Overwrite existing files. Default is False', action='store_true', default=False)
    args = parser.parse_args()

    currentDir = args.workdir
    os.chdir(currentDir)
    if not os.path.isdir(currentDir+'/tmp'): os.mkdir(currentDir+'/tmp')
    calculate_correlation = True
    df = pd.read_csv(args.polytable,sep='\t', low_memory=False) ### this is the concatenation of bern_metagenomes and ploen FinalPolyTables with the ep_merge_tables.py script
    df = df.drop_duplicates(['SampleID','PolyID'])
    console = Console()

    #filtered_polys = df.loc[(df.MeanFreq > 0.2) & (df.MousePrevalence >= 2)]['PolyID'].tolist()
    #df = df.loc[df.PolyID.isin(filtered_polys)]


    cov_df = pd.DataFrame()

    start = time.perf_counter()
    #results = [pool.starmap(make_poly_correlation_comparisons, [ (name, group) for name, group in df[['SampleID', 'Freq','PolyID', "Mouse","Day"]].groupby("Mouse")])]
    if calculate_correlation is True:
        pool = mp.Pool(args.threads)
        if args.mode == 'invivo':
            ### use only fecal samples (time points)

            df = df.loc[(df.SampleType == 'Feces') & (df.MousePrevalencePercent >= 0.8) & (df.TotalCov.ge(100))] ###  4/5 prevalent in the trajectory
            print(f'Starting analysis for {df.drop_duplicates("Mouse").shape[0]} specimens..')
            correlations = pool.starmap(make_poly_correlation_comparisons,[(name,group) for name, group in df[['SampleID', 'Freq','PolyID', "Mouse","Day"]].groupby("Mouse")])
            pool.close()
            pool.join()


    os.chdir(os.getcwd()+"/tmp")
    cor_df = pd.concat([pd.read_csv(x) for x in glob.glob("*_Poly_correlation.csv")])
    cor_df.to_csv("Poly_correlation.csv", index=False)

    print('\n\nDone!')


