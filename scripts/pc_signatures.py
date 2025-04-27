
from bacevo_fun import *


def compute_poly_clusters(name,df, trajectory_group):
    clusters = pd.DataFrame()
    c_index = 1
    #df.drop_duplicates(['polyA','polyB'], inplace=True)
    #df = df.loc[df.cor.ge(0.8)]
    df.sort_values('cor', ascending=False, inplace=True)
    polys = sorted(list(set(df.polyA.tolist() + df.polyB.tolist()))) #### the order might have an effect on the cluster !

    def create_tab(polys):
        tab = pd.DataFrame(columns=polys , index=polys, data= np.zeros(shape=(len(polys),len(polys))))
        for i,row in df.iterrows():
            if row.polyA != row.polyB:
                tab.loc[row.polyA, row.polyB] = row.cor
                tab.loc[row.polyB, row.polyA] = row.cor
        return tab

    init_tab = create_tab(polys)
    tab = init_tab.copy()
    cl = pd.DataFrame({x:[[x],[]] for x in polys})
    """First row stores polys, second row stores correlations"""
    c_index = 1
    for p in polys:
        pair, max_row = tab[p].idxmax(), tab[p].max()
        if max_row >= 0.8 and p != pair:
            cl.loc[0, pair].append(p)
            cl.loc[1, pair].append(max_row)
            cl.drop(columns=[p],inplace=True)
            if pair in polys:
                polys.remove(pair)
        tab.drop(columns=p, inplace=True)
        tab.drop(index=p, inplace=True)
    cl.rename(columns={y:f'CL-{name}-{x+1}' for x,y in enumerate(cl.columns.values)}, inplace=True)
    singletons = []
    cl_sizes = []
    for c in cl.columns:
        if len(cl.loc[0,c]) > 1:
            cl_sizes.append(len(cl.loc[0,c]))

        if len(cl.loc[0,c]) == 1:
            singletons.append(len(cl.loc[0,c]))
    print(f'Singletons: {len(singletons)}, clusters: {len(cl_sizes)} ({sum(cl_sizes)} polys)')

    # try to cluster the unclusterd:
    uncl = [c for c in cl.columns if len(cl.loc[0,c]) == 1]
    clust = list(set(cl.columns.tolist()) - set(uncl))
    tab = init_tab.copy()

    if len(uncl)>0: ## in case there are unclustered polys in the first step
        print('Clustering singletons...')
        for u in uncl:
            p1 = cl.loc[0,u][0]
            pair, max_row = tab.loc[p1].idxmax(), tab.loc[p1].max()
            for c in set(cl.columns.values)-set(uncl):
                if pair in cl.loc[0,c]:
                    cl.loc[0,c].append(p1)
                    cl.loc[1,c].append(max_row)
                    cl.drop(columns=[u],inplace=True)

    assert len([c for c in cl.columns if len(cl.loc[0,c]) == 1]) == 0

    clusters = pd.DataFrame()
    for c in cl.columns:
        c_polys = cl.loc[0,c]
        clusters = pd.concat([clusters, pd.DataFrame({'Cluster':[c for x in c_polys], "PolyID":c_polys, trajectory_group:[name for x in c_polys]})])

    return clusters


if __name__ == '__main__':

    #cluster = LocalCluster(n_workers=2, memory_limit='8 GiB')
    #client = Client(cluster)

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    console = Console()

    parser = argparse.ArgumentParser(description='Identifies shufflon signatures as covariated polys that are within a certain range. Works with a tmp/Poly_correlation.csv. Run pc_correlations.py before')
    parser.add_argument('-d', dest='workdir',help='Working directory.., Default is current working directory', default=os.getcwd())
    parser.add_argument('-m', dest='mode', help='invitro or invivo mode', choices=['invitro','invivo'], default='invivo')
    parser.add_argument('-p', dest='polytable', help='FinalPolyTable.tsv from bacevo pipeline.', required=True)
    args = parser.parse_args()

    os.chdir(args.workdir)

    df = pd.read_csv(args.polytable, sep='\t', low_memory=False, dtype={'Cage':'object'})
    if args.mode == 'invivo':
        trajectory_group = 'Mouse'
        df = df.loc[(df.SampleType == 'Feces')]

    else:
        trajectory_group = "Mix"

    """
    First define clusters of highly correlated polymorphisms.
    """
    clusters = pd.DataFrame()
    cor = pd.read_csv('tmp/Poly_correlation.csv')
    cor = cor.loc[cor.Mouse.str.contains('b')] ### filter for bern polys

    cutoff = 0.8
    cor.cor = cor.cor.apply(lambda x: abs(x)) ## transform to absolute correlation values
    cor = cor.loc[cor.cor.ge(cutoff)]
    # overwrite = True
    # if overwrite:
    #     clusters_vec = []
    #     for name,group in cor.groupby(trajectory_group):
    #         clusters.append(compute_poly_clusters2(name,group, trajectory_group))
    #     clusters = pd.concat(clusters)
    #     clusters.to_csv("tmp/clusters.csv", index=False)
    for name, group in cor.groupby(trajectory_group):
        clusters = pd.concat([clusters, compute_poly_clusters(name, group, trajectory_group)])
    #clusters = pd.read_csv('tmp/clusters.csv')

    """Column name is determined from the number of steps needed for clustering..Last step should be considered
    IMPORTANT: downstream in the df object is renamed to Cluster2
    """
    #clusterColName = np.max([int(x[-1]) for x in clusters.columns if x.startswith('Cluster')])

    df['Cluster'] = ['no cluster' for x in df.index] ## update FinalPolyTable with Cluster info

    for i, row in clusters.iterrows():
        df.loc[(df[trajectory_group] == row[trajectory_group]) & (df.PolyID == row.PolyID), 'Cluster'] = row.Cluster



    df.to_csv('tmp/Polys_in_clusters.csv',index=False) ## this is a FinalPolyTable with additional Cluster info
    clusters.to_csv("tmp/my_clusters.csv",index=False)
    go = pd.read_csv('gene_order/target_genes.csv')

    #df= pd.read_csv('tmp/my_clusters.csv')
    df = df.drop_duplicates(['SampleID','PolyID'])
    df['ClusterSize'] = [0 for x in df.index]
    df['ClusterSpread'] = [0 for x  in df.index]
    df['ClusterAveragePosition'] = [0 for x in df.index]
    df['Signature'] = [False for x in df.index]
    df['ClusterWithSignature'] = [False for x in df.index]

    """
    Now, for each cluster find the distance between polys in the cluster and group them as signature polys if found within a certain distance
    """
    def estimate_distance_cutoff(df):
        f = pd.DataFrame()
        cutoffs = [i for i in range(1,100,5)]
        for c in cutoffs:
            polys = []
            clusters = []
            genes = []
            for name,group in df.loc[df.Cluster != 'no cluster'].drop_duplicates(['Cluster','PolyID']).groupby("Cluster"):
                for x,y in itertools.combinations(group.PolyID.values, 2):
                    distance = abs(group.loc[group.PolyID == x, "Pos"].astype(int).values[0] - group.loc[group.PolyID == y, "Pos"].astype(int).values[0])
                    if distance <= c:
                        polys += [x,y]
                        clusters += [name]
                        genes += [group.loc[group.PolyID == x, "GeneName"].values[0], group.loc[group.PolyID == y, "GeneName"].values[0]]
            f = pd.concat([f, pd.DataFrame({'Distance':[c], "SignaturePolys":[len(set(polys))],"SignatureClusters":[len(set(clusters))], "SignatureGenes":[len(set(genes))]})])
        return f
    #e = estimate_distance_cutoff(df)
    #e.to_csv('tmp/estimate_distance_cutoff_for_signatures.csv')
    #g = sns.lineplot(data=e, x='Distance',y='SignatureGenes')
    #g.set_xticks(range(e.shape[0], 2))
    #g.set_xticklabels([(str(i) for i in range(e.shape[0],2))])
    #plt.show()
    distance_cutoff = 1000#50
    df['SignatureID'] = ['no' for x in df.index]
    sig_idx = 0
    for i, x  in zip(track(range(len(df.groupby("Cluster").groups))),df.groupby("Cluster")):
        name, group = x
        group = group.drop_duplicates("PolyID")
        if name != 'no cluster':
            df.loc[df['Cluster'] == name, 'ClusterSize'] = group.shape[0]
            min_pos = group.Pos.astype(int).min()
            max_pos = group.Pos.astype(int).max()
            #ave_pos = group.Pos.astype(int).mean()
            df.loc[df["Cluster"] == name, 'ClusterSpread'] = max_pos - min_pos
            #df.loc[df.Cluster == name, 'ClusterAveragePosition'] = ave_pos
            for x,y in itertools.combinations(group.PolyID.values, 2):
                if x != y: # this should be a redundant line

                    distance = abs(group.loc[group.PolyID == x, "Pos"].astype(int).values[0] - group.loc[group.PolyID == y, "Pos"].astype(int).values[0])
                    if distance <= distance_cutoff:
                        sig_idx += 1
                        #df.loc[(df.Cluster == name) & (df.PolyID.isin([x,y])), ["Signature",'SignatureID']] = [True, f'S{sig_idx}']

                        df.loc[(df.Cluster == name) & (df.PolyID.isin([x, y])), ["Signature", 'SignatureID']] = [True,
                                                                                                                 f'S{sig_idx}']
    for name,group in df.groupby("Cluster"):
        if group.Signature.any():
            df.loc[df["Cluster"] == name, "ClusterWithSignature"] = True



    df.sort_values(['Cluster','GeneName']).to_csv('tmp/my_clusters_stats.csv', index=False) ## this now exports a drop_duplicates(['SampleID','PolyID']) version of FinalPolyTable
    print(f'{df.loc[df.Cluster != "no cluster"].drop_duplicates("Cluster").shape[0]} clusters (CL) identified')
    #print(f'{df.loc[df.Cluster2 != "no cluster"].drop_duplicates("Cluster2").shape[0]} broad clusters (BCL) identified')

    df = df.loc[(df.Cluster.str.startswith("CL")) & (df.Signature == True)]
    print(f"Shufflon related clusters (CL): {df.drop_duplicates('Cluster').shape[0]}")
    #print(f"Shufflon related broad clusters (BCL): {df.drop_duplicates('Cluster2').shape[0]}")

    print(f"Total number of polymorphisms with shufflon signature:{df.drop_duplicates('PolyID').shape[0]}")
    print(f"Total number of genes with shufflon signatures:{df.drop_duplicates('GeneName').shape[0]}")
    print("\nList of genes with shufflon signatures:")
    print(df.drop_duplicates("GeneName")[['GeneName','oldGeneName', 'Signature','Cluster']].sort_values('GeneName'))


    print('\n\nDone!')







