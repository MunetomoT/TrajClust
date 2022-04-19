import scanpy as sc
import random
import pandas as pd
import numpy as np
import anndata as ad
from sklearn import metrics
import matplotlib as mpl
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram, single, complete, average, ward, set_link_color_palette
from tqdm import tqdm
from skmisc.loess import loess
from sktime.distances.elastic_cython import dtw_distance
import matplotlib.pyplot as plt

def TrajClust(testData, cloneID, pseudotime, supervised=False, cloneLabel=None, showDendogram=False, maxClusters=None, linkageCriterion="ward", pseudotimeSteps=0.5, impute=False, geneList=None, seed=0):
    """\
    Implementation of TrajClust. TrajClust expects a global cell ordering (pseudotime)
    to be assigned for each cell. Returns the cluster assignments of each clone and the clustering performance (Rand and NMI) if correct clusters are known.
    Parameters
    ----------
    testData
        Annotated data matrix (AnnData) with raw counts
    cloneID
        Identifier for each clone in AnnData observations
    pseudotime
        Identifier for pseudotime value of each cell in AnnData observations
    supervised
        True if clustering results known 
    cloneLabel 
        Identifier for clonal cluster identity of each cell in AnnData observations if cluster identities of clones known beforehand
    showDendogram
        Plot results of clustering as dendogram if True
    maxCluster
        Number of maximum clusters to find. If left undefined, determines optimum cluster number by maximising the silhouette score of the resulting clusters
    linkageCriterion
        Linkage criterion to use in final clustering
    pseudotimeSteps
        Step at which pseudotime values are incremented
    impute
        Impute genes of provided gene list if True
    geneList
        Genes to impute if impute is True
    seed
        Random seed for Louvain clustering
    -------
    Returns
        Dictionary containing cluster assignment of clones, calculated distance matrix and Rand and NMI (Normalised Mutualised Information) scores of the clusters if cluster identities of clones known beforehand
    """   
    
    results = {}

    # Create copy of data
    s = testData.copy()
    
    # Normalize, logarithmize and scale data
    sc.pp.normalize_per_cell(s, counts_per_cell_after=1e4)
    sc.pp.log1p(s)
    
    # Impute genes (optional)
    if impute:
        s = sc.external.pp.magic(s, random_state=seed, name_list=geneList)

    sc.pp.scale(s)
    
    # Get list of genes and clones
    genes = list(s.var.index)
    clones = list(s.obs[cloneID].unique())
    
    # Temporary data structures to hold loess results
    csv = []
    time = []
    
    # Iterate over each clone
    for clone in tqdm(clones):
        
        # Lay clones along pseudotime
        sample = s[s.obs[cloneID]==clone]
        x = sample.obs["pseudotime"]
        
        # Define times to sample
        times = np.arange(min(x),max(x),pseudotimeSteps)
        time.append(times)
        
        clone_values = []
        
        # For each gene
        for gene in genes:
            
            # Fit to gene trajectory using loess and compute values for sampling time
            try:
                y = sample[:,gene].X.flatten()
                l = loess(x,y)
                l.fit()
                pred = l.predict(times, stderror=True)
                lowess = list(pred.values)
                clone_values.append([lowess])
            except:
                print("Clone: "+ clone + " excluded")
                
        csv.append(clone_values)
    
    # Rearrange data 
    test = [0]*len(clones)
    for i,clone in enumerate(clones):
        value = np.array(csv[i]).reshape(len(genes), len(csv[i][0][0])).T
        test[i] = value
    
    feature = [s[s.obs[cloneID]==x].obs[cloneID][0] for x in clones]

    # Iterate over all gene trajectories and compute DTW distance
    # Save in distance matrix
    series_list = test
    n_series = len(series_list)
    distance_matrix = np.zeros(shape=(n_series, n_series))
    for i in tqdm(range(n_series)):
        for j in range(n_series):
            x = series_list[i]
            y = series_list[j]
            if i != j:
                dist = dtw_distance(x, y)
                distance_matrix[i, j] = dist

    # Compute linkage matrix based on distance matrix depending on linkage criterion
    if linkageCriterion == 'complete':
        linkage_matrix = complete(distance_matrix)
    if linkageCriterion == 'single':
        linkage_matrix = single(distance_matrix)
    if linkageCriterion == 'average':
        linkage_matrix = average(distance_matrix)
    if linkageCriterion == 'ward':
        linkage_matrix = ward(distance_matrix)
    
    # If cluster number not user-defined, find optimal cluster number
    if maxClusters is None:
        Z = linkage(distance_matrix, linkageCriterion, optimal_ordering=True)
        tempResults = []
        nmiResults = []
        for i in range(2, len(distance_matrix)):
            nodes = fcluster(Z, i, criterion='maxclust')
            tempResults.append(metrics.silhouette_score(distance_matrix, nodes, metric='precomputed'))
            nmiResults.append(metrics.normalized_mutual_info_score(nodes, feature))
        maxClusters = tempResults.index(max(tempResults))+2
        print("Optimal cluster number is " + str(maxClusters))
        
    # Get cluster labels
    cluster_labels = fcluster(linkage_matrix, maxClusters, criterion='maxclust')
    results["clusterLabels"] = cluster_labels
    results["distanceMatrix"] = distance_matrix
    feature = [s[s.obs[cloneID]==x].obs[cloneID][0] for x in clones]

    # If the clustering is supervised, calculate Rand and NMI
    if supervised:
        feature = [s[s.obs[cloneID]==x].obs[cloneLabel][0] for x in clones]
        rand = metrics.adjusted_rand_score(cluster_labels, feature)
        nmi = metrics.normalized_mutual_info_score(cluster_labels, feature)
        results["rand"] = rand
        results["nmi"] = nmi
        
    # Plot dendogram results if required
    if showDendogram:
        set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in mpl.cm.tab20.colors])
        color = linkage_matrix[:,2][-maxClusters+1]
        dendrogram(linkage_matrix, labels=feature, color_threshold=color)
        
    # Return dictionary containing cluster assignment of clones and calculated distance matrix (and Rand and NMI scores)
    return results

