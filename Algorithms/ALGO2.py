import scanpy as sc
import random
import pandas as pd
import numpy as np
import anndata as ad
from sklearn import metrics
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram, single, complete, average, ward
from tqdm import tqdm
from skmisc.loess import loess
from sktime.distances.elastic_cython import dtw_distance

def ALGO2(testData, cloneID, pseudotime, supervised=False, cloneLabel=None, showDendogram=False, maxClusters=5, linkageCriterion="ward", pseudotimeSteps=0.5, seed=0):
    """\
    Implementation of Algorithm 2. Algorithm expects global pseudotime
    to be computed. Returns cluster_labels and metrics if labels known.
    Parameters
    ----------
    testData
        Annotated data matrix.
    cloneID
        Identifier for each clone
    pseudotime
        Identifier for pseudotime value of each cell
    supervised
        True if clustering results known 
    cloneLabel 
        Identifier for ClonalType if known (supervised)
    showDendogram
        Plot results of clustering as dendogram if True
    maxCluster
        Number of maximum clusters to find
    linkageCriterion
        Linkage criterion to use in final clustering
    pseudotimeSteps
        Step at which pseudotime values are incremented
    seed
        Random seed for Louvain clustering
    -------
    Returns
        Dictionary containing cluster assignments and metrics (if supervised)
    """   
    
    results = {}

    # Create copy of data
    s = testData.copy()
    
    # Get list of genes and clones
    genes = list(s.var.index)
    clones = list(s.obs["CloneID"].unique())
    
    # Temporary data structures to hold loess results
    csv = []
    time = []
    
    # Iterate over each clone
    for clone in tqdm(clones):
        
        # Lay clones along pseudotime
        sample = s[s.obs["CloneID"]==clone]
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
    
    # Get cluster labels
    cluster_labels = fcluster(linkage_matrix, maxClusters, criterion='maxclust')
    results["clusterLabels"] = cluster_labels
    
    # If the clustering is supervised, calculate rand and nmi
    if supervised:
        feature = [s[s.obs[cloneID]==x].obs[cloneLabel][0] for x in clones]
        rand = metrics.adjusted_rand_score(cluster_labels, feature)
        nmi = metrics.normalized_mutual_info_score(cluster_labels, feature)
        results["rand"] = rand
        results["nmi"] = nmi
        
    # Plot dendogram results if required
    if showDendogram:
        dendrogram(linkage_matrix, labels=feature)

    # Return dictionary containing cluster assignments and metrics (if supervised)
    return results

