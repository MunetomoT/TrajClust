import scanpy as sc
import random
import pandas as pd
import numpy as np
import anndata as ad
from sklearn import metrics
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram


def ALGO1(testData, cloneID, supervised=False, cloneLabel=None, showDendogram=False, maxClusters=5, linkageCriterion="ward", resolution=50, seed=0):
    """\
    Naive implementation of Algorithm 1 with resolution parameter
    set at default 50. Returns clusterLabels and metrics if labels known.
    Parameters
    ----------
    testData
        Annotated data matrix.
    cloneID
        Identifier for each clone
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
    resolution
        Resolution paramater for clustering
    seed
        Random seed for Louvain clustering
    -------
    Returns
        Dictionary containing cluster assignments and metrics (if supervised)
    """    
    sc.settings.verbosity = 0
    results = {}
    
    # Take a copy of the data
    b = testData.copy()
    
    # Cluster the data
    sc.tl.louvain(b, resolution=resolution, random_state=seed)
    
    # Create matrix of cluster distribution for each clone
    matrix = pd.DataFrame(pd.DataFrame(b.obs[["louvain", cloneID]]).groupby(["louvain", cloneID]).size()/pd.DataFrame(b.obs[cloneID]).groupby([cloneID]).size()).unstack()
    
    # Use matrix to cluster results using hierarchical clustering
    linkage_matrix = linkage(matrix.T, linkageCriterion)
    cluster_labels = fcluster(linkage_matrix, maxClusters, criterion='maxclust')
    results["clusterLabels"] = cluster_labels
    feature = [b[b.obs[cloneID]==x].obs[cloneID][0] for x in matrix[0].columns] 
    
    # If the clustering is supervised, calculate rand and nmi
    if supervised:
        feature = [b[b.obs[cloneID]==x].obs[cloneLabel][0] for x in matrix[0].columns] 
        rand = metrics.adjusted_rand_score(cluster_labels, feature)
        nmi = metrics.normalized_mutual_info_score(cluster_labels, feature)
        results["rand"] = rand
        results["nmi"] = nmi
        
    # Plot dendogram results if required
    if showDendogram:
        dendrogram(linkage_matrix, labels=feature)
    
    # Return dictionary containing cluster assignments and metrics (if supervised)
    return results