import scanpy as sc
import random
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
from sklearn import metrics
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram
from tqdm import tqdm

def ALGO1STAR(testData, cloneID, cloneLabel, showCurves=False, showDendogram=False, maxClusters=5, linkageCriterion="ward", maxResolution=60, steps=1, seed=0):
    """\
    Supervised implementation to show theoretical maximum of Algorithm 1 
    with resolution parameter chosen for optimal performance.
    Returns clusterLabels, best resolution and metrics.
    Parameters
    ----------
    testData
        Annotated data matrix.
    cloneID
        Identifier for each clone
    cloneLabel 
        Identifier for ClonalType 
    showDendogram
        Plot results of clustering as dendogram if True
    showCurves
        Plot NMI and Rand curves if True
    maxCluster
        Number of maximum clusters to find
    linkageCriterion
        Linkage criterion to use in final clustering
    maxResolution
        Maximum resolution paramater to try for clustering
    steps
        Adjust increments that resolution is increased by
    seed
        Random seed for Louvain clustering
    -------
    Returns
        Dictionary containing cluster assignments and metrics (if supervised)
    """    
    sc.settings.verbosity = 0
    results = {}
    rand = []
    nmi = []
    
    # Take a copy of the data
    b = testData.copy()
    
    # Set resolution parameters to iterate over
    res = list(np.arange(0, maxResolution, steps))

    # Iterate over these resolutions
    for i,x in tqdm(enumerate(res)):

        # Cluster the data
        sc.tl.louvain(b, resolution=x, random_state=seed)

        # Create matrix of cluster distribution for each clone
        matrix = pd.DataFrame(pd.DataFrame(b.obs[["louvain", cloneID]]).groupby(["louvain", cloneID]).size()/pd.DataFrame(b.obs[cloneID]).groupby([cloneID]).size()).unstack()

        # Use matrix to cluster results using hierarchical clustering
        linkage_matrix = linkage(matrix.T, linkageCriterion)
        cluster_labels = fcluster(linkage_matrix, maxClusters, criterion='maxclust')
        feature = [b[b.obs[cloneID]==x].obs[cloneLabel][0] for x in matrix[0].columns] 
        
        # Save Rand and NMI for each resolution
        rand.append(metrics.adjusted_rand_score(cluster_labels, feature))
        nmi.append(metrics.normalized_mutual_info_score(cluster_labels, feature))
    
    
    # Find optimal resolution (based on rand)
    max_rand = max(rand)
    max_index = rand.index(max_rand)
    max_nmi = nmi[max_index]
    max_res = res[max_index]
    
    # Rerun with maximum resolution
    
    # Cluster the data
    sc.tl.louvain(b, resolution=max_res, random_state=seed)

    # Create matrix of cluster distribution for each clone
    matrix = pd.DataFrame(pd.DataFrame(b.obs[["louvain", cloneID]]).groupby(["louvain", cloneID]).size()/pd.DataFrame(b.obs[cloneID]).groupby([cloneID]).size()).unstack()

    # Use matrix to cluster results using hierarchical clustering
    linkage_matrix = linkage(matrix.T, linkageCriterion)
    cluster_labels = fcluster(linkage_matrix, maxClusters, criterion='maxclust')
    feature = [b[b.obs[cloneID]==x].obs[cloneLabel][0] for x in matrix[0].columns] 
    
    results["clusterLabels"] = cluster_labels
    results["rand"] = max_rand
    results["nmi"] = max_nmi

    # Plot NMI and Rand results if required
    if showCurves:
        plt.plot(res, rand)
        plt.plot(res, nmi)
        plt.show()
        
    # Plot dendogram results if required
    if showDendogram:
        dendrogram(linkage_matrix, labels=feature)
    
    # Return dictionary containing cluster assignments and metrics (if supervised)
    return results