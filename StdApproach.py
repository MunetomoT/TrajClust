import scanpy as sc
import random
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
from sklearn import metrics
import matplotlib as mpl
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram, set_link_color_palette
from tqdm import tqdm

def StdApproach(testData, cloneID, supervised=False, cloneLabel=None, showDendogram=False, maxClusters=None, linkageCriterion="ward", resolution=50, impute=False, geneList=None, seed=0):
    """\
    Demonstration of performance achieveable by a current standard approach. The resolution parameter can be modified to achieve higher performance. Returns the cluster assignments of each clone and the clustering performance (Rand and NMI) if correct clusters are known.
    Parameters
    ----------
    testData
        Annotated data matrix with raw counts.
    cloneID
        Identifier for each clone
    supervised
        True if clustering results known
    cloneLabel 
        Identifier for ClonalType 
    showDendogram
        Plot results of clustering as dendogram if True
     maxCluster
        Number of maximum clusters to find. If left undefined, determines optimum cluster number by maximising the silhouette score of the resulting clusters
    linkageCriterion
        Linkage criterion to use in final clustering
    resolution
        Resolution paramater for clustering
    impute
        Impute genes of provided gene list if True
    geneList
        Genes to impute if impute is True
    seed
        Random seed for Louvain clustering
    -------
    Returns
        Dictionary containing cluster assignment of clones, calculated feature matrix and Rand and NMI (Normalised Mutualised Information) scores of the clusters if cluster identities of clones known beforehand
    """
    sc.settings.verbosity = 0
    results = {}
    rand = []
    nmi = []
    
    # Take a copy of the data
    b = testData.copy()
    
    # Normalize, logarithmize and scale data
    sc.pp.normalize_per_cell(b, counts_per_cell_after=1e4)
    sc.pp.log1p(b)
    
    # Impute genes (optional)
    if impute:
        b = sc.external.pp.magic(b, random_state=seed, name_list=geneList)
        
    sc.pp.scale(b)
    
    # Cluster the data
    sc.pp.neighbors(b)
    sc.tl.louvain(b, resolution=resolution, random_state=seed)
    
    # Create matrix of cluster distribution for each clone
    matrix = pd.DataFrame(pd.DataFrame(b.obs[["louvain", cloneID]]).groupby(["louvain", cloneID]).size()/pd.DataFrame(b.obs[cloneID]).groupby([cloneID]).size()).unstack()
    
    # Use matrix to cluster results using hierarchical clustering
    linkage_matrix = linkage(matrix.T, linkageCriterion)

    # If cluster number not user-defined, find optimal cluster number
    if maxClusters is None:
        Z = linkage(matrix.T, linkageCriterion, optimal_ordering=True)
        tempResults = []
        for i in range(2, len(matrix.T)):
            nodes = fcluster(Z, i, criterion='maxclust')
            tempResults.append(metrics.silhouette_score(matrix.T, nodes))
        maxClusters = tempResults.index(max(tempResults))+2
        print("Optimal cluster number is " + str(maxClusters))
    
    cluster_labels = fcluster(linkage_matrix, maxClusters, criterion='maxclust')
    results["clusterLabels"] = cluster_labels
    results["featureMatrix"] = matrix.T
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
        set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in mpl.cm.tab20.colors])
        color = linkage_matrix[:,2][-maxClusters+1]
        dendrogram(linkage_matrix, labels=feature, color_threshold=color)

    # Return dictionary containing cluster assignment of clones and calculated feature matrix (and Rand and NMI scores)
    return results
