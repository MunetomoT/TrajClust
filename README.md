# TrajClust

* Trajectory clustering of clone differentiation - TrajClust

Any questions and comments please contact
   
   mt739@mrc-tox.cam.ac.uk
   
   jedt2@mrc-tox.cam.ac.uk

Name
----

  TrajClust - Trajectory clustering of clone differentiation 

Synopsis
--------

  TrajClust(AnnData, clone_label, clone_axis)
   - TrajClust algorithm
  
  StdApproach(AnnData, clone_label, clone_axis)
   - Demonstration of maximum performance achieveable by current standard approach to cluster clone differentiations
   
Description
-----------

The TrajClust algorithm aims to cluster groups of similar clone differentiations from single cell RNA datasets using information about each cell’s clone membership (matching TCR alpha and beta chains in the case of T cells), and a global ordering of each clone's cells.

The user will pass in a Annotated data matrix (AnnData) containing raw gene counts, cell ordering (pseudotime) and clone membership features of individual cells. This data is analyzed to find similar clone differentiations.

TrajClust returns a dictionary containing the cluster assignment of clones, with the number of clusters either user-defined, or determined by maximising the silhouette score of the resulting clusters. If the cluster identities of the clones are known beforehand, TrajClust returns the Rand and NMI (Normalised Mutualised Information) scores of the clusters found by the algorithm.
