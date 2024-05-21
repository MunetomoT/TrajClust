<img width="72" alt="image" src="https://github.com/MunetomoT/TrajClust/assets/30675477/ac7bb56f-fd6b-420d-a674-5e434052c84b"># TrajClust

* Trajectory clustering of clone differentiation - TrajClust

Any questions and comments please contact
   
   mt739@cam.ac.uk
   
   jedt2@cam.ac.uk

Name
----

  TrajClust - Trajectory clustering of clone differentiation 

Synopsis
--------

  TrajClust(AnnData, cloneLabel, pseudotime)
   - TrajClust algorithm
  
  StdApproach(AnnData, cloneLabel, pseudotime)
   - Demonstration of performance achieveable by current standard approach to cluster clone differentiations
   
Description
-----------

The TrajClust algorithm aims to cluster groups of similar clone differentiations from single cell RNA datasets using information about each cell’s clone membership (matching TCR alpha and beta chains in the case of T cells), and a global ordering of each clone's cells.

The user will pass in a Annotated data matrix (AnnData) containing raw gene counts, cell ordering (pseudotime) and clone membership features of individual cells. This data is analyzed to find similar clone differentiations.

TrajClust returns a dictionary containing the cluster assignment of clones, with the number of clusters either user-defined, or determined by maximising the silhouette score of the resulting clusters. If the cluster identities of the clones are known beforehand, TrajClust returns the Rand and NMI (Normalised Mutualised Information) scores of the clusters found by the algorithm.

Schematic for toy dataset generation
<img width="824" alt="image" src="https://github.com/MunetomoT/TrajClust/assets/30675477/9427a9a6-c4ea-42e3-b6a4-41e334f330f4">

Schematic for StdApproach
<img width="832" alt="image" src="https://github.com/MunetomoT/TrajClust/assets/30675477/ff253d47-5558-4aef-98ba-e4c9607dcbe1">

Schematic for TrajClust
<img width="828" alt="image" src="https://github.com/MunetomoT/TrajClust/assets/30675477/cedfa3dd-21a5-4208-a01b-3ea40cf8100e">


