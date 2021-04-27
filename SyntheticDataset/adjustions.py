import scanpy as sc
import random
import pandas as pd
import numpy as np
import anndata as ad

def adjust(adata, genes, CloneLabel, pU, pL, grad, size, seed):
    """\
    Adjust the gene trajectory of cells so that cells from the same clone
    have similar gene trajectory behaviours. Adjustions made in-place
    Parameters
    ----------
    adata
        Annotated data matrix.
    genes
        Genes to be adjusted
    CloneLabel 
        Identifier for ClonalType
    pU
        Proportion of basis functions that are u-basis
    pL
        Proportion of basis functions that are l-basis
    grad
        If adjustion is increase then True else False
    size
        Size of clone
    seed
        Random seed
    -------
    """    
        
    # Set random seed
    random.seed(seed)
    
    # Subset data to only include clones of interest and find clone names
    temp = adata[adata.obs["CloneLabel"]==CloneLabel]
    clones = list(temp.obs["CloneID"].unique())
    
    # Calculate proportion of basis functions used
    pE = 1-pU-pL
    funcs = ["e"]*round(pE*len(genes))+["u"]*round(pU*len(genes))+["l"]*round(pL*len(genes))
    print(funcs)
    
    # For each gene that needs to be adjusted
    for s,g in enumerate(genes):
        
        # Find an example clone and calculate parameters
        temp2 = temp[temp.obs["CloneID"]==clones[0]]
        index = list(temp2.obs["pseudotime"].sort_values().index)
        reset_index = list(temp2.obs["pseudotime"].reset_index().sort_values('pseudotime').index)
        changedX = np.array(temp2.X)[[int(x) for x in reset_index]][:,g]
        times = list(temp2.obs["pseudotime"].sort_values())
        low = random.uniform(min(changedX), max(changedX)/2)
        high = random.uniform(max(changedX)/2, max(changedX))
        dropout = sum(i <= 0 for i in changedX)/len(changedX)
        t0, t1 = min(times), max(times)
        dt = t1-t0
        dy = high-low
        pos = max(min(np.random.normal(0.5, 0.25), 0), 1)*dt+t0
        
        # For each clone
        for i, clone in enumerate(clones):
            
            # Subset data and lay out gene trajectory
            temp2 = temp[temp.obs["CloneID"]==clones[i]]
            index = list(temp2.obs["pseudotime"].sort_values().index)
            reset_index = list(temp2.obs["pseudotime"].reset_index().sort_values('pseudotime').index)
            changedX = np.array(temp2.X)[[int(x) for x in reset_index]][:,g]
            times = list(temp2.obs["pseudotime"].sort_values())
            
            # Adjust gene counts based on parameters
            newX = [max(f(y, t, dt, dy, t0, low, high, grad, pos, funcs, s),0) 
                    for t,y in zip(times, changedX)]
            newX = [0 if random.uniform(0, 1)<dropout else x for x in newX]
            
            # Assign new values in-place
            for i in range(0, size):
                adata.X[int(index[i]), g] = newX[i]
                
def f(y, t, dt, dy, t0, low, high, grad, pos, funcs, s):
    """\
    Computes basis functions using given parameters
    Parameters
    ----------
    y
        Original gene counts
    t
        Pseudotime values
    dt
        Range of pseudotime values
    dy
        Range of gene counts
    t0
        Minimum pseudotime value
    low
        Minimum gene count
    high
        Maximum gene count
    grad
        If adjustion is increase then True else False
    pos
        Location of inflex
    funcs
        List of basis function required
    s
        Index into funcs
    -------
    Returns
        Basis function output for given value
    """    
    # If e-basis function required
    if funcs[s]=="e":
        if grad==True:
            return int(high*(1-np.exp(-(t-t0)/(pos*dt)))*0.5+y*0.5)
        else:
            return int(high*(1-np.exp((t-t0)/(pos*dt)))*0.5+y*0.5)
    
    # Id u-basis function required
    elif funcs[s]=="u":
        a=(low-high)/(pos**2)
        if grad==True:
            return int((a*(t-pos)**2+high)*0.5+y*0.5)
        else:
            return int((-a*(t-pos)**2+low)*0.5+y*0.5)
    
    # Id l-basis function required
    else:
        if grad==True:
            return int(((dy/dt)*(t-t0)+low)*0.5+y*0.5)
        else:
            return int(((-dy/dt)*(t-t0)+low)*0.5+y*0.5)