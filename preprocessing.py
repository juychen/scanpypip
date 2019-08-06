import scanpy.api as sc
import numpy as np
from matplotlib import pyplot as plot

def concat(datats):

def cal_ncount_ngenes(adata,sparse=False):
    
    mito_genes = adata.var_names.str.startswith('mt-')

    if sparse == False:    
        adata.obs['n_counts'] = adata.X.sum(axis=1)
        adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
    else:
        adata.obs['n_counts'] = adata.X.sum(axis=1).A1
        adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    return adata

def receipe_my(adata,l_n_genes = 500, r_n_genes= 5000, percent_mito = 0.05, normalize = False,log = False,sparse = False,plotinfo= False):

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    adata = cal_ncount_ngenes(adata)

    # if sparse == False:    
    #     adata.obs['n_counts'] = adata.X.sum(axis=1)
    #     adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
    # else:
    #     adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    #     adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1

    adata = adata[
        np.logical_and(
        (adata.obs['n_genes'] > l_n_genes), 
        (adata.obs['n_genes'] < r_n_genes)),:]
    adata = adata[adata.obs['percent_mito'] < percent_mito, :]

    if(plotinfo!=False):
        sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, save=True)
        #plt.savefig(plotinfo)



    print(adata.shape)
    
    if normalize == True:
        sc.pp.normalize_total(adata)
    #sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    adata.raw = adata

    if log == True:
        sc.pp.log1p(adata)

    #sc.pp.scale(adata, max_value=10)

    return adata
