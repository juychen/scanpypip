import scanpy.api as sc
import numpy as np
from matplotlib import pyplot as plot

def concat(adata_dict,join='inner', batch_key='sample', batch_categories=None, index_unique='-', combat = False, combat_key = "batch",covariates=None):
    '''
    This is a function to concat a dictonary of AnnData objects. Please predefine your own batch keys before the concatination

    Params:
    -------
    adata_dict: `dict`, {"identity":AnnData from scanpy}
        A dictionary: 
        Values should AnnData objects: https://anndata.readthedocs.io/en/latest/api.html#module-anndata
        Keys should be the identitiy of the adata object.
    
    join: `str`, optional (default: `"inner"`)
        Use intersection (``'inner'``) or union (``'outer'``) of variables.
    
    batch_key: `str`, (default: `"sample"`)
        Add the batch annotation to :attr:`obs` using this key.
    
    batch_categories: `str`,optional (default: `None`)
        Use these as categories for the batch annotation. By default, use increasing numbers.
    
    index_unique, `str`, optional (default: `"-"`)
        Make the index unique by joining the existing index names with the
        batch category, using ``index_unique='-'``, for instance. Provide
        ``None`` to keep existing indices.
    
    combat: `bool`, optional (defalut: `False`)
        Decide if to to the batch effect correction
    
    combat_key: `str`, optional (default: `"batch"`)
        Key to a categorical annotation from adata.obs that will be used for batch effect removal

    covariates
        Additional covariates such as adjustment variables or biological condition. Note that
        not including covariates may introduce bias or lead to the removal of biological signal 
        in unbalanced designs.

    inplace: bool, optional (default: `True`)
        Wether to replace adata.X or to return the corrected data

    
    Returns
    -------
    adata: AnnData
        The concatenated AnnData, where "adata.obs[batch_key]"
        stores a categorical variable labeling the original file identity.

    '''
    ada_keys = list(adata_dict.keys())

    adata = adata_dict[ada_keys[0]].concatenate([adata_dict[k] for k in ada_keys[1:]]
                                      ,join=join
                                      ,batch_key=batch_key
                                      ,batch_categories=batch_categories
                                      ,index_unique=index_unique)
    
    if(combat == False):
        return adata
    
    else:
        sc.pp.combat(adata, key=combat_key, covariates=covariates)

    return adata

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
