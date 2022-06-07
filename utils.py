import pandas as pd
import scanpy as sc
import torch
import numpy as np
from scipy import stats
from sklearn.preprocessing import StandardScaler

def get_de_dataframe(adata,index):
        df_result = pd.DataFrame({
            "names":[item[index] for item in adata.uns["rank_genes_groups"]["names"]]
            ,"score":[item[index] for item in adata.uns["rank_genes_groups"]["scores"]]
            ,"logfoldchanges":[item[index] for item in adata.uns["rank_genes_groups"]["logfoldchanges"]]
            , "pvals":[item[index] for item in adata.uns["rank_genes_groups"]["pvals"]]
            , "pvals_adj":[item[index] for item in adata.uns["rank_genes_groups"]["pvals_adj"]]
             })
        
        return df_result

def self_annotate(adata,dict_marker):
    """Auto annoatate the cell types by ther marker genes according to the perasonr beatween the average\
        expression of each marker gene in each cell type and the marker gene identitiy matrix.

    Parameters
    ----------
    adata : AnnData
        The SCANPY AnnData object.
    
    dict_marker : dict
        A dictionary of marker genes. The key is the cell type name, the value is a list of marker genes.
        
    Raises
    ------
    NotImplementedError
        If no sound is set for the animal or passed in as a
        parameter.

    Returns
    -------
    adata : AnnData
        The SCANPY AnnData object. With the cell type annotations, named as `leiden_auto`.

    df_result : DataFrame
        The dataframe of the pearsonr .
    """
    # get the marker gene identity matrix
    mk_types = dict_marker.keys()
    mk_genes = []
    for v in dict_marker.values():
        mk_genes = mk_genes+v
    mat = np.zeros(shape=(len(mk_genes),len(mk_types)))
    df_mat = pd.DataFrame(mat,columns=mk_types,index=mk_genes)
    for k in dict_marker:
        df_mat.loc[dict_marker[k],k]=1
    
    # get the average expression of each marker gene in each cell type
    try:
        mk_louvain = list(set(adata.obs.leiden))
    except Exception as e:
        print(e)
        print("No leiden information found in the adata object, please run leiden first.")
        return adata,None      
    mat_meanexp = np.zeros(shape=(len(mk_genes),len(mk_louvain)))
    df_mean_exp = pd.DataFrame(mat_meanexp,columns=mk_louvain,index=mk_genes)
    for l in mk_louvain:
        temp=adata.raw[adata.obs.leiden==l]
        temp = temp[:,mk_genes]
        mean =temp.X.mean(axis=0)
        df_mean_exp[l] = mean

    # Scale the mean expression matrix
    scaler = StandardScaler()
    scaled_data=scaler.fit_transform(df_mean_exp.values.T)
    df_mean_exp = pd.DataFrame(scaled_data.T,columns=mk_louvain,index=mk_genes)
    resutl_mat = np.zeros(shape=(len(mk_louvain),len(mk_types)))
    df_result = pd.DataFrame(resutl_mat,columns=mk_types,index=mk_louvain)

    # Calculate the pearsonr between the mean expression of each marker gene and the marker gene identity matrix
    for ct in mk_types:
        for l in mk_louvain:
            ctg = df_mat[ct].values
            lg = df_mean_exp[l].values
            df_result.loc[l,ct]=stats.pearsonr(ctg, lg)[0]

    # User identity matrix to annotate the cell type       
    result_idmax = (df_result.idxmax(axis=1)).to_frame()
    ct_map = {}
    for index, row in result_idmax.iterrows():
        ct_map[index]=row[0]

    # Store the cell type annotation in the leiden_auto column
    auto_annotate = [ct_map[t] for t in adata.obs['leiden']]
    adata.obs['leiden_auto'] = auto_annotate

    return adata,df_result