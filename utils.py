import pandas as pd
import scanpy as sc
import numpy as np
from scipy import stats
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import fdrcorrection

def get_de_dataframe(adata,index):
        df_result = pd.DataFrame({
            "names":[item[index] for item in adata.uns["rank_genes_groups"]["names"]]
            ,"score":[item[index] for item in adata.uns["rank_genes_groups"]["scores"]]
            ,"logfoldchanges":[item[index] for item in adata.uns["rank_genes_groups"]["logfoldchanges"]]
            , "pvals":[item[index] for item in adata.uns["rank_genes_groups"]["pvals"]]
            , "pvals_adj":[item[index] for item in adata.uns["rank_genes_groups"]["pvals_adj"]]
             })
        
        return df_result

def self_annotate(adata,dict_marker,scale="Transpose",scaler=StandardScaler(),clustering_key="leiden",annotation_key="leiden_auto"):
    """Auto annoatate the cell types by ther marker genes according to the perasonr beatween the average\
        expression of each marker gene in each cell type and the marker gene identitiy matrix.

    Parameters
    ----------
    adata : AnnData
        The SCANPY AnnData object.
    
    dict_marker : dict
        A dictionary of marker genes. The key is the cell type name, the value is a list of marker genes.
    
    scaler : sklearn scaler
        A sklearn scaler object. Should have fit and transform methods.
        
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
        mk_louvain = list(set(adata.obs[clustering_key]))
    except Exception as e:
        print(e)
        print("No leiden information found in the adata object, please run leiden first.")
        return adata,None      
    mat_meanexp = np.zeros(shape=(len(mk_genes),len(mk_louvain)))
    df_mean_exp = pd.DataFrame(mat_meanexp,columns=mk_louvain,index=mk_genes)
    for l in mk_louvain:
        temp=adata.raw[adata.obs[clustering_key]==l]
        temp = temp[:,mk_genes]
        mean =temp.X.mean(axis=0)
        df_mean_exp[l] = np.array(mean).ravel()

    # Scale the mean expression matrix
    # scaler = StandardScaler()
    if(scale=="Transpose"):
        scaled_data=scaler.fit_transform(df_mean_exp.values.T)
        df_mean_exp = pd.DataFrame(scaled_data.T,columns=mk_louvain,index=mk_genes)
    else:
        scaled_data=scaler.fit_transform(df_mean_exp.values)
        df_mean_exp = pd.DataFrame(scaled_data,columns=mk_louvain,index=mk_genes)
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
    auto_annotate = [ct_map[t] for t in adata.obs[clustering_key]]
    adata.obs[annotation_key] = auto_annotate

    return adata,df_result

def cal_enrich_score(df,eps=np.finfo(float).eps,celltype_label='leiden_auto',virus_label='ebv'):
    """Calculate the enrichment score give a dataframe of the celltype and the virus.

    Parameters
    ----------
    df : DataFrame
        A dataframe having the cell type and virus labels.
    eps : float
        The epsilon value (a small float) for the enrichment score calculation.        
    celltype_label : str
        The column name of cell type label in the dataframe.
    virus_label : str
        The column name of virus label in the dataframe.

    Raises
    ------
    Error
        TBD

    Returns
    -------        
    df_result : DataFrame
        The dataframe result. Rows are the cell type, columns are the virus enrichment score.
    """

    df_obs = df
    if(df_obs[virus_label].dtype == "float32"):
        df_obs[virus_label] =  (df_obs[virus_label]>0).astype('string').astype('category')
    else:
        df_obs[virus_label] =  (df_obs[virus_label]).astype('string').astype('category')
        
    df_obs_group_ct_ebv = df_obs.groupby([celltype_label,virus_label])[virus_label].count().to_frame()
    df_obs_group_ebv = df_obs.groupby([virus_label])[virus_label].count().to_frame()
    df_obs_group_ct = df_obs.groupby([celltype_label])[celltype_label].count().to_frame()
    df_score = np.log(df_obs_group_ct_ebv + eps)/ (df_obs_group_ebv.loc["True"] + eps)
    df_result = df_score.iloc[df_score.index.get_level_values(virus_label) == "True"]
    df_result = df_result.reset_index(level=[0]).set_index(celltype_label)
    #df_result.columns[0] = "enrichment_score"
    return df_result

def cal_enrich_pval(adata,permutations=100,eps=np.finfo(float).eps,celltype_label='leiden_auto',virus_label='ebv',seed=38,alpha=0.05,\
                   method='indep',inplace=False):
    """Calculate the adjust pvale and enrichment score given an AnnData object with the celltype and the virus in the adata.obs.

    Parameters
    ----------
    adata : AnnData
        The SCANPY AnnData object.

    permutations : int
        The number times of permutations.
    
    eps : float
        The epsilon value (a small float) for the enrichment score calculation.        
    
    celltype_label : str
        The column name of cell type label in the dataframe.
    
    virus_label : str
        The column name of virus label in the dataframe.
    
    seed : int
        The seed for the random number generator.
    
    alpha : float
        The alpha value for the adjusted pvalue calculation. See the statsmodels.stats.multitest.fdrcorrection function for more details.
    
    method : str    
        The method for the adjusted pvalue calculation. See the statsmodels.stats.multitest.fdrcorrection function for more details.

    Raises
    ------
    Error
        TBD

    Returns
    -------        

    df_enrich : DataFrame  
        The dataframe result. Rows are the cell type, columns are the virus enrichment score with pval, adjpval.        
    """
    df_encirh = cal_enrich_score(adata.obs,eps=eps,celltype_label=celltype_label,virus_label=virus_label)
    df_permute = adata.obs.loc[:,[celltype_label,virus_label]]
    array_permute = np.zeros(shape=df_encirh.shape)
    list_permute = []
    for i in range(0,permutations):
        np.random.seed(i+seed)
        df_permute[virus_label] = np.random.permutation(df_permute[virus_label])
        df_permute_enrich = cal_enrich_score(df_permute,eps=eps,celltype_label=celltype_label,virus_label=virus_label)
        array_permute = array_permute + (df_encirh.values<=df_permute_enrich.values).astype('int')
        list_permute.append(df_permute_enrich.values)
    list_permute = np.array(list_permute).reshape(len(df_encirh),-1)
    pvals = array_permute/permutations
    rej,pval_adj = fdrcorrection((pvals).ravel(), alpha=alpha, method=method, is_sorted=False)
    df_encirh["pval_adj"+virus_label] = pval_adj
    df_encirh["pval"+virus_label] = pvals
    df_encirh.rename(columns={virus_label: "enrich_score"+virus_label},inplace=True)
    #adata.obs = adata.obs.merge(df_encirh,left_on=celltype_label,right_on=celltype_label)
    if(inplace==True):
        adata.obs["index"] = adata.obs.index
        df_obsmerge = adata.obs.merge(df_encirh,left_on=celltype_label,right_on=celltype_label)
        adata.obs = df_obsmerge.set_index("index").loc[adata.obs.index]
        return adata
        
    return df_encirh
