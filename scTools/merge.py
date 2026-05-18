import os
from datetime import datetime

import anndata
import pandas
import scanpy
from scipy.sparse import isspmatrix_csr


def readH5AD(file, obsKey):
    adata = scanpy.read_h5ad(file)
    for i in ['uns', 'obsm', 'obsp', 'varm', 'varp', 'layers']:
        if hasattr(adata, i):
            getattr(adata, i).clear()
    assert getattr(adata, 'raw', None) is not None, 'AnnData.raw is required but missing.'
    assert (getattr(adata.raw, 'X', None) is not None) and isspmatrix_csr(adata.raw.X), 'Only CSR sparse matrix is supported.'
    adata.X = adata.raw.X.copy()
    adata.X.eliminate_zeros()
    adata.raw = None
    adata.obs = adata.obs[obsKey]
    # adata.obs.rename(columns = dict(zip(obsKey, obsValue)), inplace = True) #

    adata.var = adata.var[['feature_name', 'feature_type', 'feature_length', 'feature_is_filtered']]
    adata.var.rename(columns = {'feature_name': 'name', 'feature_length': 'length', 'feature_type': 'type', 'feature_is_filtered': 'pass'}, inplace = True)
    return adata


def createGeneInfo(gene2info, adata):
    for gene, geneName, geneType, geneLength, geneQC in zip(adata.var_names, adata.var['name'], adata.var['type'], adata.var['length'], adata.var['pass']):
        flag = 0 if geneQC else 1
        if gene in gene2info:
            gene2info[gene][3] += flag
        else:
            gene2info[gene] = [geneName, geneLength, geneType, flag]
    return None


def addGeneInfo(adata, gene2info):
    geneInfo = pandas.DataFrame.from_dict(gene2info, orient = 'index', columns = ['name', 'length', 'type', 'pass'])
    geneInfo = geneInfo.astype({'name': 'category', 'length': 'int64', 'type': 'category', 'pass': 'int32'})
    adata.var = geneInfo
    return None


def main(parameters):
    adatas = list()
    groups = list()
    gene2info = dict()
    for file in parameters.input:
        assert file.lower().endswith('.h5ad'), 'Unsupported file format. Expected a .h5ad file.'
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading: \"{file}\".', flush = True)
        adata = readH5AD(file, parameters.obs_key)
        adatas.append(adata)
        createGeneInfo(gene2info, adata)
        groups.append(os.path.basename(file)[ : -5])
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Merging files.', flush = True)
    adata = anndata.concat(adatas, axis = 0, join = 'outer', label = 'organ', keys = groups, fill_value = 0, merge = None, uns_merge = None, pairwise = False)
    del adatas
    adata.X.eliminate_zeros()
    addGeneInfo(adata, gene2info)
    adata.obs_names_make_unique()
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Writing: \"{parameters.output}\".', flush = True)
    adata.write_h5ad(parameters.output, compression = None, compression_opts = None)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Done.', flush = True)
    return None
