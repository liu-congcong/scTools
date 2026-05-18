import warnings
from datetime import datetime
from uuid import uuid4

import numpy
import pandas
import scanpy
from pandas.errors import PerformanceWarning

from . import qc

warnings.filterwarnings("ignore", category = PerformanceWarning)

def readH5adFile(file):
    adata = scanpy.read_h5ad(file)
    return adata


def addGroup(adata, file):
    x = uuid4().hex
    y = dict()
    openFile = open(file, 'r')
    groups = openFile.readline().rstrip('\n').split('\t')[1 : ]
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        y['|'.join(lines[1 : ])] = lines[0]
    openFile.close()
    z = adata.obs[groups].astype(str).agg('|'.join, axis = 1)
    adata.obs.loc[ : , x] = z.replace(y).astype('category')
    return (x, z.isin(y))


def rankGenes(adata, group):
    groups = adata.obs[group]
    mean = numpy.empty(shape = (groups.cat.categories.size, adata.n_vars), dtype = numpy.float32)
    rank = numpy.empty(shape = (groups.cat.categories.size, adata.n_vars), dtype = numpy.float32)
    Rank = numpy.arange(adata.n_vars, 0, -1, dtype = numpy.float32) / (adata.n_vars / 100)
    for i, j in enumerate(groups.cat.categories):
        mean[i] = numpy.asarray(adata.X[(groups == j).values].mean(axis = 0)).reshape(-1)
        rank[i, numpy.argsort(mean[i])] = Rank
    mean = pandas.DataFrame(mean, index = groups.cat.categories, columns = adata.var_names, dtype = 'float32')
    rank = pandas.DataFrame(rank, index = groups.cat.categories, columns = adata.var_names, dtype = 'float32')
    return (mean, rank)


def filterGenes(adata, file):
    x = list()
    openFile = open(file, 'r')
    index = openFile.readline().rstrip('\n').split('\t').index('Gene')
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        x.append(lines[index])
    openFile.close()
    y = adata.var_names.str.replace(r'\.\d+', '', regex = True)
    adata = adata[ : , y.isin(x)].copy()
    return adata


def main(parameters):
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading: \"{parameters.input}\".', flush = True)
    adata = readH5adFile(parameters.input)
    group, mask = addGroup(adata, parameters.group)
    adata = adata[mask].copy()

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Performing quality control.', flush = True)
    adata = qc.main(adata, parameters.max_mt, parameters.min_genes, parameters.max_genes, parameters.min_cells, parameters.output + '.qc.pdf')

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Normalizing data.', flush = True)
    scanpy.pp.normalize_total(adata, target_sum = 1e4, layer = None, inplace = True)
    # adata = adata[ : , (adata.var['type'] == 'protein_coding') & (adata.var['pass'] >= 1)].copy()
    # adata = adata[ : , adata.var['pass'] >= 1].copy() #
    scanpy.pp.log1p(adata)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Ranking genes.', flush = True)
    mean, rank = rankGenes(adata, group)

    if parameters.marker is not None:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Filtering genes.', flush = True)
        adata = filterGenes(adata, parameters.marker)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Grouping cells.', flush = True)
    nGrouppedCells = adata.obs[group].value_counts()

    # openFile = open(parameters.output + 'groupCells.tsv', 'w')
    # openFile.write('Group\tCelles\n')
    # for i, j in nGrouppedCells.items():
    #     openFile.write(f'{i}\t{j}\n')
    # openFile.close()

    groups = nGrouppedCells[nGrouppedCells >= parameters.min_cells_group].index.tolist()
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Performing t-test for group comparison.', flush = True)
    scanpy.tl.rank_genes_groups(adata, groupby = group, groups = groups, reference = 'rest', method  = 't-test', corr_method = parameters.multiple_test, key_added = 'rank_genes_groups', n_genes = parameters.max_genes)
    # adata.uns['rank_genes_groups']['pvals_adj'] = numpy.minimum(adata.uns['rank_genes_groups']['pvals_adj'] * N, 1.0) #
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Writing: \"{parameters.output}\".', flush = True)
    openFile = open(parameters.output, 'w')
    openFile.write('Group\tGene ID\tGene Name\tMean (log1p)\tRank (%)\tLog2FC\tP\tQ\n')
    for i in groups:
        for gene, log2fc, p, q in zip(
            adata.uns['rank_genes_groups']['names'][i],
            adata.uns['rank_genes_groups']['logfoldchanges'][i],
            adata.uns['rank_genes_groups']['pvals'][i],
            adata.uns['rank_genes_groups']['pvals_adj'][i]
        ):
            if q <= parameters.alpha:
                geneName = adata.var['name'][gene]
                openFile.write(f'{i}\t{gene}\t{geneName}\t{mean.loc[i, gene]}\t{rank.loc[i, gene]}\t{log2fc}\t{p}\t{q}\n')
    openFile.close()
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Done.', flush = True)
    return None
