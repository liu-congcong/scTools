from datetime import datetime
from uuid import uuid4

import numpy
import scanpy
from matplotlib.pyplot import close, subplots
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
from sklearn.metrics.pairwise import cosine_distances

from . import color, qc


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


def plot(xy, y, z, output, width, height):
    colors = color.generateColors(numpy.unique(z).size)
    randomGenerator = numpy.random.default_rng(0)
    figure, axes = subplots(nrows = 1, ncols = 1, figsize = (width, height), layout = 'constrained', squeeze = False)
    axes[0, 0].scatter(xy[ : , 0], xy[ : , 1], s = 20, alpha = 0.9, facecolors = 'none', edgecolors = "#4B4B4B", marker = 'o', rasterized = False)
    for i, j, k in zip(xy, y, z):
        axes[0, 0].annotate(
            j, xy = i, xytext = (i[0], i[1]),
            arrowprops = None,
            size = 9,
            color = colors[k],
            horizontalalignment = 'center',
            verticalalignment = 'center',
            rotation = randomGenerator.integers(-45, 45, 1)[0]
        )
    deltaX = numpy.max(xy[ : , 0]) - numpy.min(xy[ : , 0])
    deltaY = numpy.max(xy[ : , 1]) - numpy.min(xy[ : , 1])
    axes[0, 0].set_xlim(numpy.min(xy[ : , 0]) - 0.1 * deltaX, numpy.max(xy[ : , 0]) + 0.1 * deltaX)
    axes[0, 0].set_ylim(numpy.min(xy[ : , 1]) - 0.1 * deltaY, numpy.max(xy[ : , 1]) + 0.1 * deltaY)
    axes[0, 0].tick_params(labelsize = 8)
    figure.savefig(output, bbox_inches = 'tight')
    close()
    return None


def main(parameters):
    if parameters.output.lower().endswith('.pdf'):
        parameters.output = parameters.output[ : -4]

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading: \"{parameters.input}\".', flush = True)
    adata = readH5adFile(parameters.input)
    group, mask = addGroup(adata, parameters.group)
    adata = adata[mask].copy()

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Performing quality control.', flush = True)
    adata = qc.main(adata, parameters.max_mt, parameters.min_genes, parameters.max_genes, parameters.min_cells, parameters.output + '.qc.pdf')

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Normalizing data.', flush = True)
    scanpy.pp.normalize_total(adata, target_sum = 1e4, layer = None, inplace = True)
    # adata = adata[ : , adata.var['type'] == 'protein_coding'].copy() #
    scanpy.pp.log1p(adata)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Calculating the mean of genes across all groups.', flush = True)
    x = numpy.empty(shape = (adata.obs[group].cat.categories.size, adata.n_vars), dtype = numpy.float32)
    y = list()
    for i, j in enumerate(adata.obs[group].cat.categories):
        x[i] = adata.X[(adata.obs[group] == j).values].mean(axis = 0)
        y.append(j)
    del adata

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Calculating distances between groups.', flush = True)
    d = cosine_distances(x)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Computing t-SNE embeddings.', flush = True)
    tsne = TSNE(n_components = 2, random_state = 0, max_iter = 9999, metric = 'precomputed', init = 'random', perplexity = parameters.tsne_perplexity)
    tsneXY = tsne.fit_transform(d)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Performing DBSCAN clustering on embeddings.', flush = True)
    dbscan = DBSCAN(eps = parameters.dbscan_eps, min_samples = 1)
    dbscanY = dbscan.fit_predict(tsneXY)
    plot(tsneXY, y, dbscanY, parameters.output + '.pdf', parameters.width, parameters.height)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Done.', flush = True)

    return None
