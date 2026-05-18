from datetime import datetime

import scanpy


def readH5adFile(file):
    adata = scanpy.read_h5ad(file, backed = 'r')
    return adata


def getRowInfo(adata):
    x = adata.obs.columns.tolist()
    y = adata.obs[x].astype(str).agg('\t'.join, axis = 1).astype('category')
    x.append('#cells')
    yield('\t'.join(x))
    for i, j in y.value_counts().items():
        yield f'{i}\t{j}'
    return None


def getColumnInfo(adata):
    x = adata.var.columns.tolist()
    y = adata.var[x].astype(str).agg('\t'.join, axis = 1).astype('category')
    yield('gene\t' + '\t'.join(x))
    for i, j in y.items():
        yield f'{i}\t{j}'
    return None


def main(parameters):
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading: \"{parameters.input}\".', flush = True)
    adata = readH5adFile(parameters.input)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Writing: \"{parameters.output}.row.tsv\".', flush = True)
    openFile = open(parameters.output + '.row.tsv', 'w')
    for i in getRowInfo(adata):
        openFile.write(i + '\n')
    openFile.close()
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Writing: \"{parameters.output}.column.tsv\".', flush = True)
    openFile = open(parameters.output + '.column.tsv', 'w')
    for i in getColumnInfo(adata):
        openFile.write(i + '\n')
    openFile.close()
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Done.', flush = True)
    return None
