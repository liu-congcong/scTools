import scanpy
from matplotlib import pyplot


class PlotQC:

    def __init__(self, width = 15, height = 6):
        self.gene2column = {'mt': 2}
        self.figure = pyplot.figure(figsize = (width, height), layout = 'constrained')
        self.axes = self.figure.subplots(nrows = 2, ncols = len(self.gene2column) + 2, squeeze = False)
        return None

    def basePlot(self, adata, flag):
        qcFlag = 'After' if flag else 'Before'
        self.axes[flag, 0].hist(adata.obs['total_counts'], bins = 50, edgecolor = 'black')
        self.axes[flag, 0].set_xlabel('Total Counts per Cell', fontsize = 10)
        self.axes[flag, 0].set_ylabel('Number of Cells', fontsize = 10)
        self.axes[flag, 0].set_title(f'Total UMI Counts Distribution ({qcFlag} QC)', fontsize = 12)
        self.axes[flag, 0].tick_params(axis = 'both', which = 'major', labelsize = 8)

        self.axes[flag, 1].hist(adata.obs['n_genes_by_counts'], bins = 50, edgecolor = 'black')
        self.axes[flag, 1].set_xlabel('Number of Genes per Cell', fontsize = 10)
        self.axes[flag, 1].set_ylabel('Number of Cells', fontsize = 10)
        self.axes[flag, 1].set_title(f'Detected Genes per Cell ({qcFlag} QC)', fontsize = 12)
        self.axes[flag, 1].tick_params(axis = 'both', which = 'major', labelsize = 8)

        for gene, column in self.gene2column.items():
            self.axes[flag, column].hist(adata.obs[f'pct_counts_{gene}'], bins = 50, edgecolor = 'black')
            self.axes[flag, column].set_xlabel(f'{gene.upper()} Gene Percentage', fontsize = 10)
            self.axes[flag, column].set_ylabel('Cell Count', fontsize = 10)
            self.axes[flag, column].set_title(f'Distribution of {gene.upper()} Gene Percentage ({qcFlag} QC)', fontsize = 12)
            self.axes[flag, column].tick_params(axis = 'both', which = 'major', labelsize = 8)
        return self

    def plotQC(self, adata, flag):
        self.basePlot(adata, flag)
        return self

    def saveFigure(self, file):
        pyplot.savefig(file)
        pyplot.close(self.figure)
        return self


def main(adata, maxMTPercentage, minNGenes, maxNGenes, minNCells, file):
    genes = adata.var['name'].str
    adata.var.loc[ : , 'mt'] = genes.startswith('MT-')
    adata.var.loc[ : , 'ribo'] = genes.startswith(('RPS', 'RPL'))
    scanpy.pp.calculate_qc_metrics(
        adata,
        expr_type = 'counts', var_type = 'genes', qc_vars = ['mt', 'ribo'],
        percent_top = None, log1p = False, inplace = True
    )

    plotQC = PlotQC()
    plotQC.plotQC(adata, 0)

    rowFlag = adata.obs.pct_counts_mt <= maxMTPercentage
    rowFlag &= adata.obs.n_genes_by_counts >= minNGenes
    rowFlag &= adata.obs.n_genes_by_counts <= maxNGenes

    adata = adata[rowFlag].copy()
    scanpy.pp.filter_genes(adata, min_cells = minNCells, inplace = True)
    adata = adata[ : , adata.X.max(axis = 0).toarray().flatten() > 0].copy()

    plotQC.plotQC(adata, 1)
    plotQC.saveFigure(file)
    return adata
