import sys
from argparse import ArgumentParser, RawTextHelpFormatter

from . import cluster, diffexp, info, marker, merge


class SCToolsParser(ArgumentParser):
    def error(self, message):
        self.print_help()
        print(f'\n\033[31mError: {message}\033[0m\n', flush = True)
        sys.exit()
        return None


def __init__():
    parser = SCToolsParser(
        formatter_class = RawTextHelpFormatter,
        description = 'A toolkit for single-cell sequencing data analysis in H5AD format.',
        epilog = 'Jing Guo, guojing@xjtu.edu.cn.\nCong-Cong Liu, congcong_liu@icloud.com.'
    )
    parser.add_argument(
        '-v', '--version', action = 'version', version = '%(prog)s 1.0.0'
    )

    subparsers = parser.add_subparsers(
        title = 'command', dest = 'command', required = True,
        description = 'Try "%(prog)s [command] -h|--help" for full help.',
    )
    mergeParser = subparsers.add_parser(
        'merge', formatter_class = RawTextHelpFormatter,
        help = 'Merge multiple H5AD files into a single file.',
        usage = '%(prog)s [options] -i H5ADs -o H5AD.'
    )
    mergeParser.add_argument(
        '-i', '--input', type = str, nargs = '+', required = True, metavar = '<str>',
        help = 'Path to the input h5ad files.\nThe basename will be treated as the annotation for all cells contained in each file.'
    )
    mergeParser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Path to the output h5ad file.'
    )
    mergeParser.add_argument(
        '--label', default = 'organ', type = str, required = False, metavar = '<str>',
        help = 'Default: organ.'
    )
    mergeParser.add_argument(
        '--obs-key', default = ['cell_type', ], type = str, nargs = '+', required = False, metavar = '<str>',
        help = 'Cell annotations to be retained in the merged file (e.g., donor_id, cell_type, anatomical_position).\nDefault: cell_type.'
    )

    diffexpParser = subparsers.add_parser(
        'diffexp', formatter_class = RawTextHelpFormatter,
        help='Perform differential analysis of genes between groups from an H5AD file.',
        usage = '%(prog)s [options] -i H5AD -g GROUP -o OUTPUT.'
    )
    diffexpParser.add_argument(
        '-i', '--input', type = str, required = True, metavar = '<str>',
        help = 'Path to the input h5ad file.'
    )
    diffexpParser.add_argument(
        '-g', '--group', type = str, required = True, metavar = '<str>',
        help = 'Path to the input group file.\nThe file must be tab-delimited with a header line.\nThe first column defines group names.\nAll other columns must match cell annotations in the H5AD file.\nUse "scTools info" to inspect available columns in the H5AD file.'
    )
    diffexpParser.add_argument(
        '-m', '--marker', type = str, required = False, metavar = '<str>',
        help = 'Path to the input marker file.\nSpecifies the subset of genes to be used for testing.\nThe file must be tab-delimited with a header line.\nA column \"Gene\" is required.'
    )
    diffexpParser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Path to the output file.'
    )
    diffexpParser.add_argument(
        '--max-mt', default = 10.0, type = float, required = False, metavar = '<float>',
        help = 'Maximum mitochondrial gene percentage allowed in a cell (%%).\nThe value must be a float between 0 and 100.\nDefault: 10.0.'
    )
    diffexpParser.add_argument(
        '--min-genes', default = 100, type = int, required = False, metavar = '<int>',
        help = 'Minimum number of detected genes in a cell.\nThe value must be a non-negative integer.\nDefault: 100.'
    )
    diffexpParser.add_argument(
        '--max-genes', default = 6000, type = int, required = False, metavar = '<int>',
        help = 'Maximum number of detected genes in a cell.\nThe value must be a non-negative integer.\nDefault: 6000.'
    )
    diffexpParser.add_argument(
        '--min-cells', default = 10, type = int, required = False, metavar = '<int>',
        help = 'Minimum number of cells expressing a gene.\nThe value must be a non-negative integer.\nDefault: 10.'
    )
    diffexpParser.add_argument(
        '--alpha', default = 0.05, type = float, required = False, metavar = '<float>',
        help = 'Significance threshold for statistical tests.\nThe value must be a float between 0 and 1.\nDefault: 0.05.'
    )
    diffexpParser.add_argument(
        '--multiple-test', default = 'benjamini-hochberg', type = str, choices = ('bonferroni', 'benjamini-hochberg'), required = False, metavar = '<str>',
        help = 'Method for multiple hypothesis testing.\nThe values can be bonferroni and benjamini-hochberg.\nDefault: benjamini-hochberg.'
    )
    diffexpParser.add_argument(
        '--min-cells-group', default = 2, type = int, required = False, metavar = '<int>',
        help = 'Minimum number of cells in a group.\nThe value must be a positive integer.\nDefault: 2.'
    )

    infoParser = subparsers.add_parser(
        'info', formatter_class = RawTextHelpFormatter,
        help = 'Extract information from an H5AD file.',
        usage = '%(prog)s [options] -i H5AD -o INFO.'
    )
    infoParser.add_argument(
        '-i', '--input', type = str, required = True, metavar = '<str>',
        help = 'Path to the input h5ad file.'
    )
    infoParser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Prefix for the output files.'
    )

    clusterParser = subparsers.add_parser(
        'cluster', formatter_class = RawTextHelpFormatter,
        help = 'Cluster groups.',
        usage = '%(prog)s [options] -i H5AD -g GROUP -o CLUSTER.'
    )
    clusterParser.add_argument(
        '-i', '--input', type = str, required = True, metavar = '<str>',
        help = 'Path to the input h5ad file.'
    )
    clusterParser.add_argument(
        '-g', '--group', type = str, required = True, metavar = '<str>',
        help = 'Path to the input group file.\nThe file must be tab-delimited with a header line.\nThe first column defines group names.\nAll other columns must match cell annotations in the H5AD file.\nUse "scTools info" to inspect available columns in the H5AD file.'
    )
    clusterParser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Path to the output file.'
    )
    clusterParser.add_argument(
        '--max-mt', default = 10.0, type = float, required = False, metavar = '<float>',
        help = 'Maximum mitochondrial gene percentage allowed in a cell (%%).\nThe value must be a float between 0 and 100.\nDefault: 10.0.'
    )
    clusterParser.add_argument(
        '--min-genes', default = 100, type = int, required = False, metavar = '<int>',
        help = 'Minimum number of detected genes in a cell.\nThe value must be a non-negative integer.\nDefault: 100.'
    )
    clusterParser.add_argument(
        '--max-genes', default = 6000, type = int, required = False, metavar = '<int>',
        help = 'Maximum number of detected genes in a cell.\nThe value must be a non-negative integer.\nDefault: 6000.'
    )
    clusterParser.add_argument(
        '--min-cells', default = 10, type = int, required = False, metavar = '<int>',
        help = 'Minimum number of cells expressing a gene.\nThe value must be a non-negative integer.\nDefault: 10.'
    )
    clusterParser.add_argument(
        '--tsne-perplexity', default = 1, type = int, required = False, metavar = '<int>',
        help = 'The perplexity parameter in t-SNE.\nThe value must be a positive integer.\nDefault: 1.'
    )
    clusterParser.add_argument(
        '--dbscan-eps', default = 5.0, type = float, required = False, metavar = '<float>',
        help = 'The eps parameter in DBSCAN.\nThe value must be a positive float.\nDefault: 5.'
    )
    clusterParser.add_argument(
        '--width', default = 7.0, type = float, required = False, metavar = '<float>',
        help = 'Figure width in inches.\nThe value must be a positive float.\nDefault: 7.0.'
    )
    clusterParser.add_argument(
        '--height', default = 7.0, type = float, required = False, metavar = '<float>',
        help = 'Figure height in inches.\nThe value must be a positive float.\nDefault: 7.0.'
    )

    markerParser = subparsers.add_parser(
        'marker', formatter_class = RawTextHelpFormatter,
        help = 'Identify group-specific marker genes.',
        usage = '%(prog)s [options] -t TEST -o MARKER.'
    )
    markerParser.add_argument(
        '-d', '--diffexp', type = str, required = True, metavar = '<str>',
        help = 'Path to the input differential analysis file.'
    )
    markerParser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Path to the output file.'
    )
    markerParser.add_argument(
        '-m', '--mapping', type = str, required = False, metavar = '<str>',
        help = 'Path to the input mapping file.\nThe file must be tab-delimited with a header line: \"Gene<tab>Protein\".'
    )
    markerParser.add_argument(
        '--log2fc', default = 2, type = float, required = False, metavar = '<float>',
        help = 'Minimum log2 fold-change of mean expression in a group compared to all other groups.\nDefault: 2.'
    )
    markerParser.add_argument(
        '--rank', default = 10, type = float, required = False, metavar = '<float>',
        help = 'Genes with expression rank greater than \"--rank\" are removed.\nThe value must be a float between 0 and 100.\nDefault: 10.'
    )
    markerParser.add_argument(
        '--mean', default = 0, type = float, required = False, metavar = '<float>',
        help = 'Genes with mean expression (log1p) lower than \"--mean\" are removed.\nDefault: 0.'
    )
    markerParser.add_argument(
        '--coverage', type = int, required = False, metavar = '<int>',
        help = 'Genes identified as markers in more than \"--coverage\" of groups will be removed.\nThe value must be a positive integer.\nDefault: inf.'
    )
    markerParser.add_argument(
        '--hits', type = int, required = False, metavar = '<int>',
        help = 'Genes are only considered marker genes in the top \"--hits\" groups.\nThe value must be a positive integer.\nDefault: inf.'
    )
    return parser.parse_args()


def main():
    parameters = __init__()
    if parameters.command == 'merge':
        merge.main(parameters)
    elif parameters.command == 'diffexp':
        diffexp.main(parameters)
    elif parameters.command == 'info':
        info.main(parameters)
    elif parameters.command == 'cluster':
        cluster.main(parameters)
    elif parameters.command == 'marker':
        marker.main(parameters)
    return None
