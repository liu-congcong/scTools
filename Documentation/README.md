# SCTools

A toolkit for single-cell sequencing data analysis in H5AD format.

## Commands

- [`merge`](./README.md#merge)     Merge multiple H5AD files into a single file.
- [`diffexp`](./README.md#diffexp) Perform differential analysis of genes between groups from an H5AD file.
- [`info`](./README.md#info)       Extract information from an H5AD file.
- [`cluster`](./README.md#cluster) Cluster groups.
- [`marker`](./README.md#marker)   Identify group-specific marker genes.

### merge

```text
usage: scTools merge [options] -i H5ADs -o H5AD.

options:
  -h, --help            show this help message and exit
  -i, --input <str> [<str> ...]
                        Path to the input h5ad files.
                        The basename will be treated as the annotation for all cells contained in each file.
  -o, --output <str>    Path to the output h5ad file.
  --label <str>         Default: organ.
  --obs-key <str> [<str> ...]
                        Cell annotations to be retained in the merged file (e.g., donor_id, cell_type, anatomical_position).
                        Default: cell_type.
```

### diffexp

```text
usage: scTools diffexp [options] -i H5AD -g GROUP -o OUTPUT.

options:
  -h, --help            show this help message and exit
  -i, --input <str>     Path to the input h5ad file.
  -g, --group <str>     Path to the input group file.
                        The file must be tab-delimited with a header line.
                        The first column defines group names.
                        All other columns must match cell annotations in the H5AD file.
                        Use "scTools info" to inspect available columns in the H5AD file.
  -m, --marker <str>    Path to the input marker file.
                        Specifies the subset of genes to be used for testing.
                        The file must be tab-delimited with a header line.
                        A column "Gene" is required.
  -o, --output <str>    Path to the output file.
  --max-mt <float>      Maximum mitochondrial gene percentage allowed in a cell (%).
                        The value must be a float between 0 and 100.
                        Default: 10.0.
  --min-genes <int>     Minimum number of detected genes in a cell.
                        The value must be a non-negative integer.
                        Default: 100.
  --max-genes <int>     Maximum number of detected genes in a cell.
                        The value must be a non-negative integer.
                        Default: 6000.
  --min-cells <int>     Minimum number of cells expressing a gene.
                        The value must be a non-negative integer.
                        Default: 10.
  --alpha <float>       Significance threshold for statistical tests.
                        The value must be a float between 0 and 1.
                        Default: 0.05.
  --multiple-test <str>
                        Method for multiple hypothesis testing.
                        The values can be bonferroni and benjamini-hochberg.
                        Default: benjamini-hochberg.
  --min-cells-group <int>
                        Minimum number of cells in a group.
                        The value must be a positive integer.
                        Default: 2.
```

### info

```text
usage: scTools info [options] -i H5AD -o INFO.

options:
  -h, --help          show this help message and exit
  -i, --input <str>   Path to the input h5ad file.
  -o, --output <str>  Prefix for the output files.
```

### cluster

```text
usage: scTools cluster [options] -i H5AD -g GROUP -o CLUSTER.

options:
  -h, --help            show this help message and exit
  -i, --input <str>     Path to the input h5ad file.
  -g, --group <str>     Path to the input group file.
                        The file must be tab-delimited with a header line.
                        The first column defines group names.
                        All other columns must match cell annotations in the H5AD file.
                        Use "scTools info" to inspect available columns in the H5AD file.
  -o, --output <str>    Path to the output file.
  --max-mt <float>      Maximum mitochondrial gene percentage allowed in a cell (%).
                        The value must be a float between 0 and 100.
                        Default: 10.0.
  --min-genes <int>     Minimum number of detected genes in a cell.
                        The value must be a non-negative integer.
                        Default: 100.
  --max-genes <int>     Maximum number of detected genes in a cell.
                        The value must be a non-negative integer.
                        Default: 6000.
  --min-cells <int>     Minimum number of cells expressing a gene.
                        The value must be a non-negative integer.
                        Default: 10.
  --tsne-perplexity <int>
                        The perplexity parameter in t-SNE.
                        The value must be a positive integer.
                        Default: 1.
  --dbscan-eps <float>  The eps parameter in DBSCAN.
                        The value must be a positive float.
                        Default: 5.
  --width <float>       Figure width in inches.
                        The value be a positive float.
                        Default: 7.0.
  --height <float>      Figure height in inches.
                        The value be a positive float.
                        Default: 7.0.
```

### marker

```text
usage: scTools marker [options] -t TEST -o MARKER.

options:
  -h, --help           show this help message and exit
  -d, --diffexp <str>  Path to the input differential analysis file.
  -o, --output <str>   Path to the output file.
  -m, --mapping <str>  Path to the input mapping file.
                       The file must be tab-delimited with a header line: "Gene<tab>Protein".
  --log2fc <float>     Minimum log2 fold-change of mean expression in a group compared to all other groups.
                       Default: 2.
  --rank <float>       Genes with expression rank greater than "--rank" are removed.
                       The value must be a float between 0 and 100.
                       Default: 10.
  --mean <float>       Genes with mean expression (log1p) lower than "--mean" are removed.
                       Default: 0.
  --coverage <int>     Genes identified as markers in more than "--coverage" of groups will be removed.
                       The value must be a positive integer.
                       Default: inf.
  --hits <int>         Genes are only considered marker genes in the top "--hits" groups.
                       The value must be a positive integer.
                       Default: inf.
```
