# scTools

A toolkit for single-cell sequencing data analysis in H5AD format.

## Install scTools

```text
python3 -m venv sctools
source sctools/bin/activate
pip3 install https://github.com/liu-congcong/scTools/releases/download/v1.0.0/sctools-1.0.0-py3-none-any.whl
deactivate
```

## Examples

### 1. Preparation

```text
wget -O kidney.h5ad https://datasets.cellxgene.cziscience.com/65ca6e36-73b0-4c88-b0f3-7b23b48844ad.h5ad
wget -O pancreas.h5ad https://datasets.cellxgene.cziscience.com/4be7951a-4ae7-4a81-86f1-e2649b304d1c.h5ad
```

### 2. Extract information from an H5AD file (< 1 min)

```text
scTools info -i kidney.h5ad -o kidney
scTools info -i pancreas.h5ad -o pancreas
```

### 3. Merge multiple files into a single file (< 1 min)

```text
scTools merge --label organ --obs-key cell_type -i kidney.h5ad pancreas.h5ad -o scTools-merge.h5ad
```

In `scTools-merge.h5ad`, only the `cell_type` annotation was retained, and a new annotation `organ` was added.

All cells originating from `kidney.h5ad` were assigned the label `kidney`, and all cells from `pancreas.h5ad` were assigned the label `pancreas`.

You can use `scTools info` to inspect the information of `scTools-merge.h5ad`.

### 4. Perform differential analysis of genes between groups from an H5AD file (< 1 min)

To perform differential analysis of genes across different organ - cell type combinations, a `group` file must be defined with the following format:

|  *  |organ|cell_type|
|:---:|:---:|:-------:|
| ... | ... |   ...   |

The first column specifies the group identifier, and all subsequent columns must match existing annotations in H5AD file.

An example `group` file is provided here: [`group.tsv`](./Examples/group.tsv).

By default, scTools performs differential analysis on all genes across groups.

To restrict the analysis to a subset of genes, a `marker` file can be provided:

| * |Gene| * |
|:-:|:--:|:-:|
|...|... |...|

The `marker` file must contain a column named `Gene` (case-sensitive).

An example `marker` file is provided here: [`marker.txt`](./Examples/marker.txt).

```text
scTools diffexp -i scTools-merge.h5ad -g group.tsv -m marker.txt -o scTools-diffexp.tsv
```

### 5. Identify group-specific marker genes (< 1 min)

```text
scTools marker -d scTools-diffexp.tsv -m mapping.tsv -o scTools-marker.tsv
```

## Documentation

[HERE](https://github.com/liu-congcong/scTools/tree/main/Documentation).

## Updates

### v1.0.0

The first release.
