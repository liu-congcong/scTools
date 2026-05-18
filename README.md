# SCTools

A toolkit for single-cell sequencing data analysis in H5AD format.

## Install SCTools

```text
python3 -m venv sctools
source sctools/bin/activate
pip3 install https://github.com/liu-congcong/SCTools/releases/download/v1.0.0/sctools-1.0.0-py3-none-any.whl
deactivate
```

## Examples

### 1. Data preparation

```text
wget https://datasets.cellxgene.cziscience.com/65ca6e36-73b0-4c88-b0f3-7b23b48844ad.h5ad
mv 65ca6e36-73b0-4c88-b0f3-7b23b48844ad.h5ad kidney.h5ad
wget https://datasets.cellxgene.cziscience.com/4be7951a-4ae7-4a81-86f1-e2649b304d1c.h5ad
mv 4be7951a-4ae7-4a81-86f1-e2649b304d1c.h5ad pancreas.h5ad
```

### 2. Extract information from an H5AD file

```text
scTools info -i kidney.h5ad -o kidney
scTools info -i pancreas.h5ad -o pancreas
```

### 3. Merge multiple files into a single file

```text
scTools merge --label organ --obs-key cell_type -i kidney.h5ad pancreas.h5ad -o scTools-merge.h5ad
```

### 4. Perform differential analysis of genes between groups from an H5AD file

```text
scTools diffexp -i scTools-merge.h5ad -g group.tsv -m marker.tsv -o scTools-diffexp
```

### 5. Identify group-specific marker genes

```text
scTools marker -d scTools-diffexp -m mapping.tsv -o scTools-marker
```

## Documentation

[HERE](https://github.com/liu-congcong/SCTools/tree/main/Documentation).

## Updates

### v1.0.0

The first release.
