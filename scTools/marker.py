def readDiffExpFile(file, minMean, maxRank, minLog2fc):
    gene2log2FCMeanRankGroups = dict()
    openFile = open(file, 'r')
    headers = openFile.readline().rstrip('\n').split('\t')
    group = headers.index('Group')
    gene = headers.index('Gene ID')
    mean = headers.index('Mean (log1p)')
    rank = headers.index('Rank (%)')
    log2fc = headers.index('Log2FC')
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        mean_ = float(lines[mean])
        rank_ = float(lines[rank])
        log2fc_ = float(lines[log2fc])
        if mean_ >= minMean and rank_ <= maxRank and log2fc_ >= minLog2fc:
            gene2log2FCMeanRankGroups.setdefault(lines[gene], list()).append((log2fc_, mean_, rank_, lines[group]))
    openFile.close()
    return gene2log2FCMeanRankGroups


def readMappingFile(file):
    gene2protein = dict()
    openFile = open(file, 'r')
    headers = openFile.readline().rstrip('\n').split('\t')
    gene = headers.index('Gene')
    protein = headers.index('Protein')
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        gene2protein[lines[gene]] = lines[protein]
    openFile.close()
    return gene2protein


def main(parameters):
    gene2log2FCMeanRankGroups = readDiffExpFile(parameters.diffexp, parameters.mean, parameters.rank, parameters.log2fc)
    if parameters.mapping is None:
        gene2protein = dict()
    else:
        gene2protein = readMappingFile(parameters.mapping)
    openFile = open(parameters.output, 'w')
    openFile.write('Group\tGene\tProtein\tMean (log1p)\tRank (%)\tLog2FC\n')
    for gene, x in gene2log2FCMeanRankGroups.items():
        if (parameters.coverage is None) or (len(x) <= parameters.coverage):
            x.sort(reverse = True) # log2FCMeanRankGroup, ..., log2FCMeanRankGroup #
            for log2FC, mean, rank, group in x[ : parameters.hits]:
                openFile.write(f'{group}\t{gene}\t{gene2protein.get(gene, 'NA')}\t{mean}\t{rank}\t{log2FC}\n')
    openFile.close()
    return None
