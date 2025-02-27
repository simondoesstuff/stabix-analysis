from typing import List


def guess_pval_col(all_columns):
    pval_col = None

    # infer pval column in a greedy way
    for col in all_columns:
        if 'pval' in col:
            pval_col = col
            break

    if not pval_col:
        raise RuntimeError("Unable to infer pval column. Didn't find a 'pval' in any column names.")

    return pval_col


def get_column_names(tsv_file) -> List[str]:
    if type(tsv_file) == str:
        f = open(tsv_file, 'r')
    else:
        f = tsv_file

    header = next(f).strip().split()
    # 'chr' -> 'chrm'
    header = list(map(lambda x: 'chrm' if x == 'chr' else x, header))
    f.close()
    return header


def get_genes(bed_file):
    genes = {}

    with open(bed_file, 'r') as f:
        for line in f:
            chrom, start, end, gene = line.strip().split()[:4]
            if gene not in genes:
                genes[gene] = {}
            if chrom not in genes[gene]:
                genes[gene][chrom] = []
            genes[gene][chrom].append((int(start), int(end)))

    f.close()
    return genes
