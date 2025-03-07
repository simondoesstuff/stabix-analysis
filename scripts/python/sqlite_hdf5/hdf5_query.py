import pandas as pd
import os
import time
from utils import get_genes
import argparse
from pathlib import Path


def build_gwas_hdf5(tsv_path, hdf5_path, pval_col, chunksize=10000):
    """Convert GWAS TSV to optimized HDF5 format with chromosome-specific tables."""
    with pd.HDFStore(hdf5_path, mode='w') as store:
        for chunk in pd.read_csv(tsv_path, sep='\t', chunksize=chunksize,
                                 na_values=['NA'], keep_default_na=False,
                                 # Keep chr as string
                                 dtype={'chr': 'category'},
                                 true_values=['true'],
                                 false_values=['false']):
            # Process each column type
            for col in chunk.columns:
                if col == 'pos':
                    chunk[col] = pd.to_numeric(
                        chunk[col], errors='coerce').astype('int32')
                elif col == pval_col:
                    chunk[col] = pd.to_numeric(
                        chunk[col], errors='coerce').astype('float32')
                else:
                    chunk[col] = chunk[col].astype(str).str.slice(0, 6)

            # Group and store by chromosome
            for chr_name, group in chunk.groupby('chr', observed=True):
                key = f"/chr{chr_name.strip().upper()}"
                group = group.drop(columns=['chr']).sort_values('pos')

                store.append(
                    key,
                    group,
                    format='table',
                    index=False,
                    # chr not included its the top level key
                    data_columns=['pos', pval_col],
                    min_itemsize={
                        col: 6 for col in group.select_dtypes('string').columns})


def query_one(store, chrom, start, end, pval_col, pval):
    chrom = str(chrom).strip().upper()
    key = f"/chr{chrom}"

    if key not in store:
        return pd.DataFrame()

    query = (
        f"pos >= {start} & "
        f"pos <= {end} & "
        f"{pval_col}")

    return store.select(
        key,
        where=query)


def query_many(db_path, bed_path, timings_path, pval_col, pval):
    """Query GWAS data with chromosome-aware filtering."""
    total = 0

    with pd.HDFStore(db_path, mode='r') as store:
        with open(timings_path, 'w') as f:
            base_name = Path(db_path).stem
            f.write(f"GWAS file: {base_name}\n")
            genes = get_genes(bed_path)

            for gene in genes:
                for chrom in genes[gene]:
                    for start, end in genes[gene][chrom]:
                        t0 = time.time()

                        results = query_one(
                            store, chrom, start, end, pval_col, pval)

                        if not results.empty:
                            total += len(results)

                        t1 = time.time()
                        duration = t1 - t0
                        f.write(f'Gene: {gene},time: {duration}\n')

    print(f"Done. {total} row(s) answered query.")


def main(tsv_file, db_path, pval_column, timings_path=None, bed_path=None, pval=None):
    if os.path.exists(db_path):
        print(f"Database already exists: {db_path}")
        exit(1)

    build_gwas_hdf5(tsv_file, db_path, pval_column)
    print("Database created")

    if not timings_path or not bed_path or not pval:
        print("Skipping query because missing timings_path, bed_path or pval")
        return

    print("Performing query")
    query_many(db_path, bed_path, timings_path, pval_column, pval)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Build HDF5 DB for GWAS data.')
    parser.add_argument('--tsv',  required=True,
                        help='Input TSV file.')
    parser.add_argument('--db', required=True, help='Output .h5 path.')
    parser.add_argument('--pval_col', required=True,
                        help='P-value column name.')
    parser.add_argument('--timings', required=False,
                        help='Output timings data path.')
    parser.add_argument('--bed', required=False,
                        help='Input bed data path.')
    parser.add_argument('--pval', required=False,
                        help='P-value condition. eg: ">= 7.3"')
    args = parser.parse_args()
    main(args.tsv, args.db, args.pval_col, args.timings, args.bed, args.pval)
