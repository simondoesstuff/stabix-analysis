import argparse
import os
import time

import h5py
import numpy as np
from typing import List, Tuple, Dict
import csv

from utils import get_column_names, get_genes, guess_pval_col


def setup_dataset(h5_path: str, columns: List[str], pval_col) -> None:
    """
    Create or open HDF5 file and set up datasets for genomic data.

    Args:
        h5_path: Path to HDF5 file
        columns: List of column names from TSV header
    """
    # Verify required columns exist
    required_cols = {'chrm', 'pos', pval_col}
    if not required_cols.issubset(set(columns)):
        missing = required_cols - set(columns)
        raise ValueError(f"Missing required columns: {missing}")

    # INFO: only storing indexed columns for performance reasons
    target_columns = list()
    for col in columns:
        if col == pval_col:
            target_columns.append((col, 'float64'))
        elif col == 'pos':
            target_columns.append((col, 'int32'))
        elif col == 'chrm':
            target_columns.append((col, h5py.string_dtype(encoding='utf-8')))

    with h5py.File(h5_path, 'w') as f:
        # Create a group for metadata
        meta = f.create_group('metadata')
        # Store column names as attribute
        meta.attrs['columns'] = columns

        # Create extendable datasets for each column
        variants = f.create_group('variants')

        for name, dtype in target_columns:
            variants.create_dataset(name, (0,), maxshape=(None,), dtype=dtype, compression='gzip')


def load_data(h5_path: str, tsv_path: str, pval_col) -> None:
    """
    Load data from TSV file into the HDF5 file.
    """
    # Get column names and set up file
    columns = get_column_names(tsv_path)
    setup_dataset(h5_path, columns, pval_col)

    # Find indices of required columns
    col_indices = {
        'chrm': columns.index('chrm'),
        'pos': columns.index('pos'),
        pval_col: columns.index(pval_col)
    }

    # INFO: skipping non-indexed columns for performance reasons
    target_columns = [ 'chrm', 'pos', pval_col ]

    # Read all data first to determine array sizes
    data: Dict[str, List] = {col: [] for col in columns}

    print("... Parsing input")
    with open(tsv_path, 'r') as tsv_file:
        next(tsv_file)  # Skip header
        tsv_reader = csv.reader(tsv_file, delimiter='\t')

        for i, row in enumerate(tsv_reader):
            if i % 10000 == 0:
                print(f"line {i}")

            pval_idx = col_indices[pval_col]

            if pval_idx >= len(row):
                continue

            pval = row[pval_idx]

            if pval == 'NA':
                continue

            try:
                # Store each column's data
                for col in target_columns:
                    idx = col_indices[col]
                    val = row[idx]

                    # Convert types for required numeric columns
                    if col == 'pos':
                        val = int(val)
                    elif col == pval_col:
                        val = val.upper().replace('EE', 'E')  # Fix scientific notation
                        val = float(val)

                    data[col].append(val)
            except IndexError as e:
                print(f"Warning: Error processing row {row}: {e}")
                continue

    # Write data to HDF5 file
    print("... Writing index")
    with h5py.File(h5_path, 'a') as tsv_file:
        variants = tsv_file['variants']
        n_rows = len(data['chrm'])

        # Resize and write each dataset
        for i, col in enumerate(target_columns):
            print(f"col: {col}")
            dataset = variants[col]
            dataset.resize((n_rows,))

            if col in {'pos', pval_col}:
                # Numeric arrays
                dataset[:] = np.array(data[col])
            else:
                # String arrays
                dataset[:] = np.array(data[col], dtype=h5py.string_dtype(encoding='utf-8'))


def query_interval(
        h5_path: str,
        chrm: str,
        start_pos: int,
        end_pos: int,
        max_pval: float,
        pval_col
) -> Tuple[List[Tuple], List[str]]:
    """
    Query variants within a genomic interval that meet the p-value threshold.

    Args:
        h5_path: Path to HDF5 file
        chrm: Chromosome name
        start_pos: Start position of interval
        end_pos: End position of interval
        max_pval: Maximum p-value threshold

    Returns:
        Tuple of (results list, column names list)
    """
    with h5py.File(h5_path, 'r') as f:
        variants = f['variants']
        columns = list(f['metadata'].attrs['columns'])
        # INFO: skipping non-indexed columns for performance reasons
        target_columns = [ 'chrm', 'pos', pval_col ]

        # Convert chromosome strings for comparison
        chrm_array = variants['chrm'][:]
        chrm_match = chrm_array.astype(str) == chrm

        # Get positions and p-values
        pos_array = variants['pos'][:]
        pval_array = variants[pval_col][:]

        # Create boolean masks for each condition
        pos_match = (pos_array >= start_pos) & (pos_array <= end_pos)
        pval_match = pval_array <= max_pval

        # Combine all conditions
        mask = chrm_match & pos_match & pval_match

        # Get indices where mask is True
        indices = np.where(mask)[0]

        # Collect results
        results = []
        for idx in indices:
            row = []
            for col in target_columns:
                val = variants[col][idx]
                if isinstance(val, np.bytes_):
                    val = val.decode('utf-8')
                row.append(val)
            results.append(tuple(row))

        # Sort results by position
        # results.sort(key=lambda x: x[columns.index('pos')])

        return results, columns


def print_database(h5_path: str) -> None:
    with h5py.File(h5_path, 'r') as f:
        variants = f['variants']
        columns = list(f['metadata'].attrs['columns'])

        # Get all data
        data = []
        for i in range(len(variants['chrm'])):
            row = []
            for col in columns:
                val = variants[col][i]
                if isinstance(val, np.bytes_):
                    val = val.decode('utf-8')
                row.append(val)
            data.append(row)

        # Sort by chromosome and position
        pos_idx = columns.index('pos')
        chrm_idx = columns.index('chrm')
        data.sort(key=lambda x: (x[chrm_idx], x[pos_idx]))

        print("\nDatabase Contents:")
        print("-" * 80)
        print("Total rows:", len(data))
        print("-" * 80)

        # Calculate column widths
        widths = {}
        for i, col in enumerate(columns):
            values = [str(row[i]) for row in data]
            values.append(col)
            widths[col] = max(10, max(len(str(val)) for val in values) + 2)

        # Print column headers
        header = "".join(col.ljust(widths[col]) for col in columns)
        print(header)
        print("-" * len(header))

        # Print rows
        for row in data:
            print("".join(str(val).ljust(widths[col]) for val, col in zip(row, columns)))

        print("-" * 80)

        # Print file structure
        print("\nHDF5 File Structure:")
        def print_structure(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f"Dataset: {name}")
                print(f"  Shape: {obj.shape}")
                print(f"  Type: {obj.dtype}")
                print(f"  Compression: {obj.compression}")
            elif isinstance(obj, h5py.Group):
                print(f"Group: {name}")
                if hasattr(obj, 'attrs') and obj.attrs:
                    print("  Attributes:")
                    for key, val in obj.attrs.items():
                        print(f"    {key}: {val}")

        f.visititems(print_structure)


def parse_args():
    parser = argparse.ArgumentParser(description='HDF5 GWAS query')
    parser.add_argument('--bed', type=str, help='gene bed file')
    parser.add_argument('--gwas', type=str, help='GWAS file (.bgz)')
    parser.add_argument('--pval_threshold', type=str, help='greater than or equal to this value')
    parser.add_argument('--out', type=str, help='test_output directory for results')
    return parser.parse_args()


def main():
    args = parse_args()
    gwas = args.gwas
    H5_PATH = gwas + '.h5'
    TSV_PATH = gwas
    pval = float(args.pval_threshold)

    columns = get_column_names(TSV_PATH)
    pval_col = guess_pval_col(columns)

    print("... Indexing.")
    load_data(H5_PATH, TSV_PATH, pval_col)
    print("Index created.")

    output_tabix_query_file = os.path.join(args.out)
    out_file = open(output_tabix_query_file, 'a')
    out_file.truncate(0)
    gwas_file_basename = os.path.basename(gwas).replace('.tsv', '')
    out_file.write('GWAS file: {}\n'.format(gwas_file_basename))

    genes = get_genes(args.bed)
    for i, gene in enumerate(genes):
        if i % 100 == 0:
            print(f"{i} / {len(genes)}")
        for chrom in genes[gene]:
            for start, end in genes[gene][chrom]:
                start_time = time.time()

                # Results are discarded
                results, columns = query_interval(
                    h5_path=H5_PATH,
                    chrm=chrom,
                    start_pos=start,
                    end_pos=end,
                    max_pval=pval,
                    pval_col=pval_col)

                end_time = time.time()
                duration = end_time - start_time
                out_file.write('Gene: {},time: {}\n'.format(gene, duration))

    os.remove(H5_PATH)

if __name__ == "__main__":
    main()
