import os
import sqlite3
import time
from typing import List, Tuple
import csv
import argparse

from utils import get_column_names, get_genes, guess_pval_col


def setup_database(db_path: str, columns: List[str], required_cols: List[str]) -> None:
    # Verify required columns exist
    if not set(required_cols).issubset(set(columns)):
        missing = required_cols - set(columns)
        raise ValueError(f"Missing required columns: {missing}")

    chrm_col = required_cols[0]
    pos_col = required_cols[1]
    pval_col = required_cols[2]

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Build CREATE TABLE statement dynamically
    # Start with required columns
    col_defs = []

    for column in columns:
        if column == pos_col:
            col_defs.append(f'{column} INTEGER')
        elif column == pval_col:
            col_defs.append(f'{column} REAL')
        else:
            col_defs.append(f'{column} TEXT')

    create_table_sql = f'''
    CREATE TABLE IF NOT EXISTS variants (
        {','.join(col_defs)},
        UNIQUE({chrm_col}, {pos_col})
    )
    '''

    cursor.execute(create_table_sql)

    # Create indexes for faster queries
    cursor.execute(f'CREATE INDEX IF NOT EXISTS chr_pos_idx ON variants({chrm_col}, {pos_col})')
    cursor.execute(f'CREATE INDEX IF NOT EXISTS pval_idx ON variants({pval_col})')

    conn.commit()
    conn.close()


def load_data(db_path: str, tsv_path: str, columns: List[str], required_cols: List[str]) -> None:
    setup_database(db_path, columns, required_cols)

    chrm_col = required_cols[0]
    pos_col = required_cols[1]
    pval_col = required_cols[2]

    chrm_idx = columns.index(chrm_col)
    pos_idx = columns.index(pos_col)
    pval_idx = columns.index(pval_col)

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    with open(tsv_path, 'r') as tsv_file:
        # Skip header
        next(tsv_file)
        tsv_reader = csv.reader(tsv_file, delimiter='\t')

        # Prepare data for batch insert
        data = []
        for i, row in enumerate(tsv_reader):
            # Convert position to int and pval to float
            row_data = list(row)  # Convert to list to allow modification
            if len(row_data) != 3: # ignore problem columns
                continue
            row_data[pos_idx] = int(row_data[chrm_idx])
            pval_literal = row_data[pval_idx]
            pval_literal = pval_literal.upper().replace('EE', 'E')  # Fix scientific notation
            pval = float(pval_literal) if pval_literal != 'NA' else pval_literal
            # pval = float(pval_literal) if pval_literal != 'NA' else 0
            row_data[pval_idx] = pval

            try:
                data.append(tuple(row_data))
            except IndexError as e:
                print(f"Warning: Error processing row {row}: {e}")
                continue

        # Create the INSERT statement dynamically
        placeholders = ','.join(['?' for _ in columns])
        insert_sql = f'''
            INSERT OR REPLACE INTO variants 
            ({','.join(columns)})
            VALUES ({placeholders})
        '''

        # Batch insert
        cursor.executemany(insert_sql, data)

    conn.commit()
    conn.close()


def query_interval(
        db_path: str,
        chrm: str,
        start_pos: int,
        end_pos: int,
        max_pval: float,
        required_columns: List[str]
) -> List[Tuple]:
    chrm_col = required_columns[0]
    pos_col = required_columns[1]
    pval_col = required_columns[2]

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute(f'''
        SELECT *
        FROM variants
        WHERE {chrm_col} = ?
        AND {pos_col} >= ?
        AND {pos_col} <= ?
        AND {pval_col} <= ?
        ORDER BY {pos_col}
    ''', (chrm, start_pos, end_pos, max_pval))

    results = cursor.fetchall()

    # Get column names for reference
    cursor.execute('PRAGMA table_info(variants)')
    columns = [col[1] for col in cursor.fetchall()]

    conn.close()
    return results, columns


def print_database(db_path: str) -> None:
    """
    Print the entire contents of the database in a formatted way.
    Useful for debugging.

    Args:
        db_path: Path to SQLite database file
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Get column names
    cursor.execute('PRAGMA table_info(variants)')
    columns = [col[1] for col in cursor.fetchall()]

    # Get all data
    cursor.execute('SELECT * FROM variants ORDER BY chrm, pos')
    rows = cursor.fetchall()

    # Print header
    print("\nDatabase Contents:")
    print("-" * 80)
    print("Total rows:", len(rows))
    print("-" * 80)

    # Calculate column widths (minimum 10 characters, or length of longest entry + 2)
    widths = {}
    for i, col in enumerate(columns):
        values = [str(row[i]) for row in rows]
        values.append(col)  # Include column name in width calculation
        widths[col] = max(10, max(len(str(val)) for val in values) + 2)

    # Print column headers
    header = "".join(col.ljust(widths[col]) for col in columns)
    print(header)
    print("-" * len(header))

    # Print rows
    for row in rows:
        print("".join(str(val).ljust(widths[col]) for val, col in zip(row, columns)))

    print("-" * 80)

    # Print table schema
    print("\nTable Schema:")
    cursor.execute('PRAGMA table_info(variants)')
    for col in cursor.fetchall():
        print(f"Column: {col[1]}, Type: {col[2]}")

    # Print indexes
    print("\nIndexes:")
    cursor.execute("SELECT name, sql FROM sqlite_master WHERE type='index'")
    for index in cursor.fetchall():
        print(f"Index: {index[0]}")
        print(f"SQL: {index[1]}")

    conn.close()


def parse_args():
    parser = argparse.ArgumentParser(description='Sqlite GWAS query')
    parser.add_argument('--bed', type=str, help='gene bed file')
    parser.add_argument('--gwas', type=str, help='GWAS file (.bgz)')
    parser.add_argument('--pval_threshold', type=str, help='greater than or equal to this value')
    parser.add_argument('--out', type=str, help='test_output directory for tabix results')
    return parser.parse_args()


def main():
    args = parse_args()
    gwas = args.gwas
    DB_PATH = gwas + '.db'
    TSV_PATH = gwas
    pval = float(args.pval_threshold)

    all_columns = get_column_names(TSV_PATH)
    required_cols = ["chrm", "pos", guess_pval_col(all_columns)]

    # prepare database
    load_data(DB_PATH, TSV_PATH, all_columns, required_cols)

    output_tabix_query_file = os.path.join(args.out)
    out_file = open(output_tabix_query_file, 'a')
    out_file.truncate(0)
    gwas_file_basename = os.path.basename(gwas).replace('.tsv', '')
    out_file.write('GWAS file: {}\n'.format(gwas_file_basename))

    genes = get_genes(args.bed)
    for gene in genes:
        for chrom in genes[gene]:
            for start, end in genes[gene][chrom]:
                start_time = time.time()

                # Results are discarded
                results, columns = query_interval(
                    db_path=DB_PATH,
                    chrm=chrom,
                    start_pos=start,
                    end_pos=end,
                    max_pval=pval,
                    required_columns=required_cols)

                end_time = time.time()
                duration = end_time - start_time
                out_file.write('Gene: {},time: {}\n'.format(gene, duration))

if __name__ == "__main__":
    main()