import sqlite3
import argparse
import csv
from pathlib import Path
import time
from utils import get_genes


def create_table(cur, header, pval_column):
    columns = []
    for col in header:
        if col == 'chr' or col == 'pos':
            dtype = 'INTEGER'
        elif col == pval_column:
            dtype = 'REAL'
        else:
            dtype = 'TEXT'
        columns.append(f'"{col}" {dtype}')
    create_sql = f"CREATE TABLE IF NOT EXISTS variants ({', '.join(columns)})"
    cur.execute(create_sql)


def process_value(value, col_name, pval_column):
    if value == 'NA':
        return None
    if col_name == 'chr' or col_name == 'pos':
        try:
            return int(value)
        except:
            return None
    elif col_name == pval_column:
        try:
            return float(value)
        except:
            return None
    return value


def import_tsv(cur, tsv_path, header, pval_column):
    with open(tsv_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # Skip header
        rows = []
        for row in reader:
            processed = [process_value(v, header[i], pval_column)
                         for i, v in enumerate(row)]
            rows.append(processed)
            if len(rows) >= 1000:
                cur.executemany("INSERT INTO variants VALUES ({})".format(
                    ','.join(['?']*len(header))), rows)
                rows = []
        if rows:
            cur.executemany("INSERT INTO variants VALUES ({})".format(
                ','.join(['?']*len(header))), rows)


def main(tsv_files, db_path, pval_column, timings_path=None, bed_path=None, pval=None):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # Get header from first file
    with open(tsv_files[0], 'r') as f:
        header = next(csv.reader(f, delimiter='\t'))

    create_table(cur, header, pval_column)

    for tsv in tsv_files:
        import_tsv(cur, tsv, header, pval_column)
        conn.commit()

    # Create index
    cur.execute(
        f"CREATE INDEX idx_chr_pos_pval ON variants (chr, pos, {pval_column})")
    conn.commit()

    print("Database created")

    try:

        # Perform query
        if not timings_path or not bed_path or not pval:
            print("Skipping query because missing timings_path, bed_path or pval")
            return

        print("Performing query")

        with open(timings_path, 'w') as f:
            base_name = Path(db_path).stem
            f.write(f"GWAS file: {base_name}\n")
            genes = get_genes(bed_path)

            for gene in genes:
                for chrom in genes[gene]:
                    for start, end in genes[gene][chrom]:
                        t0 = time.time()

                        cur.execute(
                            f"""SELECT * FROM variants
                            WHERE chr = "{chrom}"
                            AND pos >= {start}
                            AND pos <= {end}
                            AND {pval_column} {pval};""")
                        _ = cur.fetchall()  # results are discarded

                        t1 = time.time()
                        duration = t1 - t0
                        f.write(f'Gene: {gene},time: {duration}\n')

        print("Done")

    finally:
        conn.close()  # absolutely


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Build SQLite DB for GWAS data.')
    parser.add_argument('--tsv', nargs='+', required=True,
                        help='Input TSV files.')
    parser.add_argument('--db', required=True, help='Output SQLite DB path.')
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
