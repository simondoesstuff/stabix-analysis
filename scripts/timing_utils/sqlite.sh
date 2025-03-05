#!/bin/bash

get_to_work() {
  # ❯ python3 scripts/python/sqlite_hdf5/sqlite_query.py  --help
  # usage: sqlite_query.py [-h] --tsv TSV [TSV ...] --db DB --pval_col PVAL_COL
  #                        [--timings TIMINGS] [--bed BED] [--pval PVAL]
  #
  # Build SQLite DB for GWAS data.
  #
  # options:
  #   -h, --help           show this help message and exit
  #   --tsv TSV [TSV ...]  Input TSV files.
  #   --db DB              Output SQLite DB path.
  #   --pval_col PVAL_COL  P-value column name.
  #   --timings TIMINGS    Output timings data path.
  #   --bed BED            Input bed data path.
  #   --pval PVAL          P-value condition. eg: ">= 7.3"

  gwas="data/gwas"
  db="$gwas/$1.db"

  [ -f $db ] && rm $db # remove iff exists

  echo -e "\tPreparing $db..."
  python3 scripts/python/sqlite_hdf5/sqlite_query.py \
    --tsv "$gwas/$1.tsv" \
    --db $db \
    --pval_col $2 \
    --timings tabix_output/$1_sqlite_output.txt \
    --bed data/bed_files/hg19.protein_coding.bed \
    --pval ">= 7.3"
}

pval="neglog10_pval_meta_hq"
get_to_work categorical-1210-both_sexes-1210 $pval
get_to_work categorical-20096-both_sexes-2 $pval
get_to_work categorical-20116-both_sexes-0 $pval
get_to_work phecode-250-both_sexes-1210 $pval
get_to_work phecode-250-both_sexes $pval
get_to_work phecode-250.2-both_sexes $pval
get_to_work phecode-282.5-both_sexes $pval
get_to_work phecode-696.4-both_sexes $pval
pval="neglog10_pval_EAS"
get_to_work continuous-103220-both_sexes $pval
get_to_work continuous-50-both_sexes-irnt $pval
get_to_work continuous-30130-both_sexes-irnt $pval

echo "Celebrate."
