- added `--name tabix/sqlite/hdf5` to several plotting methods
  - publish.py
  - plot_tabix.py
  - table.py
  - plot_compare.py
  - ryan_hist.py did not have a CLI so I included a `name = 'tabix'` variable.

- missing `/Users/krsc0813/PycharmProjects/gwas-analysis/data/UKBB/test-ukbb-manifest.csv`
  I infer this should be the lines corresponding to "the 10 PUKBB files" from the larger manifest.
- I had to guess the bin boundaries for the stabix queries (regarding "the ten pan ukbb files").
  I borrowed the boundaries from the default stabix config files from the main stabix repo.
 
environment
- mamba env included wrong py version. Need note in README to explicitly specify py (11?).
- why are there so many different conda_env ymls? Separation, particularly when each arent minimal, is overkill.
- the given environments were definitely intended for a linux x86 system. Missing osx env details.
