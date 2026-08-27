#!/bin/bash

#SBATCH --partition long
#SBATCH --cpus-per-task 1
#SBATCH --mem 5GB
#SBATCH --output 00_refdb_download_mzgdb.out

########################################
# mzgdb
########################################

module load nextflow/23.04.1
module load r/4.3.1

# downloaded 2025-04-10

wget --no-clobber https://metazoogene.org/mzgdb/atlas/data-src/MZGdata-coi__T4000000__o00__C.psv.gz \
  --directory-prefix archive/refdb/mzgdb/

zcat archive/refdb/mzgdb/MZGdata-coi__T4000000__o00__C.psv.gz | \
  awk -F "\|" '{sub(/;$/,"",$34); print ">" $9 " " $34 "\n" toupper($31)}' | \
  gzip > archive/refdb/mzgdb/MZGdata-coi__T4000000__o00__C.fasta.gz

nextflow \
    run https://gitlab.com/metabarcoding_utils/metab-refdb-importer \
    -r v0.0.6 \
    -c input/config/refdb_import_abims.config \
    -params-file input/parameters/mzgdb_import_params.yaml
