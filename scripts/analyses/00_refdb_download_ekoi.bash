#!/bin/bash

#SBATCH --partition long
#SBATCH --cpus-per-task 1
#SBATCH --mem 5GB
#SBATCH --output 00_refdb_download_ekoi.out

########################################
# ekoi
########################################

module load nextflow/23.04.1
module load r/4.3.1

wget --no-clobber https://github.com/rubenmiguens/eKOI_taxonomy_database/raw/99dbf5c407de920056274eca1db53cffc24faef3/eKOI.fasta \
  --directory-prefix archive/refdb/ekoi/

awk '
  />/ {
    sub(/>/, "", $1)
    id = taxo = $1
    sub(/^.+;/, "", id)
    sub(/;[^;]+$/, "", taxo)
    print ">"id" "taxo
    next
  }
  {print}
' archive/refdb/ekoi/eKOI.fasta | gzip > archive/refdb/ekoi/eKOI_formated.fasta.gz

nextflow \
    run https://gitlab.com/metabarcoding_utils/metab-refdb-importer \
    -r v0.0.6 \
    -c input/config/refdb_import_abims.config \
    -params-file input/parameters/ekoi_import_params.yaml
