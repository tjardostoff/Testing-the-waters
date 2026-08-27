#!/bin/bash

#SBATCH --partition long
#SBATCH --cpus-per-task 1
#SBATCH --mem 5GB
#SBATCH --output 00_refdb_download_bold.out

########################################
# bold
########################################

module load nextflow/23.04.1
module load r/4.3.1

zcat archive/refdb/bold/BOLD_Public.04-Apr-2025.fasta.gz | \
  awk '
    /^>/ { sub(/,[^,]+$/, ""); gsub(/ /, "_"); sub(/\|.+\|/, " ") }
    {print}
    
  ' | \
  gzip > archive/refdb/bold/BOLD_Public.04-Apr-2025_preformated.fasta.gz

nextflow \
    run https://gitlab.com/metabarcoding_utils/metab-refdb-importer \
    -r v0.0.6 \
    -c input/config/refdb_import_abims.config \
    -params-file input/parameters/bold_import_params.yaml
