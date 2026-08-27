#!/bin/bash

###################################################################
# script name   : make.bash
# description   : This script runs the nf-core ampliseq workflow to
#                 generate an ASV table from metabarcoding raw reads.
#                 Extra steps run the taxonomic assignment and
#                 format the output.
# usage         : Wait for step to finish before launching another one.
# author	: Nicolas Henry
# contact	: nicolas.henry@cnrs.fr
###################################################################

########################################################
# Setup
########################################################

# create directories for outputs

[ ! -d "tmp" ] && mkdir tmp
[ ! -d "outputs/taxo_assignment/" ] && mkdir -p outputs/taxo_assignment/
[ ! -d "outputs/final_tables/" ] && mkdir -p outputs/final_tables/

ls input/manifest | \
    grep "_ampliseq.tsv" | \
    sed 's/_ampliseq.tsv//' > tmp/my_projects.txt

while read PROJECT; do

  [ ! -f "input/parameters/${PROJECT}_params_ampliseq.yaml" ] && cp input/parameters/template_params_ampliseq.yaml input/parameters/${PROJECT}_params_ampliseq.yaml

done < tmp/my_projects.txt

# list references databases from taxonomic assignment parameters files in input/parameters

ls input/parameters/ | \
    grep "_taxo_assign.yaml" | \
    sed 's/_taxo_assign.yaml//' > tmp/refdbs.txt

## bash variables to update

META_PROJECT="coastclim_spatial"
FAV_ASSIGN="vsearch"

# dowload reference databases

sbatch scripts/analyses/00_refdb_download_ekoi.bash

sbatch scripts/analyses/00_refdb_download_mzgdb.bash

sbatch scripts/analyses/00_refdb_download_bold.bash

########################################################
# Run ampliseq and format output
########################################################

sbatch \
    --output 01_asvs_per_project.out \
    scripts/analyses/01_asvs_per_project.bash

# manually check if ampliseq output is alright

sbatch \
    --output 02_format_tables_per_project.out \
    scripts/analyses/02_format_tables_per_project.bash

########################################################
# Taxonomic assignment
########################################################

sbatch \
    scripts/analyses/03_taxo_assignment.bash

########################################################
# Prepare tables for export
########################################################

# R has to be loaded before running the command below
module load r/4.4.1

sbatch \
    --mem 10GB \
    --output 04_format_table_export.out \
    --wrap="Rscript scripts/analyses/04_format_table_export.R ${META_PROJECT} ${FAV_ASSIGN}"
