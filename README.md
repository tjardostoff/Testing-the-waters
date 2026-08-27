# Testing the waters

## Presentation

This repository gathers all the code necessary to assemble the ASV table analyzes in the study "Testing the water". It also comprised the code use to perform the statistical analyses.

## How to run the code

All the scripts to assemble the ASV table are called by `main.bash`. For statistical analyses, run the script `scripts/analyses/05_script_statistical_analyses.R` independently.

## Structure

### input

#### manifest

* `spatial_ampliseq.tsv`: nf-core/ampliseq input file, indicating for each sample, the corresponding run, the forward file and reverse file.  

#### Nextflow config files

* `ampliseq.config`: configuration file to run [nf-core/ampliseq](https://github.com/nf-core/ampliseq/) on ABiMS cluster

#### Taxonomic assignment config files

* `refdb_import_abims.config`: configuration file to import the reference databases used for taxonomic assignment
* `bold_assign.config`: configuration file to run the taxonomic assignment against BOLD.

#### Nextflow parameters files

* `spatial_params_ampliseq.yaml`: list of parameters used with nf-core/ampliseq.
* `{refdb}_import_params.yaml`: list of parameters used for the import of {refdb}
* `{refdb}_taxo_assign.yaml`: list of parameters used for the taxonomic assignment using {refdb}

### scripts

#### analyses

* `00_refdb_download_bold.bash`: import BOLD reference database
* `00_refdb_download_ekoi.bash`: import EKOI reference database
* `00_refdb_download_mzgdb.bash`: import mzgdb reference database
* `01_asvs_per_project.bash`: run nf-core/ampliseq for each project (one file in `input/manifest` per project)
* `02_format_tables_per_project.bash`: script transforming counts in long format for each project
* `03_taxo_assignment.bash`: run the taxonomic assignment for all the projects together
* `04_format_table_export.R`: script formatting files for export
* `05_script_statistical_analyses.R`: run the ecological/statistical analyses

#### awk

* `asv_table_longer.awk`: awk script called by `03_long_format_counts.bash` to transform ASV counts in long format

#### bash

* `ampliseq.bash`: script called by `01_asvs_per_project.bash` to run nf-core/ampliseq on a given project raw data to generate an ASV table
* `subset_and_format_table.bash`: script called by `02_format_tables_per_project.bash` to transform a project ASV table into long format counts (see `asv_table_longer.awk`), filter out rare ASVs (< 3 reads or < 2 samples) and export filtered fasta file
* `taxo_assign.bash`: script called by `03_taxo_assignment.bash`, run the taxonomic assignment for one reference data base