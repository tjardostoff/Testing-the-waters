# Testing the waters

## Presentation

This repository gathers all the code necessary to assemble the ASV table analyzes in the study "Testing the water". It also comprised the code use to perform the statistical analyses.

## How to run the code

All the scripts to assemble the ASV table are called by `main.bash`. For statistical analyses the script 

## Structure

### input

#### manifest

* `{project_name}_ampliseq.tsv`: nf-core/ampliseq input file, indicating for each sample, the corresponding run, the forward file and reverse file for `{project_name}`.  

to remove once template filled: `input/manifest/project1_ampliseq.tsv` gives an example for sequencing replicates (`sample1_1` and `sample1_2` are sequencing replicates of sample1) and `input/manifest/project2_ampliseq.tsv` gives an example for mixed orientated reads

#### Nextflow config files

* `ampliseq.config`: configuration file to run [nf-core/ampliseq](https://github.com/nf-core/ampliseq/) on ABiMS cluster

#### Nextflow parameters files

* `template_params_ampliseq.yaml`: list of parameters used with nf-core/ampliseq. During setup, this file will be copied for each project and each copy can be modified to adapt to each project specificity
* `{refdb}_taxo_assign.yaml`: list of parameters used for the taxonomic assignment using {refdb} (multiple files when multiple reference databases)

### reports

* `readme_zenodo.qmd`: quarto file describing the resulting datasets
* `references.bib`: list of references cited in `readme_zenodo.qmd`

### scripts

#### analyses

* `01_asvs_per_project.bash`: run nf-core/ampliseq for each project (one file in `input/manifest` per project)
* `02_format_tables_per_project.bash`: script transforming counts in long format for each project
* `03_taxo_assignment.bash`: run the taxonomic assignment for all the projects together
* `04_format_table_export.R`: script formatting files for export

#### awk

* `asv_table_longer.awk`: awk script called by `03_long_format_counts.bash` to transform ASV counts in long format

#### bash

* `ampliseq.bash`: script called by `01_asvs_per_project.bash` to run nf-core/ampliseq on a given project raw data to generate an ASV table
* `subset_and_format_table.bash`: script called by `02_format_tables_per_project.bash` to transform a project ASV table into long format counts (see `asv_table_longer.awk`), filter out rare ASVs (< 3 reads or < 2 samples) and export filtered fasta file