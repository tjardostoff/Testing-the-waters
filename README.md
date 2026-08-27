# Testing the waters

---

## Presentation

This repository contains all the code required to **assemble ASV tables** and perform **statistical analyses** for the study *"Testing the waters"*.

---

## How to Run the Code

1. **Assemble ASV tables**: Run the main script to process all projects:
  ```bash
   bash main.bash
  ```
2. **Run statistical analyses**: Execute the R script independently:
  ```bash
   Rscript scripts/analyses/05_script_statistical_analyses.R
  ```

Configuration files need to be adapted to work in a different environment than the ABiMS SLURM cluster.

---

## Dependencies

- [Nextflow](https://www.nextflow.io/)
- [nf-core/ampliseq](https://github.com/nf-core/ampliseq/)
- R packages: `dplyr`, `vegan`, `tidyverse`

---

## Structure

---

### `input/`

#### Manifest Files

- `spatial_ampliseq.tsv`: TSV file (nf-core/ampliseq input) mapping each sample to its corresponding run, forward file, and reverse file.

#### Nextflow Configuration

- `ampliseq.config`: Configuration file to run [nf-core/ampliseq](https://github.com/nf-core/ampliseq/) on the ABiMS cluster.

#### Taxonomic Assignment Configuration

- `refdb_import_abims.config`: Configuration file to import reference databases for taxonomic assignment.
- `bold_assign.config`: Configuration file to run taxonomic assignment against BOLD.

#### Parameters

- `spatial_params_ampliseq.yaml`: Parameters for running nf-core/ampliseq.
- `[refdb]_import_params.yaml`: Parameters for importing the `[refdb]` reference database.
- `[refdb]_taxo_assign.yaml`: Parameters for taxonomic assignment using `[refdb]`.

---

### `scripts/`

#### `analyses/`

- `00_refdb_download_bold.bash`: Imports the BOLD reference database.
- `00_refdb_download_ekoi.bash`: Imports the EKOI reference database.
- `00_refdb_download_mzgdb.bash`: Imports the mzgdb reference database.
- `01_asvs_per_project.bash`: Runs nf-core/ampliseq for each project (one file per project in `input/manifest`).
- `02_format_tables_per_project.bash`: Transforms ASV counts into long format and filters out rare ASVs (< 3 reads or < 2 samples).
- `03_taxo_assignment.bash`: Runs taxonomic assignment for all projects.
- `04_format_table_export.R`: Formats and exports tables for downstream analyses.
- `05_script_statistical_analyses.R`: Performs ecological and statistical analyses.

#### `awk/`

- `asv_table_longer.awk`: AWK script to transform ASV counts into long format (called by `02_format_tables_per_project.bash`).

#### `bash/`

- `ampliseq.bash`: Runs nf-core/ampliseq on raw data for a given project (called by `01_asvs_per_project.bash`).
- `subset_and_format_table.bash`: Filters ASVs and exports filtered FASTA files (called by `02_format_tables_per_project.bash`).
- `taxo_assign.bash`: Runs taxonomic assignment for one reference database (called by `03_taxo_assignment.bash`).

---

## Outputs

- `results/asv_tables/`: ASV tables per project (filtered and in long format).
- `results/taxo_assignment/`: Taxonomic assignment results.
- `output/formatted_tables/`: Formatted tables for export.
- `output/figures/`: Plots and visualizations from statistical analyses.

---
