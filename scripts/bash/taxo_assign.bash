#!/bin/bash

#SBATCH --partition long
#SBATCH --cpus-per-task 1
#SBATCH --mem 5GB

module load nextflow/24.04.4

nextflow \
    run https://gitlab.com/metabarcoding_utils/taxonomic-assignment \
    -r 0.2.4 \
    -profile abims_config \
    -c input/config/bold_assign.config \
    -params-file input/parameters/${REFDB}_taxo_assign.yaml