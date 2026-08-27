#!/bin/bash

#SBATCH --partition long
#SBATCH --cpus-per-task 1
#SBATCH --mem 5GB

module load nextflow/24.04.4
module load singularity

# run metab pipeline

nextflow \
  run nf-core/ampliseq \
  -r 2.13.0 \
  -profile singularity,abims \
  -params-file "input/parameters/${PROJECT}_params_ampliseq.yaml" \
  --input "input/manifest/${PROJECT}_ampliseq.tsv" \
  --outdir "outputs/ampliseq/${PROJECT}"



