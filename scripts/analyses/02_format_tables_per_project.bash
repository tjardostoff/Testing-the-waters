#!/bin/bash

#SBATCH --partition long
#SBATCH --cpus-per-task 1
#SBATCH --mem 5GB

while read PROJECT; do

  sbatch \
    --wait \
    --export=PROJECT=${PROJECT} \
    --output 02_subset_and_format_table_${PROJECT}.out \
    scripts/bash/subset_and_format_table.bash

done < tmp/my_projects.txt


