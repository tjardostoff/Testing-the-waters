#!/bin/bash

#SBATCH --partition long
#SBATCH --cpus-per-task 1
#SBATCH --mem 5GB

while read PROJECT; do

  echo $PROJECT"is running"

  [ ! -d "outputs/ampliseq/${PROJECT}" ] && sbatch \
    --wait \
    --export=PROJECT=${PROJECT} \
    --output 01_ampliseq_${PROJECT}.out \
    scripts/bash/ampliseq.bash

  wait

  echo "done"

done < tmp/my_projects.txt


