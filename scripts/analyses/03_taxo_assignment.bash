#!/bin/bash

#SBATCH --partition long
#SBATCH --cpus-per-task 1
#SBATCH --mem 5GB
#SBATCH --output 03_taxo_assignment.out

zcat outputs/per_project/*/dada2_asvs.tsv.gz outputs/per_project/*/vsearch_otus.tsv.gz | \
    sort --uniq | \
    awk -F "\t" '{print ">"$1"\n"$2}' | \
    gzip > outputs/taxo_assignment/sequences_to_assign.fasta.gz

while read REFDB; do

  echo $REFDB" is running"

  sbatch \
    --wait \
    --export=REFDB=${REFDB} \
    --output 03_taxo_assignment_${REFDB}.out \
    scripts/bash/taxo_assign.bash

  wait

  echo "done"

done < tmp/refdbs.txt

echo "All done"
