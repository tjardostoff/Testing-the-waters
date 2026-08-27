#!/bin/bash

#SBATCH --partition long
#SBATCH --cpus-per-task 1
#SBATCH --mem 20GB

if [ ! -d "outputs/per_project/${PROJECT}" ]; then

  mkdir -p outputs/per_project/${PROJECT}

  ## ASV table

  # Input files
  ASV_TABLE="outputs/ampliseq/${PROJECT}/dada2/ASV_table.tsv"

  # Output files
  COUNTS_TABLE="outputs/per_project/${PROJECT}/dada2_counts.tsv.gz"
  ASV_DESC="outputs/per_project/${PROJECT}/dada2_asvs.tsv.gz"

  # TMP file
  TMPFILE=$(mktemp -p tmp/)

  # Long format and subset
  head -n 1 $ASV_TABLE | \
    tr '\t' '\n' | \
    awk 'NR >1 {printf $1"\t"; sub(/_([0-9]+_)?orientation_[12]/,""); print}' > $TMPFILE

  awk \
    -v mintotal=3 \
    -v minspread=2 \
    -v asv_id_index=1 \
    -f scripts/awk/asv_table_longer.awk \
    $TMPFILE \
    $ASV_TABLE | \
    gzip > $COUNTS_TABLE

  # Export sequences as table
  cat outputs/ampliseq/${PROJECT}/dada2/ASV_seqs.fasta | \
      paste - - | \
      sed 's/>//' | \
      sort | \
      join - <(zcat $COUNTS_TABLE | awk -F "\t" '{print $1}' | sort --uniq) -t $'\t' | \
      gzip > $ASV_DESC

  rm $TMPFILE

  ## OTU table

  # Input files
  OTU_TABLE="outputs/ampliseq/${PROJECT}/vsearch_cluster/ASV_post_clustering_filtered.table.tsv"

  # Output files
  COUNTS_TABLE="outputs/per_project/${PROJECT}/vsearch_counts.tsv.gz"
  OTU_DESC="outputs/per_project/${PROJECT}/vsearch_otus.tsv.gz"

  # TMP file
  TMPFILE=$(mktemp -p tmp/)

  # Long format and subset
  head -n 1 $OTU_TABLE | \
    tr '\t' '\n' | \
    awk 'NR >1 {printf $1"\t"; sub(/_([0-9]+_)?orientation_[12]/,""); print}' > $TMPFILE

  awk \
    -v mintotal=3 \
    -v minspread=2 \
    -v asv_id_index=1 \
    -f scripts/awk/asv_table_longer.awk \
    $TMPFILE \
    $OTU_TABLE | \
    gzip > $COUNTS_TABLE

  # Export sequences as table
  cat outputs/ampliseq/${PROJECT}/vsearch_cluster/ASV_post_clustering_filtered.fna | \
      awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' | \
      paste - - | \
      sed 's/>//' | \
      sort | \
      join - <(zcat $COUNTS_TABLE | awk -F "\t" '{print $1}' | sort --uniq) -t $'\t' | \
      gzip > $OTU_DESC

  rm $TMPFILE

fi




     
