library(tidyverse)
library(readxl)

####################################################
# Input
####################################################

args <- commandArgs(trailingOnly = TRUE)
print(args)
metaproject <- args[1]
favourite_assign <- args[2]
refdbs <- readLines("tmp/refdbs.txt")

counts_asvs_files <-
  dir(
    "outputs/per_project",
    recursive = TRUE,
    full.names = TRUE,
    pattern = "dada2_counts\\.tsv\\.gz$"
  )

counts_asvs <- map(counts_asvs_files, read_tsv) |> bind_rows()

asvs_in_counts <- counts_asvs |>
  pull(asv_id) |>
  unique()

counts_otus_files <-
  dir(
    "outputs/per_project",
    recursive = TRUE,
    full.names = TRUE,
    pattern = "vsearch_counts\\.tsv\\.gz$"
  )

counts_otus <- map(counts_otus_files, read_tsv) |> bind_rows()

otus_in_counts <- counts_otus |>
  pull(asv_id) |>
  unique()

sequences <- Biostrings::readDNAStringSet("outputs/taxo_assignment/sequences_to_assign.fasta.gz")

sequences <- tibble(
  asv_id = names(sequences),
  sequence = unname(as.character(sequences))
)

####################################################
# Taxonomy
####################################################

taxo_res_importer <- function(my_refdb) {

  taxo_res <- list()

  taxo_dir <- "outputs/taxo_assignment/"

  dada2_res_path <- str_c(taxo_dir, my_refdb, "/dada2_res.tsv")
  idataxa_res_path <- str_c(taxo_dir, my_refdb, "/idtaxa_res.tsv")
  vsearch_res_path <- str_c(taxo_dir, my_refdb, "/vsearch_best_hit.tsv")

  if (file.exists(dada2_res_path)) {

    taxo_res[["dada2"]] <-
      read_tsv(
        file = dada2_res_path
      ) |>
      filter(asv_id %in% c(asvs_in_counts, otus_in_counts))

  }

  if (file.exists(idataxa_res_path)) {

    taxo_res[["idtaxa"]] <-
      read_tsv(
        file = idataxa_res_path
      ) |>
      filter(asv_id %in% c(asvs_in_counts, otus_in_counts))

  }

  if (file.exists(vsearch_res_path)) {

    taxo_res[["vsearch"]] <-
      read_tsv(
        file = vsearch_res_path,
        col_select = -lower_threshold_score
      ) |>
      rename(
        similarity = best_hit_score,
        taxonomy = lca_taxonomy
      ) |>
      relocate(taxonomy, .after = asv_id) |>
      right_join(tibble(asv_id = unique(c(asvs_in_counts, otus_in_counts))), by = "asv_id")

  }

  names(taxo_res) <- str_c(my_refdb, names(taxo_res), sep = "_")

  return(taxo_res)

}

taxo_res <- map(refdbs, taxo_res_importer) |>
  reduce(c)

renaming_func <- function(tax_table, tax_name) {
  rename_with(
    tax_table,
    ~ str_c(tax_name, "_", .x),
    .cols = -asv_id
  )
}

taxo_res <-
  map(
    names(taxo_res),
    ~ taxo_res[[.x]] |> renaming_func(.x)
  ) |>
  reduce(full_join, by = "asv_id")

write_tsv(
  taxo_res,
  str_c("outputs/final_tables/",
  metaproject,
  "_dada2_taxo_extra.tsv.gz")
)

taxo_in_asv_desc <- taxo_res |>
  select(asv_id, matches(str_c("_", favourite_assign, "_"))) |>
  filter(asv_id %in% asvs_in_counts) |>
  left_join(sequences, by = "asv_id") |>
  relocate(sequence, .after = "asv_id")

taxo_in_otu_desc <- taxo_res |>
  select(asv_id, matches(str_c("_", favourite_assign, "_"))) |>
  filter(asv_id %in% otus_in_counts) |>
  left_join(sequences, by = "asv_id") |>
  relocate(sequence, .after = "asv_id")

#####################################################################
# counts
#####################################################################

write_tsv(
  counts_asvs,
  str_c("outputs/final_tables/", metaproject, "_dada2_asvs_counts.tsv.gz")
)

write_tsv(
  counts_asvs,
  str_c("outputs/final_tables/", metaproject, "_dada2_otus_counts.tsv.gz")
)

#####################################################################
# overall summary
#####################################################################

overall_summary_files <- dir(
  "outputs/ampliseq",
  recursive = TRUE,
  full.names = TRUE,
  pattern = "overall_summary\\.tsv$"
)

overall_summary <- map(overall_summary_files, read_tsv) |> bind_rows()

# sum orientations

overall_summary <- overall_summary |> 
  mutate(sample = str_replace(sample,"^([^_]+_[^_]+)_.+$", "\\1")) |> 
  group_by(sample) |> 
  summarise(
    cutadapt_total_processed = unique(cutadapt_total_processed),
    cutadapt_passing_filters = sum(cutadapt_passing_filters),
    cutadapt_passing_filters_percent = (cutadapt_passing_filters / cutadapt_total_processed) * 100,
    cutadapt_passing_filters_percent = paste0(round(cutadapt_passing_filters_percent, digits = 1),"%"),
    DADA2_input = sum(DADA2_input),
    filtered = sum(filtered),
    denoisedF = sum(denoisedF),
    denoisedR = sum(denoisedR),
    merged = sum(merged),
    nonchim = sum(nonchim)
  ) |>
  ungroup() |>
  mutate(sample = str_remove(sample, "_.+")) |>
  group_by(sample) |>
  summarise(
    n_replicates = n(),
    cutadapt_total_processed = sum(cutadapt_total_processed),
    cutadapt_passing_filters = sum(cutadapt_passing_filters),
    cutadapt_passing_filters_percent = (cutadapt_passing_filters / cutadapt_total_processed) * 100,
    cutadapt_passing_filters_percent = paste0(round(cutadapt_passing_filters_percent, digits = 1),"%"),
    DADA2_input = sum(DADA2_input),
    filtered = sum(filtered),
    denoisedF = sum(denoisedF),
    denoisedR = sum(denoisedR),
    merged = sum(merged),
    nonchim = sum(nonchim)
  )

summary_otus <- counts_otus |>
  group_by(sample) |>
  summarise(
    final_nreads_otus = sum(nreads),
    final_notus = n_distinct(asv_id)
  )

overall_summary <-
  counts_asvs |>
  group_by(sample) |>
  summarise(
    final_nreads = sum(nreads),
    final_nasvs = n_distinct(asv_id)
  ) |>
  inner_join(summary_otus, ., by = "sample") |>
  inner_join(overall_summary, ., by = "sample")

write_tsv(overall_summary, str_c("outputs/final_tables/", metaproject, "_dada2_overall_summary.tsv.gz"))

####################################################
# ASVs statistics + selected taxo
####################################################

asvs_desc <- counts_asvs |>
  group_by(asv_id) |>
  summarise(total = sum(nreads), spread = n()) |>
  ungroup() |>
  left_join(taxo_in_asv_desc, by = "asv_id") |>
  mutate(
    sequence_length = nchar(sequence),
  ) |>
  arrange(desc(total))


write_tsv(asvs_desc, str_c("outputs/final_tables/", metaproject, "_dada2_asvs.tsv.gz"))


sequences_to_export_asvs <- setNames(asvs_desc$sequence, asvs_desc$asv_id)
sequences_to_export_asvs <- Biostrings::DNAStringSet(sequences_to_export_asvs, use.names = TRUE)

Biostrings::writeXStringSet(
  sequences_to_export_asvs,
  str_c("outputs/final_tables/", metaproject, "_dada2_asvs.fasta.gz"),
  compress = TRUE,
  format = "fasta"
)

####################################################
# OTUs statistics + selected taxo
####################################################

otus_desc <- counts_otus |>
  group_by(asv_id) |>
  summarise(total = sum(nreads), spread = n()) |>
  ungroup() |>
  left_join(taxo_in_otu_desc, by = "asv_id") |>
  mutate(
    sequence_length = nchar(sequence),
  ) |>
  arrange(desc(total))


write_tsv(otus_desc, str_c("outputs/final_tables/", metaproject, "_dada2_otus.tsv.gz"))


sequences_to_export_otus <- setNames(otus_desc$sequence, otus_desc$asv_id)
sequences_to_export_otus <- Biostrings::DNAStringSet(sequences_to_export_otus, use.names = TRUE)

Biostrings::writeXStringSet(
  sequences_to_export_otus,
  str_c("outputs/final_tables/", metaproject, "_dada2_otus.fasta.gz"),
  compress = TRUE,
  format = "fasta"
)

####################################################
# nextflow summary report
####################################################

projects <- list.dirs(
  "outputs/ampliseq",
  recursive = FALSE,
  full.names = FALSE
)

copy_report <- function(my_project) {

  file.copy(
    from = str_c("outputs/ampliseq/", my_project, "/summary_report/summary_report.html"),
    to = str_c("outputs/final_tables/", my_project, "_report.html"),
    overwrite = TRUE
  )

}

map(projects, copy_report)

####################################################
# Context
####################################################

context <- read_xlsx("archive/context/SpatialCampaign_phyloseq.xlsx")

context <- context |>
  mutate(
    date = as_date(date),
    Seq_ID = str_remove(Seq_ID, "SC-"),
    Seq_ID = str_replace(Seq_ID, "ED0?", "EDNA"),
    Sample_Plot = str_c(Sample_Type, str_pad(map_Id, 2,pad = "0", side = "left"), Replicate_ID, sep = "_")
  ) |>
  relocate(Sample_Plot, .after = Seq_ID)

write_tsv(context, str_c("outputs/final_tables/", metaproject, "_context.tsv"))


####################################################
# phyloseq export
####################################################

context <- context |>
  column_to_rownames("Seq_ID")

export_phyloseq <- function(counts, asvs_or_otus, refdb) {

  taxo_column <- str_c(refdb, "_", favourite_assign,"_taxonomy")
  sim_column <- str_c(refdb, "_", favourite_assign,"_similarity")

  counts_for_phyloseq <- 
    counts |>
    pivot_wider(names_from = "sample", values_from = "nreads", values_fill = 0) |> 
    column_to_rownames("asv_id")

  taxo_phyloseq <- taxo_res %>%
    select(asv_id, all_of(c(taxo_column, sim_column))) %>%
    separate_wider_delim(all_of(taxo_column), delim = ";", names_sep = "_rank", too_few = "align_start") |>
    column_to_rownames("asv_id")

  phyloseq_data <- phyloseq::phyloseq(
    phyloseq::otu_table(counts_for_phyloseq, taxa_are_rows = TRUE),
    phyloseq::tax_table(as.matrix(taxo_phyloseq)),
    phyloseq::sample_data(context),
    phyloseq::refseq(switch(asvs_or_otus, asvs = sequences_to_export_asvs, otus = sequences_to_export_otus))
  )

  saveRDS(
    phyloseq_data,
    file = file.path("outputs/final_tables", str_c(metaproject, "dada2",asvs_or_otus, refdb, "phyloseq.rds", sep = "_"))
  )

}

export_phyloseq(counts = counts_asvs, asvs_or_otus = "asvs", refdb = "ekoi")
export_phyloseq(counts = counts_asvs, asvs_or_otus = "asvs", refdb = "bold")
export_phyloseq(counts = counts_asvs, asvs_or_otus = "asvs", refdb = "mzgdb")

export_phyloseq(counts = counts_otus, asvs_or_otus = "otus", refdb = "ekoi")
export_phyloseq(counts = counts_otus, asvs_or_otus = "otus", refdb = "bold")
export_phyloseq(counts = counts_otus, asvs_or_otus = "otus", refdb = "mzgdb")

####################################################
# session info
####################################################

sessionInfo()