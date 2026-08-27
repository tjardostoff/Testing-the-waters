library(phyloseq)
library(tidyverse)

###repeat until species table for each assignment database###
####phyloseq ASV ekoi only metazoans####

#load the phyloseq object
ekoi_asv <- readRDS("outputs/final_tables/coastclim_spatial_dada2_asvs_ekoi_phyloseq.rds")
ekoi_asv<- subset_taxa(ekoi_asv, ekoi_vsearch_taxonomy_rank4 =="Metazoa")
ekoi_asv@otu_table # to extract the ASV table
ekoi_asv@tax_table # to extract the taxonomy table
ekoi_asv@sam_data # to extract the sample table

#overview commands
sample_sums(ekoi_asv)
rank_names(ekoi_asv)

get_taxa_unique(ekoi_asv, taxonomic.rank = "ekoi_vsearch_taxonomy_rank4")

#filter data to remove blanks
ekoi_asv_small <-
  ekoi_asv |> 
  subset_samples(Sample_Type %in% c("BULK","EDNA")) |>
  filter_taxa(function(x) sum(x) > 0, prune = TRUE)

plot_bar(ekoi_asv_small)

#transform to relative number of reads
ekoi_asv_perc <- phyloseq::transform_sample_counts(
  ekoi_asv_small,
  function(x) x/sum(x) * 100
)

plot_bar(ekoi_asv_perc)

#change taxrank for different taxonomic level

ekoi_asv_perc_division <- phyloseq::tax_glom(
  ekoi_asv_perc,
  taxrank = "ekoi_vsearch_taxonomy_rank10",
  NArm = FALSE
)

plot_bar(
  ekoi_asv_perc_division,
  fill = "ekoi_vsearch_taxonomy_rank10"
) 

#remove rare ranks
ekoi_asv_perc_division_mean1 <- phyloseq::filter_taxa(
  ekoi_asv_perc_division,
  flist = function(x) mean(x) >= 1,
  prune = TRUE
)

#Plot at Division taxonomic rank, with colour
barplot_ekoi_asv <- plot_bar(ekoi_asv_perc_division_mean1, fill = "ekoi_vsearch_taxonomy_rank10", x="Sample_Plot") 
ggsave("barplot_ekoi.pdf")

#subset only crustaceans
cop <- subset_taxa(ekoi_asv_perc, ekoi_vsearch_taxonomy_rank6 =="Copepoda")
#heatmap with only copepods
plot_heatmap(cop, sample.label = "Sample_Plot", sample.order="Sample_Plot", taxa.label = "ekoi_vsearch_taxonomy_rank10")
ggsave("test.pdf", width = 6, height = 20, units = "in" )

##get species table
taxo <- read_tsv("outputs/final_tables/coastclim_spatial_dada2_taxo_extra.tsv.gz")
ekoi_asv_counts <- read_tsv("outputs/final_tables/coastclim_spatial_dada2_asvs_counts.tsv.gz")

# for each species/sample the number of OTUs detected.
# NA means no OTU detected for that species

ekoi_asv_species <- taxo |> 
  select(asv_id, taxonomy = ekoi_vsearch_taxonomy) |>
  mutate(taxonomy = str_remove(taxonomy, ".+;")) |> 
  filter(!is.na(taxonomy) & taxonomy != "*") |> 
  inner_join(ekoi_asv_counts, by = "asv_id") |> 
  pivot_wider(names_from = sample, id_cols = taxonomy, values_from = asv_id, values_fn = n_distinct)

write.csv(ekoi_asv_species, "outputs/final_tables/species_asv_ekoi.csv")

####phyloseq ASV bold####

#load the phyloseq object
bold_asv <- readRDS("outputs/final_tables/coastclim_spatial_dada2_asvs_bold_phyloseq.rds")

bold_asv@otu_table # to extract the ASV table
bold_asv@tax_table # to extract the taxonomy table
bold_asv@sam_data # to extract the sample table

#overview commands
sample_sums(bold_asv)
rank_names(bold_asv)

get_taxa_unique(bold_asv, taxonomic.rank = "bold_vsearch_taxonomy_rank1")

#filter data to remove blanks
bold_asv_small <-
  bold_asv |> 
  subset_samples(Sample_Type %in% c("BULK","EDNA")) |>
  filter_taxa(function(x) sum(x) > 0, prune = TRUE)

plot_bar(bold_asv_small)

#transform to relative number of reads
bold_asv_perc <- phyloseq::transform_sample_counts(
  bold_asv_small,
  function(x) x/sum(x) * 100
)

plot_bar(bold_asv_perc)

#change taxrank for different taxonomic level

bold_asv_perc_division <- phyloseq::tax_glom(
  bold_asv_perc,
  taxrank = "bold_vsearch_taxonomy_rank1",
  NArm = FALSE
)

plot_bar(
  bold_asv_perc_division,
  fill = "bold_vsearch_taxonomy_rank1"
) 

#remove rare ranks
bold_asv_perc_division_mean1 <- phyloseq::filter_taxa(
  bold_asv_perc_division,
  flist = function(x) mean(x) >= 1,
  prune = TRUE
)

#Plot at Division taxonomic rank, with colour
barplot_bold_asv <- plot_bar(bold_asv_perc_division_mean1, fill = "bold_vsearch_taxonomy_rank1", x="Sample_Plot", legend=FALSE) 
ggsave("barplot_bold.pdf")

#subset only crustaceans
cop <- subset_taxa(bold_asv_perc, bold_vsearch_taxonomy_rank6 =="Copepoda")
#heatmap with only copepods
plot_heatmap(cop, sample.label = "Sample_Plot", sample.order="Sample_Plot", taxa.label = "bold_vsearch_taxonomy_rank10")
ggsave("test.pdf", width = 6, height = 20, units = "in" )

##get species table (change counts table for otus)
taxo <- read_tsv("outputs/final_tables/coastclim_spatial_dada2_taxo_extra.tsv.gz")
asv_counts <- read_tsv("outputs/final_tables/coastclim_spatial_dada2_asvs_counts.tsv.gz")

# for each species/sample the number of OTUs detected.
# NA means no OTU detected for that species

bold_asv_species <- taxo |> 
  select(asv_id, taxonomy = bold_vsearch_taxonomy) |>
  mutate(taxonomy = str_remove(taxonomy, ".+;")) |> 
  filter(!is.na(taxonomy) & taxonomy != "*") |> 
  inner_join(asv_counts, by = "asv_id") |> 
  pivot_wider(names_from = sample, id_cols = taxonomy, values_from = asv_id, values_fn = n_distinct)

write.csv(bold_asv_species, "outputs/final_tables/species_asv_bold.csv")

####phyloseq ASV mzgdb####

#load the phyloseq object
mzgdb_asv <- readRDS("outputs/final_tables/coastclim_spatial_dada2_asvs_mzgdb_phyloseq.rds")

mzgdb_asv@otu_table # to extract the ASV table
mzgdb_asv@tax_table # to extract the taxonomy table
mzgdb_asv@sam_data # to extract the sample table

#overview commands
sample_sums(mzgdb_asv)
rank_names(mzgdb_asv)

get_taxa_unique(mzgdb_asv, taxonomic.rank = "mzgdb_vsearch_taxonomy_rank20")

#filter data to remove blanks
mzgdb_asv_small <-
  mzgdb_asv |> 
  subset_samples(Sample_Type %in% c("BULK","EDNA")) |>
  filter_taxa(function(x) sum(x) > 0, prune = TRUE)

plot_bar(mzgdb_asv_small)

#transform to relative number of reads
mzgdb_asv_perc <- phyloseq::transform_sample_counts(
  mzgdb_asv_small,
  function(x) x/sum(x) * 100
)

plot_bar(mzgdb_asv_perc)

#change taxrank for different taxonomic level

mzgdb_asv_perc_division <- phyloseq::tax_glom(
  mzgdb_asv_perc,
  taxrank = "mzgdb_vsearch_taxonomy_rank20",
  NArm = FALSE
)

plot_bar(
  mzgdb_asv_perc_division,
  fill = "mzgdb_vsearch_taxonomy_rank20"
) 

#remove rare ranks
mzgdb_asv_perc_division_mean1 <- phyloseq::filter_taxa(
  mzgdb_asv_perc_division,
  flist = function(x) mean(x) >= 1,
  prune = TRUE
)

#Plot at Division taxonomic rank, with colour
barplot_mzgdb_asv <- plot_bar(mzgdb_asv_perc_division_mean1, fill = "mzgdb_vsearch_taxonomy_rank20", x="Sample_Plot") 
ggsave("barplot_ekoi.pdf")

#subset only crustaceans
cop <- subset_taxa(mzgdb_asv_perc, mzgdb_vsearch_taxonomy_rank6 =="Copepoda")
#heatmap with only copepods
plot_heatmap(cop, sample.label = "Sample_Plot", sample.order="Sample_Plot", taxa.label = "mzgdb_vsearch_taxonomy_rank10")
ggsave("test.pdf", width = 6, height = 20, units = "in" )

##get species table (change counts table for otus)
taxo <- read_tsv("outputs/final_tables/coastclim_spatial_dada2_taxo_extra.tsv.gz")
asv_counts <- read_tsv("outputs/final_tables/coastclim_spatial_dada2_asvs_counts.tsv.gz")

# for each species/sample the number of OTUs detected.
# NA means no OTU detected for that species

mzgdb_asv_species <- taxo |> 
  select(asv_id, taxonomy = mzgdb_vsearch_taxonomy) |>
  mutate(taxonomy = str_remove(taxonomy, ".+;")) |> 
  filter(!is.na(taxonomy) & taxonomy != "*") |> 
  inner_join(asv_counts, by = "asv_id") |> 
  pivot_wider(names_from = sample, id_cols = taxonomy, values_from = asv_id, values_fn = n_distinct)

write.csv(mzgdb_asv_species, "outputs/final_tables/species_asv_mzgdb.csv")

###repeat until species table for each assignment database###
####phyloseq OTU ekoi####

#load the phyloseq object
ekoi_otu <- readRDS("outputs/final_tables/coastclim_spatial_dada2_otus_ekoi_phyloseq.rds")
ekoi_otu <- subset_taxa(ekoi_otu, ekoi_vsearch_taxonomy_rank4 =="Metazoa")
ekoi_otu@otu_table # to extract the ASV table
ekoi_otu@tax_table # to extract the taxonomy table
ekoi_otu@sam_data # to extract the sample table

#overview commands
sample_sums(ekoi_otu)
rank_names(ekoi_otu)

get_taxa_unique(ekoi_otu, taxonomic.rank = "ekoi_vsearch_taxonomy_rank6")

#filter data to remove blanks
ekoi_otu_small <-
  ekoi_otu |> 
  subset_samples(Sample_Type %in% c("BULK","EDNA")) |>
  filter_taxa(function(x) sum(x) > 0, prune = TRUE)

plot_bar(ekoi_otu_small)

#transform to relative number of reads
ekoi_otu_perc <- phyloseq::transform_sample_counts(
  ekoi_otu_small,
  function(x) x/sum(x) * 100
)

plot_bar(ekoi_otu_perc)

#change taxrank for different taxonomic level

ekoi_otu_perc_division <- phyloseq::tax_glom(
  ekoi_otu_perc,
  taxrank = "ekoi_vsearch_taxonomy_rank10",
  NArm = FALSE
)

plot_bar(
  ekoi_otu_perc_division,
  fill = "ekoi_vsearch_taxonomy_rank10"
) 

#remove rare ranks
ekoi_otu_perc_division_mean1 <- phyloseq::filter_taxa(
  ekoi_otu_perc_division,
  flist = function(x) mean(x) >= 1,
  prune = TRUE
)

#Plot at Division taxonomic rank, with colour
barplot_ekoi_otu <- plot_bar(ekoi_otu_perc_division_mean1, fill = "ekoi_vsearch_taxonomy_rank10", x="Sample_Plot") 


#subset only crustaceans
cop <- subset_taxa(ekoi_otu_perc, ekoi_vsearch_taxonomy_rank6 =="Copepoda")
#heatmap with only copepods
plot_heatmap(cop, sample.label = "Sample_Plot", sample.order="Sample_Plot", taxa.label = "ekoi_vsearch_taxonomy_rank10")
ggsave("test.pdf", width = 6, height = 20, units = "in" )

##get species table
taxo <- read_tsv("outputs/final_tables/coastclim_spatial_dada2_taxo_extra.tsv.gz")
otu_counts <- read_tsv("outputs/final_tables/coastclim_spatial_dada2_otus_counts.tsv.gz")

# for each species/sample the number of OTUs detected.
# NA means no OTU detected for that species

ekoi_otu_species <- taxo |> 
  select(asv_id, taxonomy = ekoi_vsearch_taxonomy) |>
  mutate(taxonomy = str_remove(taxonomy, ".+;")) |> 
  filter(!is.na(taxonomy) & taxonomy != "*") |> 
  inner_join(otu_counts, by = "asv_id") |> 
  pivot_wider(names_from = sample, id_cols = taxonomy, values_from = asv_id, values_fn = n_distinct)

write.csv(ekoi_otu_species, "outputs/final_tables/species_otu_ekoi.csv")

####phyloseq OTU bold####

#load the phyloseq object
bold_otu <- readRDS("outputs/final_tables/coastclim_spatial_dada2_otus_bold_phyloseq.rds")

bold_otu@otu_table # to extract the OTU table
bold_otu@tax_table # to extract the taxonomy table
bold_otu@sam_data # to extract the sample table

#overview commands
sample_sums(bold_otu)
rank_names(bold_otu)

get_taxa_unique(bold_otu, taxonomic.rank = "bold_vsearch_taxonomy_rank1")

#filter data to remove blanks
bold_otu_small <-
  bold_otu |> 
  subset_samples(Sample_Type %in% c("BULK","EDNA")) |>
  filter_taxa(function(x) sum(x) > 0, prune = TRUE)

plot_bar(bold_otu_small)

#transform to relative number of reads
bold_otu_perc <- phyloseq::transform_sample_counts(
  bold_otu_small,
  function(x) x/sum(x) * 100
)

plot_bar(bold_otu_perc)

#change taxrank for different taxonomic level

bold_otu_perc_division <- phyloseq::tax_glom(
  bold_otu_perc,
  taxrank = "bold_vsearch_taxonomy_rank1",
  NArm = FALSE
)

plot_bar(
  bold_otu_perc_division,
  fill = "bold_vsearch_taxonomy_rank1"
) 

#remove rare ranks
bold_otu_perc_division_mean1 <- phyloseq::filter_taxa(
  bold_otu_perc_division,
  flist = function(x) mean(x) >= 1,
  prune = TRUE
)

#Plot at Division taxonomic rank, with colour
barplot_bold_otu <- plot_bar(bold_otu_perc_division_mean1, fill = "bold_vsearch_taxonomy_rank1", x="Sample_Plot") 


#subset only crustaceans
cop <- subset_taxa(bold_otu_perc, bold_vsearch_taxonomy_rank6 =="Copepoda")
#heatmap with only copepods
plot_heatmap(cop, sample.label = "Sample_Plot", sample.order="Sample_Plot", taxa.label = "bold_vsearch_taxonomy_rank10")
ggsave("test.pdf", width = 6, height = 20, units = "in" )

##get species table (change counts table for otus)
taxo <- read_tsv("outputs/final_tables/coastclim_spatial_dada2_taxo_extra.tsv.gz")
otu_counts <- read_tsv("outputs/final_tables/coastclim_spatial_dada2_otus_counts.tsv.gz")

# for each species/sample the number of OTUs detected.
# NA means no OTU detected for that species

bold_otu_species <- taxo |> 
  select(asv_id, taxonomy = bold_vsearch_taxonomy) |>
  mutate(taxonomy = str_remove(taxonomy, ".+;")) |> 
  filter(!is.na(taxonomy) & taxonomy != "*") |> 
  inner_join(asv_counts, by = "asv_id") |> 
  pivot_wider(names_from = sample, id_cols = taxonomy, values_from = asv_id, values_fn = n_distinct)

write.csv(bold_otu_species, "outputs/final_tables/species_otu_bold.csv")

####phyloseq OTU mzgdb####

#load the phyloseq object
mzgdb_otu <- readRDS("outputs/final_tables/coastclim_spatial_dada2_otus_mzgdb_phyloseq.rds")

mzgdb_otu@otu_table # to extract the OTU table
mzgdb_otu@tax_table # to extract the taxonomy table
mzgdb_otu@sam_data # to extract the sample table

#overview commands
sample_sums(mzgdb_otu)
rank_names(mzgdb_otu)

get_taxa_unique(mzgdb_otu, taxonomic.rank = "mzgdb_vsearch_taxonomy_rank20")

#filter data to remove blanks
mzgdb_otu_small <-
  mzgdb_otu |> 
  subset_samples(Sample_Type %in% c("BULK","EDNA")) |>
  filter_taxa(function(x) sum(x) > 0, prune = TRUE)

plot_bar(mzgdb_otu_small)

#transform to relative number of reads
mzgdb_otu_perc <- phyloseq::transform_sample_counts(
  mzgdb_otu_small,
  function(x) x/sum(x) * 100
)

plot_bar(mzgdb_otu_perc)

#change taxrank for different taxonomic level

mzgdb_otu_perc_division <- phyloseq::tax_glom(
  mzgdb_otu_perc,
  taxrank = "mzgdb_vsearch_taxonomy_rank20",
  NArm = FALSE
)

plot_bar(
  mzgdb_otu_perc_division,
  fill = "mzgdb_vsearch_taxonomy_rank20"
) 

#remove rare ranks
mzgdb_otu_perc_division_mean1 <- phyloseq::filter_taxa(
  mzgdb_otu_perc_division,
  flist = function(x) mean(x) >= 1,
  prune = TRUE
)

#Plot at Division taxonomic rank, with colour
barplot_mzgdb_otu <- plot_bar(mzgdb_otu_perc_division_mean1, fill = "mzgdb_vsearch_taxonomy_rank20", x="Sample_Plot") 


#subset only crustaceans
cop <- subset_taxa(mzgdb_otu_perc, mzgdb_vsearch_taxonomy_rank6 =="Copepoda")
#heatmap with only copepods
plot_heatmap(cop, sample.label = "Sample_Plot", sample.order="Sample_Plot", taxa.label = "mzgdb_vsearch_taxonomy_rank10")
ggsave("test.pdf", width = 6, height = 20, units = "in" )

##get species table (change counts table for otus)
taxo <- read_tsv("outputs/final_tables/coastclim_spatial_dada2_taxo_extra.tsv.gz")
otu_counts <- read_tsv("outputs/final_tables/coastclim_spatial_dada2_otus_counts.tsv.gz")

# for each species/sample the number of OTUs detected.
# NA means no OTU detected for that species

mzgdb_otu_species <- taxo |> 
  select(asv_id, taxonomy = mzgdb_vsearch_taxonomy) |>
  mutate(taxonomy = str_remove(taxonomy, ".+;")) |> 
  filter(!is.na(taxonomy) & taxonomy != "*") |> 
  inner_join(asv_counts, by = "asv_id") |> 
  pivot_wider(names_from = sample, id_cols = taxonomy, values_from = asv_id, values_fn = n_distinct)

write.csv(mzgdb_otu_species, "outputs/final_tables/species_otu_mzgdb.csv")

library(cowplot)

plot_grid(barplot_bold_asv, barplot_ekoi_asv, barplot_mzgdb_asv)

plot_grid(barplot_bold_otu, barplot_ekoi_otu, barplot_mzgdb_otu)

#rarefaction for ASVs (doesn't matter if ekoi, bold or mzgdb)
get_raref_plot <- function(ekoi_asv, repeats, threads) {

  require(rtk)
  
  asv_table <- ekoi_asv@otu_table |>
    data.frame()

  get_raref_plot_values <- function(my_sample) {

    table_subset <- asv_table[my_sample] 

    read_depths <- c(1,((seq(0.1,1,0.05))^8) * sum(table_subset[[1]])) |>
      ceiling() |>
      unique()

    rtk_output <- rtk(
      table_subset,
      depth = read_depths,
      repeats = repeats,
      threads = threads
    )

    tibble(
      sample = my_sample,
      nreads = read_depths,
      nasvs = get.median.diversity(rtk_output) |> unlist()
    )

  }

  samples <- names(asv_table)

  raref_plot_values <- map(samples, get_raref_plot_values) |>
    bind_rows()

  raref_plot_values <- inner_join(
    raref_plot_values,
    phyloseq_object@sam_data |> data.frame() |> rownames_to_column(var = "sample"),
    by = "sample"
  )
  
  ggplot(raref_plot_values) +
    aes(x = nreads, y = nasvs, group = sample, colour = sample) +
    geom_line() +
    labs(x = "# reads", y = "# ASVs")

}

raref_plot <- get_raref_plot(ekoi_asv_small, repeats = 3, threads = 2)

raref_plot

raref_plot +
  facet_grid(map_Id~Sample_Type) +
  scale_x_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  theme(legend.position = "none")


plot_richness(phyloseq_data_rar, measures = c("Observed", "Shannon"))

plot_richness(phyloseq_data_rar, x = "map_Id", measures = "Shannon")

##prepare nmds
phyloseq_data_BC <- phyloseq::distance(ekoi_asv_perc, method = "bray", type = "samples")

ekoi_asv_nmds <- phyloseq::ordinate(
  ekoi_asv_perc,
  "NMDS",
  distance = phyloseq_data_BC
)

ekoi_asv_nmds<-phyloseq::plot_ordination(
  ekoi_asv_perc, 
  ekoi_asv_nmds, 
  type = "samples",
  label = "map_Id",
  color = "Sample_Type"
)
plot(ekoi_asv_nmds)

##dendogram
phyloseq_data_BC_hclust_average <- hclust(phyloseq_data_BC, method = "average")

plot(phyloseq_data_BC_hclust_average)


#rarefaction for OTUs (doesn't matter if ekoi, bold or mzgdb)
get_raref_plot <- function(ekoi_otu, repeats, threads) {
  
  require(rtk)
  
  asv_table <- ekoi_otu@otu_table |>
    data.frame()
  
  get_raref_plot_values <- function(my_sample) {
    
    table_subset <- asv_table[my_sample] 
    
    read_depths <- c(1,((seq(0.1,1,0.05))^8) * sum(table_subset[[1]])) |>
      ceiling() |>
      unique()
    
    rtk_output <- rtk(
      table_subset,
      depth = read_depths,
      repeats = repeats,
      threads = threads
    )
    
    tibble(
      sample = my_sample,
      nreads = read_depths,
      nasvs = get.median.diversity(rtk_output) |> unlist()
    )
    
  }
  
  samples <- names(asv_table)
  
  raref_plot_values <- map(samples, get_raref_plot_values) |>
    bind_rows()
  
  raref_plot_values <- inner_join(
    raref_plot_values,
    phyloseq_object@sam_data |> data.frame() |> rownames_to_column(var = "sample"),
    by = "sample"
  )
  
  ggplot(raref_plot_values) +
    aes(x = nreads, y = nasvs, group = sample, colour = sample) +
    geom_line() +
    labs(x = "# reads", y = "# ASVs")
  
}

raref_plot <- get_raref_plot(ekoi_asv_small, repeats = 3, threads = 2)

raref_plot

raref_plot +
  facet_grid(map_Id~Sample_Type) +
  scale_x_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  theme(legend.position = "none")


plot_richness(phyloseq_data_rar, measures = c("Observed", "Shannon"))

plot_richness(phyloseq_data_rar, x = "map_Id", measures = "Shannon")

##prepare nmds
phyloseq_data_BC <- phyloseq::distance(ekoi_asv_perc, method = "bray", type = "samples")

ekoi_asv_nmds <- phyloseq::ordinate(
  ekoi_asv_perc,
  "NMDS",
  distance = phyloseq_data_BC
)

ekoi_asv_nmds<-phyloseq::plot_ordination(
  ekoi_asv_perc, 
  ekoi_asv_nmds, 
  type = "samples",
  label = "map_Id",
  color = "Sample_Type"
)
plot(ekoi_asv_nmds)

##dendogram
phyloseq_data_BC_hclust_average <- hclust(phyloseq_data_BC, method = "average")

plot(phyloseq_data_BC_hclust_average)

