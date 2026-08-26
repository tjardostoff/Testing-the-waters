#load required packages

library(tidyverse)
library(vegan)
library(corrplot)
library("UpSetR")
library(magrittr)
library(gtools)


#####Molecular data table preparation#####

asv_counts <- read_tsv("coastclim_spatial_dada2_asvs_counts.tsv")
asv <- read_tsv("coastclim_spatial_dada2_asvs.tsv")

#prepare count table for merge

#put in wide format
asv_counts <- asv_counts %>%
  pivot_wider(names_from = "sample",
              values_from = "nreads")

#sort the columns by name
asv_counts <- asv_counts %>%
  select(1, mixedsort(names(.)[-1]))

#replace NA with 0
asv_counts[is.na(asv_counts)] <- 0

asv_counts <- asv_counts[,c(1,5:131)]

asv_table <- left_join(asv, asv_counts, by = "asv_id")

bold <- asv_table[,c(1:9,21:ncol(asv_table))]

bold <- bold %>%
  separate(
    bold_vsearch_taxonomy,
    into = c("kingdom","phylum","class","order","family","subfamily","genus","species"),
    sep = ","
  )
bold[bold == "None"] <- "Unknown"

stations <- c(
  "exposed_1","blank1","exposed_2","exposed_3","exposed_4","exposed_5","blank2",
  "semi_3","semi_1","blank3","exposed_6","semi_2","blank4","shel_1","shel_2",
  "shel_3","blank5","semi_4","shel_4","shel_5","semi_5","blank6","pojo_1",
  "pojo_2","eke_1","blank7","eke_2","eke_3","extra_1","blank8","extra_2",
  "extra_3","extra_4","blank9"
)
mbstations<- c(
  "exposed_1","exposed_2","exposed_3","exposed_4","exposed_5",
  "semi_3","semi_1","exposed_6","semi_2","shel_1","shel_2",
  "shel_3","semi_4","shel_4","shel_5","semi_5","pojo_1",
  "pojo_2","eke_1","eke_2","eke_3","extra_1","extra_2",
  "extra_3","extra_4")
edna_new_names <- paste0(rep(stations, each = 3), "_", 1:3)
mb_new_names <- paste0(mbstations)

bold_cols <- names(bold)

# find EDNA columns
bold_edna_idx <- grep("^EDNA", bold_cols)

# find MB columns
bold_mb_idx <- grep("^MB", bold_cols)

# replace names
bold_cols[bold_edna_idx] <- edna_new_names

bold_cols[bold_mb_idx] <- mb_new_names


# assign back
names(bold) <- bold_cols

#create data frames for edna and metabarcoding
edna_bold <- bold[,c(1:2,5:14,17:118)]

mb_bold <- bold[,c(1:2,5:14,119:143)]

####Molecular Statistics####

#####edna#####

##create pooled dataframe for eDNA
long_edna_bold <- edna_bold %>%
  pivot_longer(
    cols = matches("_(1|2|3)$"),  # only triplicate columns
    names_to = "sample",
    values_to = "reads"
  ) %>%
  mutate(sample = str_remove(sample, "_[123]$")) %>%  # remove replicate suffix
  group_by(across(-c(sample, reads))) %>%             # keep taxonomy columns
  group_by(sample, .add = TRUE) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop")

#remove blanks
long_edna_bold <- long_edna_bold %>%
  filter(!startsWith(sample, "blank")) 

long_edna_bold$sample <- factor(long_edna_bold$sample, levels = c("exposed_1","exposed_2","exposed_3",
                                                                  "exposed_4","exposed_5","exposed_6",
                                                                  "semi_1","semi_2","semi_3","semi_4",
                                                                  "semi_5","shel_1","shel_2","shel_3",
                                                                  "shel_4","shel_5","extra_1","extra_2",
                                                                  "extra_3","extra_4","eke_1","eke_2",
                                                                  "eke_3", "pojo_1","pojo_2"))

long_edna_bold <- long_edna_bold %>%
  mutate(method = "eDNA") %>%
  filter(reads > 0)



####bulk####

##create pooled dataframe for bulk
long_mb_bold <- mb_bold %>%
  pivot_longer(
    cols = -c(asv_id,total,kingdom, phylum, class, order, family, subfamily, genus, species, bold_vsearch_similarity,bold_vsearch_min_query_coverage),
    names_to = "sample",
    values_to = "reads"
  )
#replace NA with 0
long_mb_bold$phylum[is.na(long_mb_bold$phylum)] <- "unknown"

#remove blanks
long_mb_bold <- long_mb_bold %>%
  filter(!startsWith(sample, "blank"))  

long_mb_bold$sample <- factor(long_mb_bold$sample, levels = c("exposed_1","exposed_2","exposed_3",
                                                              "exposed_4","exposed_5","exposed_6",
                                                              "semi_1","semi_2","semi_3","semi_4",
                                                              "semi_5","shel_1","shel_2","shel_3",
                                                              "shel_4","shel_5","extra_1","extra_2",
                                                              "extra_3","extra_4","eke_1","eke_2",
                                                              "eke_3", "pojo_1","pojo_2"))

long_mb_bold <- long_mb_bold %>%
  mutate(method = "Bulk") %>%
  filter(reads > 0)

####Combine Methods####

molecular_stats <- rbind(long_edna_bold, long_mb_bold)


molecular_stats <- molecular_stats %>%
  mutate(
    species = case_when(
      species != "Unknown" ~ species,
      
      genus != "Unknown" ~ paste0(genus, "_sp."),
      family != "Unknown" ~ paste0(family, "_sp."),
      order != "Unknown" ~ paste0(order, "_sp."),
      class != "Unknown" ~ paste0(class, "_sp."),
      phylum != "Unknown" ~ paste0(phylum, "_sp."),
      kingdom != "Unknown" ~ paste0(kingdom, "_sp."),
      
      TRUE ~ "Unknown sp."
    )
  )

#create dataframe with targets only including higher taxonomic levels

#define contamination and non-targets
cont_bold <- c("Idaea_carvalhoi","Chamobates_cuspidatus","Rhacognathus_americanus",
               "Tanytarsus_usmaensis","Cladopelma_virescens","Orthocladius_oblidens","Mellitidia_sp.",
               "Primates_sp.","Coleoptera_sp.","Diptera_sp.","Entomobryomorpha_sp.",
               "Hemiptera_sp.","Hymenoptera_sp.","Insecta","Arachnida")
non_target_bold <- c("Unknown","Scopalina_ruetzleri","Agalma_elegans",
                     "Tethya_minuta","Rotaria_mento","Rotaria_rotatoria", "Amphibalanus_amphitrite",
                     "Isopoda_sp.", "Placozoa_sp.", "Sulcospira_sp.","Prayidae_sp.",
                     "Macrothrix_sp.","Alpheidae_sp.","Protista","Plantae","Fungi")
unassigned <- c("Mixture","*",NA)



molecular_stats2 <- molecular_stats %>%
  mutate(category = case_when(
    class %in% cont_bold ~ "contamination",
    species %in% cont_bold ~ "contamination",
    species %in% non_target_bold ~ "non_target",
    kingdom %in% non_target_bold ~ "non_target",
    kingdom %in% unassigned ~ "unassigned",
    TRUE ~ "target"
  ))

molecular_stats3 <- molecular_stats2 %>%
  filter(category == "target")

molecular_stats_target <- molecular_stats3 %>%
  group_by(sample,  species,genus,class,order,phylum) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop")

unique(molecular_stats_target$species)

#summary on asvs
asv_summary_molecular_stats <- molecular_stats2 %>%
  filter(!is.na(reads) & reads > 0) %>%
  distinct(method, category, asv_id) %>%
  count(method, category, name = "n_asvs")
asv_summary_molecular_stats

#summary on reads across all samples
reads_summary_molecular_stats <- molecular_stats2 %>%
  group_by(method, category) %>%
  summarise(total_reads = sum(reads, na.rm = TRUE), .groups = "drop")
reads_summary_molecular_stats


####Dissimilarity Analysis####

#####Microscopy Biomass#####

#load data and prepare data frame

biomass_microscopy_species <- read_csv2("Long_Microscopy_Biomass_Species.csv")
names(biomass_microscopy_species)[names(biomass_microscopy_species) == 'Microscopy'] <- 'Species'
biomass_microscopy_species <- biomass_microscopy_species %>%
  mutate(Method = "Microscopy")

biomass_microscopy_species$Site <- factor(biomass_microscopy_species$Site, levels = c("exposed_1","exposed_2","exposed_3",
                                                                                      "exposed_4","exposed_5","exposed_6",
                                                                                      "semi_1","semi_2","semi_3","semi_4",
                                                                                      "semi_5","shel_1","shel_2","shel_3",
                                                                                      "shel_4","shel_5","extra_1","extra_2",
                                                                                      "extra_3","extra_4","eke_1","eke_2",
                                                                                      "eke_3", "pojo_1","pojo_2"))

biomass_microscopy_species2 <- biomass_microscopy_species %>%
  filter(!biomass == 0)

#put table in wide format

biomass_mic_species_diss <- biomass_microscopy_species2 %>%
  distinct() %>% # Removes duplicates
  pivot_wider(names_from = Species, 
              values_from = biomass, 
              values_fill = 0) %>% 
  select(-Method)

biomass_mic_species_diss <- biomass_mic_species_diss %>% remove_rownames %>% column_to_rownames(var='Site')

#transform data

dec_biomass_mic_species_diss <- decostand(biomass_mic_species_diss, method = "total")

#calculate dissimilarity

biomass_mic_dist <- vegdist(dec_biomass_mic_species_diss, method = "bray")

#####Imaging Biomass#####

#load data and prepare data frame

biomass_imaging_species <- read_csv2("Long_Imaging_Biomass_Species.csv")

names(biomass_imaging_species)[names(biomass_imaging_species) == 'Imaging'] <- 'Species'
biomass_imaging_species <- biomass_imaging_species %>%
  mutate(Method = "Imaging")

biomass_imaging_species$Site <- factor(biomass_imaging_species$Site, levels = c("exposed_1","exposed_2","exposed_3",
                                                                                "exposed_4","exposed_5","exposed_6",
                                                                                "semi_1","semi_2","semi_3","semi_4",
                                                                                "semi_5","shel_1","shel_2","shel_3",
                                                                                "shel_4","shel_5","extra_1","extra_2",
                                                                                "extra_3","extra_4","eke_1","eke_2",
                                                                                "eke_3", "pojo_1","pojo_2"))

biomass_imaging_species2 <- biomass_imaging_species %>%
  filter(!biomass == 0)

#put table in wide format
biomass_img_species_diss <- biomass_imaging_species2 %>%
  distinct() %>% # Removes duplicates
  pivot_wider(names_from = Species, 
              values_from = biomass, 
              values_fill = 0) %>% 
  select( -Method)

biomass_img_species_diss <- biomass_img_species_diss %>% remove_rownames %>% column_to_rownames(var='Site')

biomass_img_species_diss <- biomass_img_species_diss[rownames(biomass_mic_species_diss), ]

#transform data

dec_biomass_img_species_diss <- decostand(biomass_img_species_diss, method = "total")

#calculate dissimilarity

biomass_img_dist <- vegdist(dec_biomass_img_species_diss, method = "bray")

#####Molecular#####

#load data and prepare data frame

molecular_species <- read_csv2("Species_BOLD_Long.csv")

names(molecular_species)[names(molecular_species) == 'sample'] <- 'Site'
names(molecular_species)[names(molecular_species) == 'method'] <- 'Method'
names(molecular_species)[names(molecular_species) == 'species'] <- 'Species'
names(molecular_species)[names(molecular_species) == 'reads'] <- 'biomass'

molecular_species$Site <- factor(molecular_species$Site, levels = c("exposed_1","exposed_2","exposed_3",
                                                                    "exposed_4","exposed_5","exposed_6",
                                                                    "semi_1","semi_2","semi_3","semi_4",
                                                                    "semi_5","shel_1","shel_2","shel_3",
                                                                    "shel_4","shel_5","extra_1","extra_2",
                                                                    "extra_3","extra_4","eke_1","eke_2",
                                                                    "eke_3", "pojo_1","pojo_2"))

##Subsets Molecular Methods##
bulk_species <- subset(molecular_species, Method == "bulk")
edna_species <- subset(molecular_species, Method == "eDNA")

#####Bulk#####

bulk_species <- bulk_species[,c(1:3,8)]

bulk_species2 <- bulk_species %>%
  filter(!biomass == 0)

#put table in wide format
bulk_species_diss <- bulk_species2 %>%
  distinct() %>% # Removes duplicates
  pivot_wider(names_from = Species, 
              values_from = biomass, 
              values_fill = 0) %>% 
  select( -Method)

bulk_species_diss <- bulk_species_diss %>% remove_rownames %>% column_to_rownames(var='Site')

bulk_species_diss <- bulk_species_diss[rownames(biomass_mic_species_diss), ]

#transform data

dec_bulk_species_diss <- decostand(bulk_species_diss, method = "total")

#calculate dissimilarity

bulk_dist <- vegdist(dec_bulk_species_diss, method = "bray")


#####eDNA#####

edna_species <- edna_species[,c(1:3,8)]

edna_species2 <- edna_species %>%
  filter(!biomass == 0)

#put table in wide format
edna_species_diss <- edna_species2 %>%
  group_by(Site, Species) %>%
  summarise(biomass = sum(biomass), .groups = "drop") %>%
  pivot_wider(
    names_from = Species,
    values_from = biomass,
    values_fill = 0
  )

edna_species_diss <- edna_species_diss %>% remove_rownames %>% column_to_rownames(var='Site')

edna_species_diss <- edna_species_diss[rownames(biomass_mic_species_diss), ]

#transform data

dec_edna_species_diss <- decostand(edna_species_diss, method = "total")

#calculate dissimilarity

edna_dist <- vegdist(dec_edna_species_diss, method = "bray")

#####compare all dissimilarities#####



biomass_dists <- list(
  Microscopy = biomass_mic_dist,
  Imaging = biomass_img_dist,
  Bulk = bulk_dist,
  eDNA = edna_dist
)

n <- length(biomass_dists)

#create matrix

biomass_results <- matrix(NA, n, n,
                          dimnames = list(names(biomass_dists),
                                          names(biomass_dists)))

for(i in 1:(n-1)){
  for(j in (i+1):n){
    biomass_results[i,j] <- mantel(
      biomass_dists[[i]],
      biomass_dists[[j]]
    )$statistic
  }
}

biomass_results


#plot correlation

corrplot(biomass_results, type = 'upper', method = "number")

#Mantel's test

mantel(biomass_mic_dist, biomass_img_dist)
mantel(biomass_mic_dist, bulk_dist)
mantel(biomass_mic_dist, edna_dist)
mantel(biomass_img_dist, bulk_dist)
mantel(biomass_img_dist, edna_dist)
mantel(bulk_dist, edna_dist)


####biomass plots####

#combine microscopy and imaging data

biomass_morpho_species <- rbind(biomass_microscopy_species, biomass_imaging_species)

######plot top 10 species#####
# Find the 10 species with the highest total biomass morphological
biomass_morpho_species2 <- biomass_morpho_species

biomass_morpho_species2$Species <- recode(biomass_morpho_species2$Species, 
                                          "biomass_ac_bifilosa_ad" = "Acartia_bifilosa",
                                          'biomass_ac_sp_total' = 'Acartia_sp',
                                          'biomass_ac_tonsa_ad' = 'Acartia_tonsa',
                                          'biomass_all_cop_unident' = 'Copepoda',
                                          'biomass_all_cyclo' = 'Cyclopidae',
                                          'biomass_Asplancha_spp' = 'Asplanchna_spp',
                                          'biomass_Balanidae_LV' = 'Balanidae',
                                          'biomass_Bivalvia' = 'Bivalvia',
                                          'biomass_bos_coregoni_total'= 'Bosmina_sp',
                                          'biomass_Cercopagis pengoi'= 'Cercopagis_pengoi',
                                          'biomass_Ceriodaphnia_spp' = 'Ceriodaphnia_spp',
                                          'biomass_Chydorus_spp' = 'Chydorus_spp',
                                          'biomass_Daphnia_christata' = 'Daphnia_christata',
                                          'biomass_Daphnia_cucullata' = 'Daphnia_cucullata',
                                          'biomass_Daphnia_spp' = 'Daphnia_spp',
                                          'biomass_eury_total' = 'Eurytemora_spp',
                                          'biomass_Evadne_spp' = 'Evadne_spp',
                                          'biomass_Gastropoda' = 'Gastropoda',
                                          'biomass_Harpacticoida' = 'Harpacticoida',
                                          'biomass_Helicostomella_subulata' = 'Helicostomella_subulata',
                                          'biomass_Keratella_cochlearis' = 'Keratella_cochlearis',
                                          'biomass_Keratella_cruciformis' = 'Keratella_cruciformis',
                                          'biomass_Keratella_quadrata' = 'Keratella_quadrata',
                                          'biomass_Podon_spp' = 'Podon_spp',
                                          'biomass_Polychaeta' = 'Polychaeta',
                                          'biomass_Sididae' = 'Sididae',
                                          'biomass_Synchaeta_spp' = 'Synchaeta_spp',
                                          'biomass_tem_total' = 'Temora_spp',
                                          'biomass_Tintinnopsis_brandti' = 'Tintinnopsis_brandti',
                                          'biomass_Tintinnopsis_fimbriata' = 'Tintinnopsis_fimbriata',
                                          'biomass_Tintinnopsis_lobiancoi' = 'Tintinnopsis_lobiancoi',
                                          'img_biomass_ac_total' = 'Acartia_sp',
                                          'img_biomass_all_cal' = 'Calanoida',
                                          'img_biomass_all_cop' = 'Copepoda',
                                          'img_biomass_all_cyclo' = 'Cyclopidae',
                                          'img_biomass_Asplancha' = 'Asplanchna_spp',
                                          'img_biomass_Balanidae' = 'Balanidae',
                                          'img_biomass_Bivalvia' = 'Bivalvia',
                                          'img_biomass_Bosmina' = 'Bosmina_sp',
                                          'img_biomass_Cercopagis_pengoi' = 'Cercopagis_pengoi',
                                          'img_biomass_Ceriodaphnia' = 'Ceriodaphnia_spp',
                                          'img_biomass_Daphnia' = 'Daphnia_spp',
                                          'img_biomass_Daphnia_cucullata' = 'Daphnia_cucullata',
                                          'img_biomass_Diplostraca' = 'Diplostraca',
                                          'img_biomass_Eurytemora' = 'Eurytemora_spp',
                                          'img_biomass_Evadne' = 'Evadne_spp',
                                          'img_biomass_Gastropoda' = 'Gastropoda',
                                          'img_biomass_Harpacticoida' = 'Harpacticoida',
                                          'img_biomass_Helicostomella_subulata' = 'Helicostomella_subulata',
                                          'img_biomass_Keratella' = 'Keratella',
                                          'img_biomass_Keratella_cochlearis' = 'Keratella_cochlearis',
                                          'img_biomass_Keratella_cruciformis' = 'Keratella_cruciformis',
                                          'img_biomass_Keratella_quadrata' = 'Keratella_quadrata',
                                          'img_biomass_Podon' = 'Podon_spp',
                                          'img_biomass_Polychaeta' = 'Polychaeta',
                                          'img_biomass_Rotifera' = 'Rotifera',
                                          'img_biomass_Sididae' = 'Sididae',
                                          'img_biomass_Synchaeta' = 'Synchaeta_spp',
                                          'img_biomass_Tintinnopsis' = 'Tintinnopsis_spp')


top10_species_morpho <- biomass_morpho_species2 %>%
  group_by(Species) %>%
  summarise(total_biomass = sum(biomass, na.rm = TRUE)) %>%
  arrange(desc(total_biomass)) %>%
  slice_head(n = 10) %>%
  pull(Species)

# Recode all other species as "Other"
df_top10_morpho <- biomass_morpho_species2 %>%
  mutate(
    Species = ifelse(Species %in% top10_species_morpho,
                     as.character(Species),
                     "Other")
  ) %>%
  group_by(Site, Species, Method) %>%
  summarise(biomass = sum(biomass), .groups = "drop")

#plot top 10 morphological
df_top10_morpho$Species <- factor(df_top10_morpho$Species, 
                                  levels = c("Acartia_tonsa","Acartia_sp","Eurytemora_affinis","Eurytemora_spp",
                                             "Calanoida","Thermocyclops_oithonoides","Cyclopidae",
                                             "Copepoda","Bosmina_sp","Podon_spp","Diaphanosoma_brachyurum",
                                             "Daphnia_cucullata","Diplostraca","Keratella_cochlearis","Keratella_quadrata",
                                             "Synchaeta_spp","Polyarthra_dolichoptera","Hydra_oligactis","Other"
                                  ))

df_top10_morpho$Method <- recode(df_top10_morpho$Method,
                                 "Imaging" = "ZooScan + FlowCam")

top10_morph <- ggplot(df_top10_morpho, aes(fill=Species, y=biomass, x=Site)) + 
  geom_bar(position="stack", stat="identity")+
  ylab("Absolute Biomass [µg/m³]")+
  theme(text = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank())+ #remove x axis ticks
  scale_fill_manual(values  = c("#FF0010","#FFA8BB","#C20088","#4C005C","#E47C96","#FFE100",
                                "#100AFF","#19A405","#FFCC99","#426600","#A8AFB7",
                                "#78E8BD","#78DBE8",
                                "gold",
                                "#F0A3FF","#CE92B7","#E0FF66","#003380","#191919",
                                "#6AC9F6","#6AB2F6","#5EF1F2","#B26EF2","grey",
                                "#94FFB5","#005C31","green",
                                "#2BCE48","#9DCC00","#8F7C00","#00998F","#89D7BA",
                                "#808080","grey","#990000"
  ))+
  labs(fill = "Species (Morphological)")+
  facet_wrap(~Method)

#####plot top 5 molecular####

combine bulk and edna
molec_species <- rbind(bulk_species, edna_species)

top5_species_molec <- molec_species %>%
  group_by(Species) %>%
  summarise(total_biomass = sum(biomass, na.rm = TRUE)) %>%
  arrange(desc(total_biomass)) %>%
  slice_head(n = 5) %>%
  pull(Species)

# Recode all other species as "Other"
df_top5_molec <- molec_species %>%
  mutate(
    Species = ifelse(Species %in% top5_species_molec,
                     as.character(Species),
                     "Other")
  ) %>%
  group_by(Site, Species, Method) %>%
  summarise(biomass = sum(biomass), .groups = "drop")





df_top5_molec$Species <- factor(df_top5_molec$Species, 
                                levels = c("Acartia_tonsa","Acartia_sp","Eurytemora_affinis","Eurytemora_spp",
                                           "Calanoida","Thermocyclops_oithonoides","Cyclopidae",
                                           "Copepoda","Bosmina_sp","Podon_spp","Diaphanosoma_brachyurum",
                                           "Daphnia_cucullata","Diplostraca","Keratella_cochlearis","Keratella_quadrata",
                                           "Synchaeta_spp","Polyarthra_dolichoptera","Hydra_oligactis","Other"
                                ))
df_top5_molec$Method <- recode(df_top5_molec$Method,
                               "bulk" = "Bulk")
top5_mol <- ggplot(df_top5_molec, aes(fill=Species, y=biomass, x=Site)) + 
  geom_bar(position="fill", stat="identity")+
  ylab("Relative number of reads")+
  scale_x_discrete(labels=c('1', '2', '3', '4','5','6','7','8','9','10','11','12',
                            '13','14','15','16','17','18','19','20',
                            '21','22','23','24','25'))+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.3),)+
  scale_fill_manual(values  = c("#990000","#FF0010","#78E8BD","#78DBE8","#FFCC99","#A8AFB7","#4C005C","#100AFF","gold",
                                "#F0A3FF","#CE92B7","#E0FF66","#003380","#191919","grey",
                                "#6AC9F6","#6AB2F6","#5EF1F2","#B26EF2","#FFE100",
                                "#19A405","#94FFB5","#005C31","green",
                                "#2BCE48","#9DCC00","#8F7C00","#426600","#00998F","#89D7BA",
                                "#808080","grey","#A8AFB7","#C20088","#E47C96","#FFA8BB"
  ))+
  labs(fill = "Species (Molecular)")+
  facet_wrap(~Method)

#put together
top15_biomass_combined <- rbind(df_top10_morpho,df_top5_molec)

unique(top15_biomass_combined$Species)

top15_biomass_combined$Species <- factor(top15_biomass_combined$Species, 
                                         levels = c("Acartia_tonsa","Acartia_sp","Eurytemora_spp",
                                                    "Calanoida","Cyclopidae","Copepoda","Bosmina_sp","Podon_spp",
                                                    "Diaphanosoma_brachyurum","Daphnia_cucullata",
                                                    "Diplostraca","Keratella_cochlearis","Keratella_quadrata",
                                                    "Other"
                                         ))

ggplot(top15_biomass_combined, aes(fill=Species, y=biomass, x=Site)) + 
  geom_bar(position="stack", stat="identity")+
  ylab("Number of reads/Biomass [µg/m³]")+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.3),)+
  scale_fill_manual(values  = c("#990000","#FF0010","#78E8BD","#78DBE8","#FFCC99","#4C005C","#100AFF","gold",
                                "#F0A3FF","#CE92B7","#E0FF66","#003380","#191919","grey",
                                "#6AC9F6","#6AB2F6","#5EF1F2","#B26EF2","#FFE100",
                                "#19A405","#94FFB5","#005C31","green",
                                "#2BCE48","#9DCC00","#8F7C00","#426600","#00998F","#89D7BA",
                                "#808080","grey","#A8AFB7","#C20088","#E47C96","#FFA8BB"
  ))+
  facet_wrap(~Method)

#####plot together#####
plot_grid(top10_morph, top5_mol, nrow = 2, align = "v")
top10_morph+top5_mol


# Remove legends from the plots themselves
p_top_noleg <- top10_morph +
  theme(legend.position = "none")

p_bottom_noleg <- top5_mol +
  theme(legend.position = "none")

# Extract each legend
legend_top <- get_legend(top10_morph)
legend_bottom <- get_legend(top5_mol)

# Combine the two legends vertically
combined_legend <- plot_grid(
  legend_top,
  legend_bottom,
  ncol = 1,
  align = "v"
)

# Combine the plots and put the combined legend on the right
final_plot <- plot_grid(
  plot_grid(
    p_top_noleg,
    p_bottom_noleg,
    ncol = 1,
    align = "v"
  ),
  combined_legend,
  ncol = 2,
  rel_widths = c(1, 0.25)
)

final_plot


####Species Accumulation Curves####

#create function for species accumulation curve
tidy_specaccum <- function(x) {
  data.frame(
    site = x$sites,
    richness = x$richness,
    sd = x$sd
  )
}


#####Microscopy#####

#load data and modify data frame

microscopy_species <- read_csv2("Species_Microscopy_Long.csv")
names(microscopy_species)[names(microscopy_species) == 'Microscopy'] <- 'Species'
microscopy_species <- microscopy_species %>%
  mutate(Method = "Microscopy")

microscopy_species$Site <- factor(microscopy_species$Site, levels = c("exposed_1","exposed_2","exposed_3",
                                                                      "exposed_4","exposed_5","exposed_6",
                                                                      "semi_1","semi_2","semi_3","semi_4",
                                                                      "semi_5","shel_1","shel_2","shel_3",
                                                                      "shel_4","shel_5","extra_1","extra_2",
                                                                      "extra_3","extra_4","eke_1","eke_2",
                                                                      "eke_3", "pojo_1","pojo_2"))

#remove all entries without observations

microscopy_species2 <- microscopy_species %>%
  filter(!abundance == 0)

#convert data into wide format

mic_species_pa <- microscopy_species2 %>%
  mutate(Presence = 1) %>%
  distinct(Site, Species, .keep_all = TRUE) %>%
  select(Site, Species, Presence) %>%
  pivot_wider(
    names_from = Species,
    values_from = Presence,
    values_fill = 0
  )

mic_species_pa <- mic_species_pa %>% remove_rownames %>% column_to_rownames(var='Site')

#calculate species accumulation curve

mic_accum <- specaccum(mic_species_pa, method = "random")
plot(mic_accum, ci.type = "polygon", ci.col = "#000004FF",
     main = "Species Accumulation Curve Microscopy",
     xlab = "Sites", ylab = "Species Richness")

tidy_mic_accum <- tidy_specaccum(mic_accum) %>%
  mutate(Method = "Microscopy")


#####Imaging Abundance#####

#load data and modify data frame

imaging_species <- read_csv2("Species_Imaging_Long.csv")

names(imaging_species)[names(imaging_species) == 'Imaging'] <- 'Species'
imaging_species <- imaging_species %>%
  mutate(Method = "Imaging")

imaging_species$Site <- factor(imaging_species$Site, levels = c("exposed_1","exposed_2","exposed_3",
                                                                "exposed_4","exposed_5","exposed_6",
                                                                "semi_1","semi_2","semi_3","semi_4",
                                                                "semi_5","shel_1","shel_2","shel_3",
                                                                "shel_4","shel_5","extra_1","extra_2",
                                                                "extra_3","extra_4","eke_1","eke_2",
                                                                "eke_3", "pojo_1","pojo_2"))

#remove all entries without observations

imaging_species2 <- imaging_species %>%
  filter(!abundance == 0)

#convert data into wide format

img_species_pa <- imaging_species2 %>%
  mutate(Presence = 1) %>%
  distinct(Site, Species, .keep_all = TRUE) %>%
  select(Site, Species, Presence) %>%
  pivot_wider(
    names_from = Species,
    values_from = Presence,
    values_fill = 0
  )

img_species_pa <- img_species_pa %>% remove_rownames %>% column_to_rownames(var='Site')

#calculate species accumulation curve

img_accum <- specaccum(img_species_pa, method = "random")
plot(img_accum, ci.type = "polygon", ci.col = "#721F81FF",
     main = "Species Accumulation Curve Imaging",
     xlab = "Sites", ylab = "Species Richness")

tidy_img_accum <- tidy_specaccum(img_accum) %>%
  mutate(Method = "Imaging")


####Molecular

#load data and modify data frame

molecular_species <- read_csv2("Species_BOLD_Long.csv")

names(molecular_species)[names(molecular_species) == 'sample'] <- 'Site'
names(molecular_species)[names(molecular_species) == 'method'] <- 'Method'
names(molecular_species)[names(molecular_species) == 'species'] <- 'Species'
names(molecular_species)[names(molecular_species) == 'reads'] <- 'abundance'

molecular_species$Site <- factor(molecular_species$Site, levels = c("exposed_1","exposed_2","exposed_3",
                                                                    "exposed_4","exposed_5","exposed_6",
                                                                    "semi_1","semi_2","semi_3","semi_4",
                                                                    "semi_5","shel_1","shel_2","shel_3",
                                                                    "shel_4","shel_5","extra_1","extra_2",
                                                                    "extra_3","extra_4","eke_1","eke_2",
                                                                    "eke_3", "pojo_1","pojo_2"))

##Subsets Molecular Methods##
bulk_species <- subset(molecular_species, Method == "bulk")
edna_species <- subset(molecular_species, Method == "eDNA")

#####Bulk#####

bulk_species <- bulk_species[,c(1:3,8)]

#remove all entries without observations

bulk_species2 <- bulk_species %>%
  filter(!abundance == 0)

#convert data into wide format

bulk_species_pa <- bulk_species2 %>%
  mutate(Presence = 1) %>%
  distinct(Site, Species, .keep_all = TRUE) %>%
  select(Site, Species, Presence) %>%
  pivot_wider(
    names_from = Species,
    values_from = Presence,
    values_fill = 0
  )

bulk_species_pa <- bulk_species_pa %>% remove_rownames %>% column_to_rownames(var='Site')

#calculate species accumulation curve

bulk_accum <- specaccum(bulk_species_pa, method = "random")
plot(bulk_accum, ci.type = "polygon", ci.col = "#F1605DFF",
     main = "Species Accumulation Curve Bulk",
     xlab = "Sites", ylab = "Species Richness")

tidy_bulk_accum <- tidy_specaccum(bulk_accum)%>%
  mutate(Method = "Bulk")


#####eDNA#####

edna_species <- edna_species[,c(1:3,8)]

#remove all entries without observations

edna_species2 <- edna_species %>%
  filter(!abundance == 0)

#convert data into wide format

edna_species_pa <- edna_species2 %>%
  mutate(Presence = 1) %>%
  distinct(Site, Species, .keep_all = TRUE) %>%
  select(Site, Species, Presence) %>%
  pivot_wider(
    names_from = Species,
    values_from = Presence,
    values_fill = 0
  )

edna_species_pa <- edna_species_pa %>% remove_rownames %>% column_to_rownames(var='Site')

#calculate species accumulation curve

edna_accum <- specaccum(edna_species_pa, method = "random")
plot(edna_accum, ci.type = "polygon", ci.col = "#FCFDBFFF",
     main = "Species Accumulation Curve eDNA",
     xlab = "Sites", ylab = "Species Richness")

tidy_edna_accum <- tidy_specaccum(edna_accum) %>%
  mutate(Method = "eDNA")



#####plot all together#####
tidy_accum <- rbind(tidy_mic_accum, tidy_img_accum, tidy_bulk_accum, tidy_edna_accum)

tidy_accum$Method <- factor(tidy_accum$Method, 
                            levels = c("Microscopy","Imaging", "Bulk","eDNA"))

ggplot(tidy_accum, aes(fill=Method, y=richness, x=site)) + 
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(x = site, ymin = richness - 2*sd, ymax = richness + 2*sd), alpha = 0.5)+
  scale_fill_manual(values = c("#000004FF","#721F81FF","#F1605DFF", "#FCFDBFFF"))+
  ylab("Species Richness")+
  xlab("Sites")+
  theme_classic()


####UpSet Plots####

#####All taxa#####

#load data

venn_lists <- read_csv2("Species_Lists_all_methods.csv")

#create subsets for each method

venn_mic <- unique(venn_lists$Microscopy)
venn_mic <- venn_mic[!is.na(venn_mic)]
venn_img <- unique(venn_lists$Imaging)
venn_img <- venn_img[!is.na(venn_img)]
venn_bulk <- unique(venn_lists$Bulk)
venn_bulk <- venn_bulk[!is.na(venn_bulk)]
venn_edna <- unique(venn_lists$eDNA)

#define input list

upsetinput <- list(Microscopy = venn_mic, Imaging = venn_img,
                   Bulk = venn_bulk, eDNA = venn_edna)

#create UpSet Plot

upset_all <- upset(fromList(upsetinput), order.by = "freq", point.size = 2.5, line.size = 1.25, 
                   mainbar.y.label = "Shared Species", sets.x.label = "Species per Method",
                   text.scale=c(2.5), sets.bar.color = c("#F4CAE4","#CBD5E8","#B3E2CD", "#FDCDAC"),
                   set_size.show = FALSE)


#####Pure Species#####

#load data

PURE_venn_lists <- read_csv2("PURE_Species_Lists_all_methods.csv")

#create subsets for each method

PURE_venn_mic <- unique(PURE_venn_lists$Microscopy)
PURE_venn_mic <- PURE_venn_mic[!is.na(PURE_venn_mic)]
PURE_venn_img <- unique(PURE_venn_lists$Imaging)
PURE_venn_img <- PURE_venn_img[!is.na(PURE_venn_img)]
PURE_venn_bulk <- unique(PURE_venn_lists$Bulk)
PURE_venn_bulk <- PURE_venn_bulk[!is.na(PURE_venn_bulk)]
PURE_venn_edna <- unique(PURE_venn_lists$eDNA)

#define input list

PUREupsetinput <- list(Microscopy = PURE_venn_mic, Imaging = PURE_venn_img,
                       Bulk = PURE_venn_bulk, eDNA = PURE_venn_edna)

#create UpSet Plot

upset_spec <- upset(fromList(PUREupsetinput), order.by = "freq", point.size = 2.5, line.size = 1.25, 
                    mainbar.y.label = "Shared Species", sets.x.label = "Species per Method",
                    text.scale=c(2.5), sets.bar.color = c("#F4CAE4","#CBD5E8","#B3E2CD", "#FDCDAC"),
                    set_size.show = FALSE)