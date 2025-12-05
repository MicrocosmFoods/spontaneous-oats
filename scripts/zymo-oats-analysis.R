library(tidyverse)
library(colorspace)
library(ggridges)

#################################
# Prep metadata and sylph profiles files
#################################

# MAG metadata from curation repo
mag_metadata_url <- "https://raw.githubusercontent.com/MicrocosmFoods/fermentedfood_metadata_curation/refs/heads/main/data/2025-05-21-genome-metadata-food-taxonomy.tsv"

mag_metadata <- read_tsv(mag_metadata_url) %>% 
  mutate(genome_accession = mag_id) %>% 
  select(genome_accession, completeness, contamination, contigs, taxonomy, species, rep_95id, food_name, main_ingredient, ingredient_group, origin, food_type)

rep_mags_metadata <- mag_metadata %>% 
  filter(genome_accession == rep_95id) %>% 
  select(-rep_95id) %>% 
  mutate(species = case_when(
    is.na(species) | str_to_lower(species) == "unknown" ~ str_c(
      str_extract(taxonomy, "[^;]+$"),
      " spp."
    ),
    TRUE ~ species
  )) %>% 
  select(genome_accession, completeness, contamination, contigs, taxonomy, species)

# sylph profiling results
sylph_profiles <- read_tsv("results/2025-12-02-profiling/combined_sylph_profiles.tsv") %>%
  mutate(accession_name = gsub("_trimmed_1.fastq.gz", "", Sample_file)) %>% 
  mutate(genome_accession = gsub(".fa", "", Genome_file)) %>% 
  select(accession_name, genome_accession, Sequence_abundance, Adjusted_ANI, Eff_cov, Contig_name)

# sample metadata
sample_metadata <- read.csv("metadata/2025-12-01-zymo-oat-sequencing-metadata.csv") %>% 
  mutate(accession_name = gsub("_R1.fastq.gz", "", fastq_1)) %>% 
  select(sample, accession_name, oat, day)

# merge with genome and sample metadata
sylph_profiles_metadata <- left_join(sylph_profiles, rep_mags_metadata) %>% 
  left_join(sample_metadata)

#################################
# Basic summary stats
#################################

# summary stats per sample
sylph_profiles_stats <- sylph_profiles_metadata %>% 
  group_by(sample, oat) %>%
  summarise(
    n_genomes = n_distinct(genome_accession),
    percent_mapped = round(sum(Sequence_abundance, na.rm = TRUE), 3),
    percent_unmapped = round(100 - sum(Sequence_abundance, na.rm = TRUE), 3),
    .groups = "drop"
  )

# prep df for showing abundance of top species
abundance_df_labelled <- sylph_profiles_metadata %>% 
  group_by(sample) %>% 
  arrange(desc(Sequence_abundance), .by_group = TRUE) %>% 
  mutate(
    rank_in_sample = row_number(),
    species_label = if_else(rank_in_sample <=5, species, "Other Species")
  ) %>% 
  ungroup() %>% 
  select(sample, oat, Sequence_abundance, species_label)

## plot of abundance of top species

# color palette prep

okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

pal17 <- colorRampPalette(okabe_ito)(17)

oat_abundance_plot <- abundance_df_labelled %>% 
  ggplot(aes(x=sample, y=Sequence_abundance, fill=species_label)) +
  geom_col() +
  facet_wrap(~ fct_rev(oat), scales = "free_x",) +
  theme_bw() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values = pal17) +
  theme(axis.text.x = element_text(size=14), axis.text.y=element_text(size=14), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), plot.title=element_text(face="bold", size=16), legend.text=element_text(size=14), legend.title=element_text(size=15), strip.text=element_text(size=15), strip.background=element_rect(fill="white", color="black")) +
  labs(
    x="Oat Sample",
    y="% Sequence Abundance",
    fill="Species",
    title="Sequence Abundance of Top Species in Spontaneous Fermented Oat Samples"
  )

ggsave("figures/zymo-oat-sequencing-abundance-plot.png", oat_abundance_plot, width=12, height=7, units=c("in"))

## plot ANI distribution of hits along with colored abundance categories
# first categorize by high, medium, and low sequence abundance to color by
abundance_cat_df <- sylph_profiles_metadata %>%
  mutate(
    abundance_cat = case_when(
      Sequence_abundance > 10 ~ "High (> 10%)",
      Sequence_abundance > 1  ~ "Medium (10% > 1%)",
      TRUE                    ~ "Low (< 1%)"
    ),
    abundance_cat = factor(abundance_cat, levels = c("High (> 10%)", "Medium (10% > 1%)", "Low (< 1%)"))
  )

abundance_cat_df %>% 
  ggplot(aes(x= Adjusted_ANI, y=fct_rev(sample), fill=abundance_cat)) +
  geom_density_ridges(alpha = 0.8, scale = 1.1, color="white", linewidth=0.3) +
  facet_wrap(~ fct_rev(oat), nrow=2, scales="free_y") +
  theme_bw() +
  scale_y_discrete(expand=c(0,0))
