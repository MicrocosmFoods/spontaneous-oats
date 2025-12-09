library(tidyverse)
library(colorspace)
library(ggridges)
library(scales)
library(tidytext)

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

#################################
# Random subsampling at 5M, 10M, 15M, 20M, and 25M read depth
# Analyze how detection levels of species change with different read depths
# The maximum read depth you can select at Zymo is 20M, and it could be higher than that, but still subsampled at 25M as well
# Our samples were selected for 20M read depth, they came back ranging from min 60M reads to 160M reads
# Now want to ask minimum coverage needed to retain detection of key species to sequence the remaining samples

# Ask how many genomes detected drop out/remain per read depth
# Ask how coverage/ANI changes for the detected genomes, probably need to require minimum 10X depth for the more abundant genomes for accuracy purposes downstream
#################################

subsampling_sylph_profiles <- read_tsv("results/2025-12-08-subsampled-profiling/combined_sylph_profiles.tsv") %>% 
  mutate(accession_name = str_remove(Sample_file, "_n.*$")) %>% 
  mutate(depth = str_match(Sample_file, "_n([^_]+)_")[,2]) %>% 
  mutate(genome_accession = gsub(".fa", "", Genome_file)) %>% 
  select(accession_name, depth, genome_accession, Sequence_abundance, Adjusted_ANI, Eff_cov, Contig_name)
  
subsampled_sylph_profiles_metadata <- left_join(subsampling_sylph_profiles, rep_mags_metadata) %>% 
  left_join(sample_metadata)

# summary stats per sample, depth 
subsampled_sylph_profiles_stats <- subsampled_sylph_profiles_metadata %>% 
  group_by(sample, depth, oat) %>%
  summarise(
    n_genomes = n_distinct(genome_accession),
    percent_mapped = round(sum(Sequence_abundance, na.rm = TRUE), 3),
    percent_unmapped = round(100 - sum(Sequence_abundance, na.rm = TRUE), 3),
    .groups = "drop"
  )

# all stats for the original samples and the subsampled samples
all_stats <- sylph_profiles_stats %>% 
  mutate(depth = "original") %>% 
  rbind(subsampled_sylph_profiles_stats)

# plot for dropout of # of genomes detected 
depth_order  <- c("original", "5000000", "10000000", "15000000", "20000000", "25000000")
sample_order <- c("oat_4_d1", "oat_4_d6", "oat_10_d1", "oat_10_d6")

all_stats_ordered <- all_stats %>% 
  mutate(
    depth  = factor(depth, levels = depth_order, ordered = TRUE),
    sample = factor(sample, levels = sample_order, ordered = TRUE)
  )

all_stats_ordered %>% 
  ggplot(aes(x=depth, y=n_genomes)) +
  geom_col() +
  facet_grid(~ sample) +
  scale_y_continuous(expand = c(0,0))

# heatmap of coverage for detected genomes with changing depth
all_sylph_profiles_metadata <- sylph_profiles_metadata %>% 
  mutate(depth = "original") %>% 
  rbind(subsampled_sylph_profiles_metadata) %>% 
  mutate(
    depth  = factor(depth, levels = depth_order, ordered = TRUE),
    sample = factor(sample, levels = sample_order, ordered = TRUE)
  ) %>% 
  mutate(genome_id = paste0(genome_accession, "_", species))

all_sylph_profiles_metadata_ordered <- all_sylph_profiles_metadata %>%
  group_by(sample, genome_id) %>%
  # compute "original" coverage summary used for ordering
  mutate(orig_cov = max(Eff_cov[depth == "original"], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    genome_id_ord = reorder_within(genome_id, orig_cov, sample)
  )


coverage_profiles <- all_sylph_profiles_metadata_ordered %>% 
  ggplot(aes(x = depth, y = genome_id_ord)) +
  geom_tile(aes(fill = Eff_cov)) +
  facet_wrap(~ sample, scales = "free_y") +
  scale_fill_viridis_c(
    trans = "log1p",
    breaks = pretty_breaks(6)
  ) +
  scale_y_reordered(expand = c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 80, hjust=1))

ggsave("figures/oat-sequencing-random-subsampling-covg-profiles.png", coverage_profiles, width=15, height=7, units=c("in"))
