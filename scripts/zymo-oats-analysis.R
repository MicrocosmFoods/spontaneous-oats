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
round1_sylph_profiles <- read_tsv("results/2025-12-02-profiling/combined_sylph_profiles.tsv") %>%
  mutate(accession_name = gsub("_trimmed_1.fastq.gz", "", Sample_file)) %>% 
  mutate(genome_accession = gsub(".fa", "", Genome_file)) %>% 
  select(accession_name, genome_accession, Sequence_abundance, Adjusted_ANI, Eff_cov, Contig_name)

round2_sylph_profiles <- read_tsv("results/2026-01-16-profiling/combined_sylph_profiles.tsv") %>% 
  mutate(accession_name = gsub("_trimmed_1.fastq.gz", "", Sample_file)) %>% 
  mutate(genome_accession = gsub(".fa", "", Genome_file)) %>% 
  select(accession_name, genome_accession, Sequence_abundance, Adjusted_ANI, Eff_cov, Contig_name)

combined_sylph_profiles <- rbind(round1_sylph_profiles, round2_sylph_profiles)

# sample metadata
round1_sample_metadata <- read.csv("metadata/2025-12-01-zymo-oat-sequencing-metadata.csv") %>% 
  mutate(accession_name = gsub("_R1.fastq.gz", "", fastq_1)) %>% 
  select(sample, accession_name, oat, day)

round2_sample_metadata <- read.csv("metadata/2026-01-16-zymo-spontaneous-oats-metadata.csv") %>% 
  mutate(accession_name = gsub("_R1.fastq.gz", "", fastq_1)) %>% 
  mutate(accession_name = gsub(".*/", "", accession_name)) %>% 
  select(sample, accession_name, oat, day)

combined_sample_metadata <- rbind(round1_sample_metadata, round2_sample_metadata)

# merge with genome and sample metadata
sylph_profiles_metadata <- left_join(combined_sylph_profiles, rep_mags_metadata) %>% 
  left_join(combined_sample_metadata) %>% 
  mutate(genus = str_extract(taxonomy, "[^;]+$")) 

genus_map <- tibble::tribble(
  ~pattern,                         ~genus_group,
  "^Bacillus(_.*)?$",               "Bacillus",
  "^Enterococcus(_.*)?$",           "Enterococcus")

sylph_profiles_metadata <- sylph_profiles_metadata %>%
  mutate(
    genus_group = purrr::map_chr(
      genus,
      \(g) {
        hit <- genus_map %>% filter(str_detect(g, pattern))
        if (nrow(hit) > 0) hit$genus_group[[1]] else g
      }
    )
  )

# write out table
write_tsv(sylph_profiles_metadata, "results/combined-spontaneous-oat-profiles.tsv")

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

# prep df for showing abundance of top genera
abundance_df_labelled <- sylph_profiles_metadata %>% 
  group_by(sample) %>% 
  arrange(desc(Sequence_abundance), .by_group = TRUE) %>% 
  mutate(
    rank_in_sample = row_number(),
    
    # robust "missing genus" flag (handles real NA + common string-y NA forms)
    genus_missing = is.na(genus_group) |
      str_trim(as.character(genus_group)) %in% c("", "NA", "N/A", "<NA>", "na", "n/a"),
    
    genus_label = case_when(
      genus_missing        ~ "Other Genera",
      rank_in_sample <= 5  ~ as.character(genus_group),
      TRUE                 ~ "Other Genera"
    )
  ) %>% 
  ungroup() %>% 
  select(sample, oat, Sequence_abundance, genus_label, day)

genus_levels <- abundance_df_labelled %>%
  distinct(genus_label) %>%
  pull(genus_label) %>%
  as.character()

genus_levels <- c(sort(setdiff(genus_levels, "Other Genera")), "Other Genera")

abundance_df_labelled <- abundance_df_labelled %>%
  mutate(genus_label = factor(genus_label, levels = genus_levels))

oat_order <- c("oat_4", "oat_9", "oat_10", "oat_16", "oat_17", "oat_18", "oat_19", "oat_22", "oat_23", "oat_26")

abundance_df_labelled <- abundance_df_labelled %>%
  mutate(
    oat = factor(oat, levels = oat_order)
  )

## plot of abundance of top genera

# color palette prep

okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

okabe_ito_palette <- colorRampPalette(okabe_ito)(22)

non_other <- setdiff(genus_levels, "Other Genera")

genus_colors <- c(
  setNames(okabe_ito_palette[seq_along(non_other)], non_other),
  "Other Genera" = "grey70"
)

# oat sample metadata
oat_order <- c("oat_4", "oat_9", "oat_10", "oat_16", "oat_17", "oat_18", "oat_19", "oat_22", "oat_23", "oat_26")


# plot
oat_abundance_plot <- abundance_df_labelled %>% 
  mutate(day = gsub("day", "", day)) %>% 
  ggplot(aes(x=day, y=Sequence_abundance, fill=genus_label)) +
  geom_col() +
  facet_wrap(~ oat, scales = "free_x", nrow = 2,
             labeller = as_labeller(function(x) x)) +
  theme_bw() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values = genus_colors) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(size=14), axis.text.y=element_text(size=14), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), plot.title=element_text(face="bold", size=16), legend.text=element_text(size=14), legend.title=element_text(size=15), strip.text=element_text(size=15), strip.background=element_rect(fill="white", color="black")) +
  labs(
    x="Sample Day",
    y="% Sequence Abundance",
    fill="Genus",
    title="Sequence Abundance of Top Genera in Spontaneously Fermented Oat Samples"
  )

oat_abundance_plot

ggsave("figures/zymo-oat-sequencing-abundance-plot.png", oat_abundance_plot, width=15, height=9, units=c("in"))

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


#################################
# Plot at different taxonomy levels until LAB/AAB are split out
#################################
sylph_profiles_metadata_tax <- sylph_profiles_metadata %>% 
  select(-genus) %>% 
  separate_wider_delim(
    taxonomy,
    delim=";",
    names=c("phylum", "class", "order", "family", "genus")
  )

oat_order <- c("oat_4", "oat_9", "oat_10", "oat_16", "oat_17", "oat_18", "oat_19", "oat_22", "oat_23", "oat_26")

# function for making plots at different taxonomy levels 
plot_abundance <- function(df,
                           tax_level,
                           palette = "Set3",
                           label = NULL,
                           title = NULL) {
  
  label <- label %||% tax_level
  title <- title %||% paste0(
    "Sequence Abundance of ", label,
    "-Level Taxa in Spontaneously Fermented Oat Samples"
  )
  
  n <- dplyr::n_distinct(df[[tax_level]])
  pal <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(
      min(n, RColorBrewer::brewer.pal.info[palette, "maxcolors"]),
      palette
    )
  )(n)
  
  df %>%
    mutate(day = gsub("day", "", day)) %>%
    ggplot(aes(x = day, y = Sequence_abundance, fill = .data[[tax_level]])) +
    geom_col() +
    facet_wrap(~ oat, scales = "free_x", nrow = 2,
               labeller = as_labeller(function(x) x)) +
    theme_bw() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = pal) +
    guides(fill = guide_legend(ncol = 1)) +
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(face = "bold", size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 15),
          strip.background = element_rect(fill = "white", color = "black")) +
    labs(x = "Sample Day", y = "% Sequence Abundance",
         fill = label, title = title)
}

phylum_plot <- plot_abundance(
  sylph_profiles_metadata_tax, "phylum",
  palette="Paired",
  label="Phylum"
)

class_plot <- plot_abundance(
  sylph_profiles_metadata_tax, "class",
  palette="Paired",
  label="Class"
)

order_plot <- plot_abundance(
  sylph_profiles_metadata_tax, "order",
  palette="Set3",
  label="Order"
)

family_plot <- plot_abundance(
  sylph_profiles_metadata_tax, "family",
  palette="Set3",
  label="Family"
)

ggsave("figures/phylum-level-abundance.png", phylum_plot, width=15, height=9, units=c("in"))
ggsave("figures/class-level-abundance.png", class_plot, width=15, height=9, units=c("in"))
ggsave("figures/order-level-abundance.png", order_plot, width=15, height=9, units=c("in"))
ggsave("figures/family-level-abundance.png", family_plot, width=15, height=9, units=c("in"))
