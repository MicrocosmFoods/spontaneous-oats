library(tidyverse)
library(vegan)

#################################
# PCoA (Bray-Curtis) of spontaneously fermented oat samples
#################################

profiles <- read_tsv("results/combined-spontaneous-oat-profiles.tsv")

# filter low coverage hits, renormalize to 100
filtered_profiles <- profiles %>% 
  filter(Eff_cov > 1) %>% 
  group_by(sample, species) %>% 
  summarise(abundance = sum(Sequence_abundance), .groups="drop") %>% 
  group_by(sample) %>% 
  mutate(abundance = 100 * abundance / sum(abundance)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = species, values_from = abundance, values_fill = 0)

# abundance matrix and smaple metadata for vegan prep
abundance_matrix <- filtered_profiles %>% 
  column_to_rownames("sample") %>% 
  as.matrix()

metadata <- tibble(sample = rownames(abundance_matrix)) %>% 
  left_join(distinct(profiles, sample, oat, day), by = "sample") %>% 
  mutate(day = factor(day, levels = c("day1", "day2", "day3", "day6")))

# Bray-Curtis distance and PCoA
bc <- vegdist(abundance_matrix, method = "bray")

# add = TRUE for negative eigenvalues
pcoa <- cmdscale(bc, k=2, eig=TRUE, add=TRUE)
pct <- round(100 * pcoa$eig[1:2]/ sum(pcoa$eig[pcoa$eig > 0], 1))

scores <- pcoa$points
colnames(scores) <- c("PCoA1", "PCoA2")
scores <- scores %>% 
  as_tibble(rownames = "sample") %>% 
  left_join(metadata, by="sample")

# plot colored by fermentation day
axis_labs <- labs(x = paste0("PCoA1 (", pct[1], "%)"),
                  y = paste0("PCoA2 (", pct[2], "%)"))

oat_cols <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
              "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC")


p_oat_day <- ggplot(scores, aes(PCoA1, PCoA2, color=oat, shape=day)) +
  geom_point(size=3.5) +
  scale_color_manual(values = oat_cols) +
  axis_labs +
  theme_bw()

p_oat_day_traj <- scores %>% 
  arrange(oat, day) %>% 
  ggplot(aes(PCoA1, PCoA2, color = oat)) +
  geom_path(aes(group = oat), linewidth = 0.6, alpha = 0.6,
            arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  geom_point(aes(shape = day), size = 3.5) +
  scale_color_manual(values = oat_cols) +
  axis_labs +
  theme_bw()

p_oat_day_traj

ggsave("figures/bc-pcoa-oat-day.png", p_oat_day, width=8, height=5, units=c("in"))
ggsave("figures/bc-pcoa-oat-day-traj.png", p_oat_day_traj, width=8, height=5, units=c("in"))

