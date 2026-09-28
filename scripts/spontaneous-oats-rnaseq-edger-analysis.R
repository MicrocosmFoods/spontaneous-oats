library(tidyverse)
library(edgeR)
library(splines)
library(msigdbr)

# THP-1 monocytes stimulated with LPS+ and treated with spontaneously fermented oat extracts
# 2 oat bases, 4 sample dates, one control each for PDTC, LPS+, and LPS-, no replicates for any of the samples

counts_table <- read_tsv("results/XKQB9L-expression-matrix-spontaneous-oats.tsv")
sample_map <- read.csv("metadata/spontaneous-oats-rnaseq-sample-map.csv") %>% 
  mutate(sample = paste0(sample_id, "_count")) %>% 
  separate_wider_regex(
    sample_condition,
    patterns = c(oat_type = "[^_]*_[^_]*", "_t", sample_day = ".*"),
    cols_remove = FALSE,
    too_few = "align_start"
  ) %>% 
  mutate(sample_day = as.numeric(sample_day))

counts_matrix <- counts_table %>% 
  column_to_rownames(var = names(counts_table)[1]) %>% 
  select(all_of(sample_map$sample)) %>% 
  as.matrix()

# design for EdgeR
oat <- sample_map$sample_type == "oat"
days <- sort(unique(sample_map$sample_day[oat]))

# 2 degrees of freedom natural spline based on fermentation day
basis <- ns(sample_map$sample_day[oat], df = 2)
S <- matrix(0, nrow(sample_map), 2, dimnames = list(NULL, c("t1","t2")))
S[oat, ] <- basis

design <- cbind(
  LPSneg = as.numeric(sample_map$sample_condition == "LPS-"),
  LPSpos = as.numeric(sample_map$sample_condition == "LPS+"),
  PDTC = as.numeric(sample_map$sample_condition == "PDTC"),
  oat_4 = as.numeric(sample_map$sample_type == "oat_4"),
  oat_10 = as.numeric(sample_map$sample_type == "oat_10"),
  S
)

rownames(design) <- sample_map$sample

# contrast functions
e <- function(nm) setNames(as.numeric(colnames(design) == nm), colnames(design))
# fit data for oat responses over the time series at any day d, and average over the time-series
oat_type_at <- function(d, type_col) {
  v <- setNames(numeric(ncol(design)), colnames(design))
  v[type_col] <- 1
  v[c("t1", "t2")] <- predict(basis, d)
  v
}

outdir <- "results"
save_res <- function(res, name) {
  tt <- topTags(res, n = Inf)$table
  write.csv(tt, file.path(outdir, paste0("DE_", name, ".csv")))
  cat(sprintf("%-28s FDR<0.05: %5d  (up %d / down %d)\n", name,
              sum(tt$FDR < 0.05),
              sum(tt$FDR < 0.05 & tt$logFC > 0),
              sum(tt$FDR < 0.05 & tt$logFC < 0)))
  invisible(tt)
}

# filter, normalization, dispersion
y <- DGEList(counts_matrix, samples = sample_map)
cpm_cut <- 10 / (median(y$samples$lib.size) / 1e6)
keep <- rowSums(cpm(y) > cpm_cut) >= 3
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

y   <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
logcpm <- cpm(y, log = TRUE, prior.count = 2)

# write out normalized log CPM results
write.csv(logcpm, "results/logCPM-normalized-spontaneous-oats-rnaseq.csv", quote = FALSE)

# QC plot checks
plotMDS(y, labels = sample_map$sample, main = "MDS (logCPM)")
plotBCV(y)
plotQLDisp(fit)

# gene-level tests
lps <- save_res(glmQLFTest(fit, contrast = e("LPSpos") - e("LPSneg")), "LPSpos_vs_LPSneg") 
    ## LPSpos_vs_LPSneg             FDR<0.05:  4304  (up 2439 / down 1865)

# suppression control, PDTC vs LPS+
save_res(glmQLFTest(fit, contrast = e("PDTC") - e("LPSpos")), "PDTC_vs_LPSpos")

# each oat fermentation day vs LPS+, fitted curve
for (d in days) {
  for (tc in c("oat_4", "oat_10")) { 
    save_res(glmQLFTest(fit, contrast = oat_type_at(d, tc) - e("LPSpos")),
             sprintf("%s_day%s_vs_LPSpos", tc, d))
  }
}

# genes that change with fermentation day, consistently in both series - don't have replicates so almost asking this instead, but they are very different oat bases/composition
save_res(glmQLFTest(fit, coef = c("t1", "t2")), "oat_time_trend")
