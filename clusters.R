#!/usr/bin/env Rscript
# haplotype_cluster_analysis.R

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

OUT_DIR <- "/mnt/mauricio/haploblock_statistics/"
PLOT_DIR <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE)

save_plot <- function(filename, width = 10, height = 6) {
  ggsave(file.path(PLOT_DIR, filename), width = width, height = height, dpi = 300)
  message("Saved: ", filename)
}

# ── Load data ─────────────────────────────────────────────────────────────────
summary_df <- read.delim(file.path(OUT_DIR, "clusters_summary.tsv"), stringsAsFactors = FALSE)
stats_df   <- read.delim(file.path(OUT_DIR, "block_stats.tsv"),   stringsAsFactors = FALSE)

stats_df <- stats_df %>%
  mutate(
    start          = as.numeric(start),
    end            = as.numeric(end),
    midpoint       = (start + end) / 2,
    block_length   = as.numeric(block_length),
    singleton_rate = singleton_count / n_clusters
  )

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Rank-abundance: cluster size distribution across ALL blocks
# Shows the global long-tail structure
# ══════════════════════════════════════════════════════════════════════════════
cluster_sizes <- summary_df %>%
  distinct(block, cluster_id, cluster_size) %>%
  arrange(desc(cluster_size)) %>%
  mutate(rank = row_number())

ggplot(cluster_sizes, aes(x = rank, y = cluster_size)) +
  geom_point(alpha = 0.3, size = 1, colour = "#4472C4") +
  geom_line(alpha = 0.4, colour = "#4472C4") +
  annotation_logticks(sides = "l") +
  theme_minimal(base_size = 13) +
  labs(title = "Rank-abundance of haplotype cluster sizes — chr21 (all blocks)",
       x = "Rank", y = "Cluster size (log₁₀)")
save_plot("01_rank_abundance_all_blocks.png")

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Cluster size distribution histogram (log-binned)
# ══════════════════════════════════════════════════════════════════════════════
ggplot(cluster_sizes, aes(x = cluster_size)) +
  geom_histogram(bins = 40, fill = "#4472C4", colour = "white", linewidth = 0.2) +
  annotation_logticks(sides = "b") +
  theme_minimal(base_size = 13) +
  labs(title = "Distribution of cluster sizes across all chr21 blocks",
       x = "Cluster size (log₁₀)", y = "Count")
save_plot("02_cluster_size_histogram.png")

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 3 — Shannon entropy along chromosome
# Low entropy = one dominant haplotype; high = many equally common ones
# ══════════════════════════════════════════════════════════════════════════════
ggplot(stats_df, aes(x = midpoint / 1e6, y = shannon_entropy)) +
  geom_point(alpha = 0.5, size = 1.5, colour = "#4472C4") +
  geom_smooth(method = "loess", span = 0.15, colour = "#C0392B", se = TRUE, linewidth = 0.8) +
  theme_minimal(base_size = 13) +
  labs(title = "Haplotype diversity (Shannon entropy) along chr21",
       x = "Chromosomal position (Mb)", y = "Shannon entropy (bits)")
save_plot("03_shannon_entropy_along_chr.png", width = 12)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 4 — Dominance along chromosome
# High dominance = one cluster contains most haplotypes (low diversity)
# ══════════════════════════════════════════════════════════════════════════════
ggplot(stats_df, aes(x = midpoint / 1e6, y = dominance)) +
  geom_point(alpha = 0.5, size = 1.5, colour = "#E67E22") +
  geom_smooth(method = "loess", span = 0.15, colour = "#C0392B", se = TRUE, linewidth = 0.8) +
  theme_minimal(base_size = 13) +
  labs(title = "Cluster dominance (largest cluster / total haplotypes) along chr21",
       x = "Chromosomal position (Mb)", y = "Dominance index")
save_plot("04_dominance_along_chr.png", width = 12)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 5 — Number of clusters per block along chromosome
# ══════════════════════════════════════════════════════════════════════════════
ggplot(stats_df, aes(x = midpoint / 1e6, y = n_clusters)) +
  geom_col(fill = "#4472C4", alpha = 0.7, width = 0.05) +
  theme_minimal(base_size = 13) +
  labs(title = "Number of clusters per block along chr21",
       x = "Chromosomal position (Mb)", y = "Number of clusters")
save_plot("05_n_clusters_along_chr.png", width = 12)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 6 — Block length vs. number of clusters
# ══════════════════════════════════════════════════════════════════════════════
ggplot(stats_df, aes(x = block_length / 1e3, y = n_clusters)) +
  geom_point(alpha = 0.4, size = 1.5, colour = "#4472C4") +
  geom_smooth(method = "lm", colour = "#C0392B", se = TRUE) +
  theme_minimal(base_size = 13) +
  labs(title = "Block length vs. number of clusters",
       x = "Block length (kb)", y = "Number of clusters")
save_plot("06_block_length_vs_n_clusters.png")

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 7 — Singleton rate along chromosome
# High singleton rate = many rare/unique haplotypes
# ══════════════════════════════════════════════════════════════════════════════
ggplot(stats_df, aes(x = midpoint / 1e6, y = singleton_rate)) +
  geom_point(alpha = 0.5, size = 1.5, colour = "#8E44AD") +
  geom_smooth(method = "loess", span = 0.15, colour = "#C0392B", se = TRUE, linewidth = 0.8) +
  theme_minimal(base_size = 13) +
  labs(title = "Singleton cluster rate along chr21",
       x = "Chromosomal position (Mb)", y = "Singleton rate (singletons / total clusters)")
save_plot("07_singleton_rate_along_chr.png", width = 12)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 8 — hap0 vs hap1 representation per block
# Detects phase asymmetry in clustering
# ══════════════════════════════════════════════════════════════════════════════
hap_balance <- summary_df %>%
  group_by(block, start, hap) %>%
  summarise(n = n_distinct(member), .groups = "drop") %>%
  pivot_wider(names_from = hap, values_from = n, values_fill = 0) %>%
  mutate(
    total    = `0` + `1`,
    hap1_frac = `1` / total,
    midpoint  = (as.numeric(start)) / 1e6
  )

ggplot(hap_balance, aes(x = midpoint, y = hap1_frac)) +
  geom_point(alpha = 0.4, size = 1.5, colour = "#27AE60") +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "gray50") +
  geom_smooth(method = "loess", span = 0.15, colour = "#C0392B", se = TRUE, linewidth = 0.8) +
  theme_minimal(base_size = 13) +
  labs(title = "hap1 fraction per block along chr21",
       subtitle = "Dashed line = balanced hap0/hap1",
       x = "Chromosomal position (Mb)", y = "Fraction of hap1 members")
save_plot("08_hap_balance_along_chr.png", width = 12)

# ══════════════════════════════════════════════════════════════════════════════
# Summary table
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── chr21 cluster summary ────────────────────────────────────\n")
cat(sprintf("  Blocks analysed        : %d\n",   nrow(stats_df)))
cat(sprintf("  Total clusters         : %d\n",   sum(stats_df$n_clusters)))
cat(sprintf("  Mean clusters/block    : %.1f\n", mean(stats_df$n_clusters)))
cat(sprintf("  Mean Shannon entropy   : %.2f bits\n", mean(stats_df$shannon_entropy)))
cat(sprintf("  Mean dominance index   : %.2f\n", mean(stats_df$dominance)))
cat(sprintf("  Mean singleton rate    : %.2f\n", mean(stats_df$singleton_rate)))
cat(sprintf("  Plots written to       : %s\n",   PLOT_DIR))
cat("─────────────────────────────────────────────────────────────\n")