#!/usr/bin/env Rscript
# haplotype_cluster_analysis_all_chr.R

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

OUT_DIR  <- "/mnt/mauricio/haploblocks-qc-clean"
PLOT_DIR <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE)

save_plot <- function(filename, width = 12, height = 7) {
  ggsave(file.path(PLOT_DIR, filename), width = width, height = height, dpi = 300)
  message("Saved: ", filename)
}

# Chromosome ordering helper (chr1, chr2 ... chr22, chrX, chrY, chrM)
chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

# ── Load data ─────────────────────────────────────────────────────────────────
message("Loading data...")
summary_df <- read.delim(file.path(OUT_DIR, "clusters_summary.tsv"), stringsAsFactors = FALSE)
stats_df   <- read.delim(file.path(OUT_DIR, "block_stats.tsv"),      stringsAsFactors = FALSE)

stats_df <- stats_df %>%
  mutate(
    chr            = factor(chr, levels = chr_levels),
    start          = as.numeric(start),
    end            = as.numeric(end),
    midpoint       = (start + end) / 2,
    block_length   = as.numeric(block_length),
    singleton_rate = singleton_count / n_clusters
  )

summary_df <- summary_df %>%
  mutate(chr = factor(chr, levels = chr_levels))

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Blocks per chromosome (overview)
# ══════════════════════════════════════════════════════════════════════════════
stats_df %>%
  count(chr) %>%
  ggplot(aes(x = chr, y = n)) +
  geom_col(fill = "#4472C4", alpha = 0.85) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of haplotype blocks per chromosome",
       x = NULL, y = "Number of blocks")
save_plot("01_blocks_per_chromosome.png", width = 12, height = 5)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Shannon entropy distribution per chromosome (violin + boxplot)
# Reveals chromosomes with unusual diversity patterns
# ══════════════════════════════════════════════════════════════════════════════
ggplot(stats_df, aes(x = chr, y = shannon_entropy)) +
  geom_violin(fill = "#4472C4", alpha = 0.4, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, outlier.alpha = 0.3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Shannon entropy distribution per chromosome",
       x = NULL, y = "Shannon entropy (bits)")
save_plot("02_entropy_per_chromosome.png", width = 14, height = 6)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 3 — Dominance distribution per chromosome
# ══════════════════════════════════════════════════════════════════════════════
ggplot(stats_df, aes(x = chr, y = dominance)) +
  geom_violin(fill = "#E67E22", alpha = 0.4, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, outlier.alpha = 0.3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cluster dominance distribution per chromosome",
       subtitle = "Dominance = largest cluster size / total haplotypes in block",
       x = NULL, y = "Dominance index")
save_plot("03_dominance_per_chromosome.png", width = 14, height = 6)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 4 — Genome-wide entropy (one panel per chromosome, position on x)
# ══════════════════════════════════════════════════════════════════════════════
ggplot(stats_df, aes(x = midpoint / 1e6, y = shannon_entropy)) +
  geom_point(size = 0.4, alpha = 0.4, colour = "#4472C4") +
  geom_smooth(method = "loess", span = 0.2, se = FALSE,
              colour = "#C0392B", linewidth = 0.6) +
  facet_wrap(~chr, scales = "free_x", ncol = 4) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(size = 7)) +
  labs(title = "Shannon entropy along each chromosome",
       x = "Position (Mb)", y = "Shannon entropy (bits)")
save_plot("04_entropy_genome_wide.png", width = 16, height = 20)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 5 — Genome-wide dominance (faceted)
# ══════════════════════════════════════════════════════════════════════════════
ggplot(stats_df, aes(x = midpoint / 1e6, y = dominance)) +
  geom_point(size = 0.4, alpha = 0.4, colour = "#E67E22") +
  geom_smooth(method = "loess", span = 0.2, se = FALSE,
              colour = "#C0392B", linewidth = 0.6) +
  facet_wrap(~chr, scales = "free_x", ncol = 4) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(size = 7)) +
  labs(title = "Cluster dominance along each chromosome",
       x = "Position (Mb)", y = "Dominance index")
save_plot("05_dominance_genome_wide.png", width = 16, height = 20)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 6 — Singleton rate per chromosome
# ══════════════════════════════════════════════════════════════════════════════
ggplot(stats_df, aes(x = chr, y = singleton_rate)) +
  geom_violin(fill = "#8E44AD", alpha = 0.4, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, outlier.alpha = 0.3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Singleton cluster rate per chromosome",
       x = NULL, y = "Singleton rate")
save_plot("06_singleton_rate_per_chromosome.png", width = 14, height = 6)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 7 — hap0 vs hap1 balance per chromosome
# ══════════════════════════════════════════════════════════════════════════════
hap_balance <- summary_df %>%
  group_by(chr, block, hap) %>%
  summarise(n = n_distinct(member), .groups = "drop") %>%
  pivot_wider(names_from = hap, values_from = n, values_fill = 0) %>%
  rename(hap0 = `0`, hap1 = `1`) %>%
  mutate(
    total     = hap0 + hap1,
    hap1_frac = hap1 / total
  )

ggplot(hap_balance, aes(x = chr, y = hap1_frac)) +
  geom_violin(fill = "#27AE60", alpha = 0.4, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "gray40") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "hap1 fraction per block by chromosome",
       subtitle = "Dashed = balanced hap0/hap1",
       x = NULL, y = "hap1 fraction")
save_plot("07_hap_balance_per_chromosome.png", width = 14, height = 6)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 8 — Rank-abundance global (all chromosomes, coloured by chr)
# ══════════════════════════════════════════════════════════════════════════════
cluster_sizes <- summary_df %>%
  distinct(chr, block, cluster_id, cluster_size) %>%
  arrange(desc(cluster_size)) %>%
  mutate(rank = row_number())

ggplot(cluster_sizes, aes(x = rank, y = cluster_size, colour = chr)) +
  geom_point(size = 0.5, alpha = 0.3) +
  scale_y_log10() +
  scale_colour_viridis_d(option = "turbo") +
  annotation_logticks(sides = "l") +
  theme_minimal(base_size = 12) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 2)) +
  labs(title = "Rank-abundance of cluster sizes — all chromosomes",
       x = "Rank", y = "Cluster size (log₁₀)", colour = NULL)
save_plot("08_rank_abundance_all_chromosomes.png", width = 12, height = 7)


# ══════════════════════════════════════════════════════════════════════════════
# PLOT 9 — Number of clusters per block, faceted by chromosome
# Points + LOESS smoother scales better than bars for dense chromosomes
# ══════════════════════════════════════════════════════════════════════════════
ggplot(stats_df, aes(x = midpoint / 1e6, y = n_clusters)) +
  geom_point(size = 0.4, alpha = 0.35, colour = "#4472C4") +
  geom_smooth(method = "loess", span = 0.2, se = FALSE,
              colour = "#C0392B", linewidth = 0.6) +
  facet_wrap(~chr, scales = "free_x", ncol = 4) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(size = 7)) +
  labs(title = "Number of clusters per block along each chromosome",
       x = "Position (Mb)", y = "Number of clusters")
save_plot("09_n_clusters_genome_wide.png", width = 16, height = 20)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 10 — Block length vs number of clusters, faceted by chromosome
# Log x-axis + LOESS replaces the linear model — relationship is nonlinear
# ══════════════════════════════════════════════════════════════════════════════
ggplot(stats_df, aes(x = block_length / 1e3, y = n_clusters)) +
  geom_point(size = 0.4, alpha = 0.35, colour = "#4472C4") +
  geom_smooth(method = "loess", se = FALSE,
              colour = "#C0392B", linewidth = 0.6) +
  scale_x_log10() +
  facet_wrap(~chr, scales = "free_x", ncol = 4) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(size = 7)) +
  labs(title = "Block length vs number of clusters per chromosome",
       x = "Block length (kb, log₁₀)", y = "Number of clusters")
save_plot("10_block_length_vs_clusters_genome_wide.png", width = 16, height = 20)

# ══════════════════════════════════════════════════════════════════════════════
# Summary table
# ══════════════════════════════════════════════════════════════════════════════
chr_summary <- stats_df %>%
  group_by(chr) %>%
  summarise(
    n_blocks        = n(),
    total_clusters  = sum(n_clusters),
    mean_entropy    = round(mean(shannon_entropy), 3),
    mean_dominance  = round(mean(dominance), 3),
    mean_singletons = round(mean(singleton_rate), 3),
    .groups = "drop"
  )

cat("\n── Genome-wide cluster summary ──────────────────────────────────────────\n")
print(as.data.frame(chr_summary), row.names = FALSE)
cat(sprintf("\nTotal blocks : %d\n", sum(chr_summary$n_blocks)))
cat(sprintf("Total clusters: %d\n", sum(chr_summary$total_clusters)))
cat(sprintf("Plots written to: %s\n", PLOT_DIR))
cat("─────────────────────────────────────────────────────────────────────────\n")

write.table(chr_summary,
            file = file.path(OUT_DIR, "chr_summary_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Summary table saved: chr_summary_table.tsv")