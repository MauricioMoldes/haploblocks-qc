#!/usr/bin/env Rscript
# haploblock_qc.R
# Genome-wide haploblock size and structure QC
# Reads from the block_stats.tsv produced by parse_clusters_all_chr.sh

library(ggplot2)
library(dplyr)
library(readr)

# ── Paths ─────────────────────────────────────────────────────────────────────
OUT_DIR  <- "/mnt/mauricio/haploblocks-qc-clean"
PLOT_DIR <- file.path(OUT_DIR, "plots", "qc")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(filename, plot = last_plot(), width = 10, height = 5) {
  ggsave(file.path(PLOT_DIR, paste0(filename, ".png")),
         plot, width = width, height = height, dpi = 300)
  ggsave(file.path(PLOT_DIR, paste0(filename, ".pdf")),
         plot, width = width, height = height)
  message("Saved: ", filename)
}

# ── Load & prepare ────────────────────────────────────────────────────────────
message("Loading block_stats.tsv ...")
df <- read_tsv(file.path(OUT_DIR, "block_stats.tsv"), show_col_types = FALSE)

# Chromosome factor in correct genomic order
chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

df <- df %>%
  mutate(
    chr          = factor(chr, levels = chr_levels),
    start        = as.numeric(start),
    end          = as.numeric(end),
    block_length = as.numeric(block_length),
    log_size     = log10(block_length),
    midpoint     = (start + end) / 2
  ) %>%
  filter(!is.na(block_length), block_length > 0)

message(nrow(df), " blocks loaded across ", n_distinct(df$chr), " chromosomes")

# ── Shared theme ──────────────────────────────────────────────────────────────
qc_theme <- theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# =============================================================================
# 01 — Global size distribution
# =============================================================================
p1 <- ggplot(df, aes(x = log_size)) +
  geom_histogram(bins = 80, fill = "#4472C4", colour = "white", linewidth = 0.2) +
  theme_minimal(base_size = 13) +
  labs(title = "Haploblock size distribution — all chromosomes",
       x = "log₁₀(block length, bp)", y = "Count")
save_plot("01_size_distribution", p1)

# =============================================================================
# 02 — Per-chromosome size boxplot
# =============================================================================
p2 <- ggplot(df, aes(x = chr, y = log_size)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4, fill = "#B5D4F4",
               colour = "#185FA5", linewidth = 0.4) +
  qc_theme +
  labs(title = "Haploblock size per chromosome",
       x = NULL, y = "log₁₀(block length, bp)")
save_plot("02_size_per_chr_boxplot", p2)

# =============================================================================
# 03 — ECDF of block sizes
# =============================================================================
p3 <- ggplot(df, aes(x = block_length)) +
  stat_ecdf(geom = "step", colour = "#4472C4", linewidth = 0.8) +
  scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_si("b"))) +
  theme_minimal(base_size = 13) +
  labs(title = "Cumulative distribution of haploblock sizes",
       x = "Block length (bp, log scale)", y = "Proportion of blocks")
save_plot("03_ecdf", p3)

# =============================================================================
# 04 — Number of blocks per chromosome
# =============================================================================
chr_counts <- df %>%
  count(chr, name = "n_blocks")

p4 <- ggplot(chr_counts, aes(x = chr, y = n_blocks)) +
  geom_col(fill = "#1D9E75", alpha = 0.85) +
  qc_theme +
  labs(title = "Number of haploblocks per chromosome",
       x = NULL, y = "Block count")
save_plot("04_blocks_per_chr", p4)

# =============================================================================
# 05 — Mean block size per chromosome (log scale)
# =============================================================================
chr_mean <- df %>%
  group_by(chr) %>%
  summarise(mean_size = mean(block_length),
            median_size = median(block_length), .groups = "drop")

p5 <- ggplot(chr_mean, aes(x = chr, y = mean_size)) +
  geom_col(fill = "#8E44AD", alpha = 0.85) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_si("b"))) +
  qc_theme +
  labs(title = "Mean haploblock size per chromosome",
       x = NULL, y = "Mean block length (bp, log scale)")
save_plot("05_mean_size_per_chr", p5)

# =============================================================================
# 06 — Block size vs genomic position (faceted)
# =============================================================================
p6 <- ggplot(df, aes(x = midpoint / 1e6, y = block_length)) +
  geom_point(alpha = 0.25, size = 0.4, colour = "#4472C4") +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_si("b"))) +
  facet_wrap(~chr, scales = "free_x", ncol = 4) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(size = 7)) +
  labs(title = "Haploblock size vs genomic position",
       x = "Position (Mb)", y = "Block length (bp, log scale)")
save_plot("06_size_vs_position", p6, width = 16, height = 20)

# =============================================================================
# 07 — Haploblock density (blocks per Mb) per chromosome
# =============================================================================
chr_density <- df %>%
  group_by(chr) %>%
  summarise(density = n() / (max(end) / 1e6), .groups = "drop")

p7 <- ggplot(chr_density, aes(x = chr, y = density)) +
  geom_col(fill = "#E67E22", alpha = 0.85) +
  qc_theme +
  labs(title = "Haploblock density per chromosome",
       x = NULL, y = "Blocks per Mb")
save_plot("07_density_per_chr", p7)

# =============================================================================
# 08 — NEW: Block length distribution as violin per chromosome
#      Richer than the boxplot — shows full shape including bimodality
# =============================================================================
p8 <- ggplot(df, aes(x = chr, y = log_size)) +
  geom_violin(fill = "#4472C4", alpha = 0.35, linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, outlier.alpha = 0.3,
               fill = "white", linewidth = 0.4) +
  qc_theme +
  labs(title = "Haploblock size distribution shape per chromosome",
       subtitle = "Violin shows full distribution; box shows median and IQR",
       x = NULL, y = "log₁₀(block length, bp)")
save_plot("08_size_violin_per_chr", p8, width = 14, height = 6)

# =============================================================================
# 09 — NEW: Cumulative genome coverage by block size threshold
#      "What fraction of the genome is covered by blocks ≥ X bp?"
# =============================================================================
coverage_df <- df %>%
  arrange(desc(block_length)) %>%
  mutate(
    cumulative_bp   = cumsum(block_length),
    total_bp        = sum(block_length),
    coverage_frac   = cumulative_bp / total_bp,
    block_rank      = row_number()
  )

p9 <- ggplot(coverage_df, aes(x = block_length, y = coverage_frac)) +
  geom_line(colour = "#4472C4", linewidth = 0.8) +
  scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_si("b"))) +
  scale_y_continuous(labels = scales::percent) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = "dashed",
             colour = "#C0392B", linewidth = 0.5) +
  annotate("text", x = min(coverage_df$block_length) * 2,
           y = 0.52, label = "50%", size = 3.5, colour = "#C0392B") +
  annotate("text", x = min(coverage_df$block_length) * 2,
           y = 0.92, label = "90%", size = 3.5, colour = "#C0392B") +
  theme_minimal(base_size = 13) +
  labs(title = "Genome coverage by haploblock size",
       subtitle = "What fraction of total haploblock bp is contributed by blocks ≥ X bp?",
       x = "Block length threshold (bp, log scale)",
       y = "Cumulative coverage fraction")
save_plot("09_cumulative_coverage", p9)

# =============================================================================
# Summary table to stdout + TSV
# =============================================================================
summary_tbl <- df %>%
  group_by(chr) %>%
  summarise(
    n_blocks     = n(),
    total_bp     = sum(block_length),
    mean_bp      = round(mean(block_length)),
    median_bp    = round(median(block_length)),
    min_bp       = min(block_length),
    max_bp       = max(block_length),
    density_per_mb = round(n() / (max(end) / 1e6), 2),
    .groups = "drop"
  )

cat("\n── Genome-wide haploblock QC summary ────────────────────────────────────\n")
print(as.data.frame(summary_tbl), row.names = FALSE)
cat(sprintf("\nTotal blocks   : %d\n", sum(summary_tbl$n_blocks)))
cat(sprintf("Total bp       : %s\n",  format(sum(summary_tbl$total_bp), big.mark = ",")))
cat(sprintf("Plots saved to : %s\n",  PLOT_DIR))
cat("─────────────────────────────────────────────────────────────────────────\n")

write_tsv(summary_tbl, file.path(OUT_DIR, "haploblock_qc_summary.tsv"))
message("QC summary table saved: haploblock_qc_summary.tsv")