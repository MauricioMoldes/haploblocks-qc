

#!/usr/bin/env Rscript
# haplotype_cluster_analysis.R
# Reads genome-wide files ONCE, then iterates over all chromosomes

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

OUT_DIR <- "/mnt/mauricio/haploblocks-qc-clean"

chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

# ── Load ONCE ─────────────────────────────────────────────────────────────────
message("Loading clusters_summary.tsv ...")
summary_all <- read_tsv(file.path(OUT_DIR, "clusters_summary.tsv"),
                        show_col_types = FALSE) %>%
  mutate(chr = factor(chr, levels = chr_levels))

message("Loading block_stats.tsv ...")
stats_all <- read_tsv(file.path(OUT_DIR, "block_stats.tsv"),
                      show_col_types = FALSE) %>%
  mutate(
    chr          = factor(chr, levels = chr_levels),
    start        = as.numeric(start),
    end          = as.numeric(end),
    midpoint     = (start + end) / 2,
    block_length = as.numeric(block_length),
    singleton_rate = singleton_count / n_clusters
  )

message(nrow(stats_all), " blocks loaded across ",
        n_distinct(stats_all$chr), " chromosomes")

# ── Plotting function — called once per chromosome ─────────────────────────────
plot_chromosome <- function(chr_name, summary_df, stats_df) {

  PLOT_DIR <- file.path(OUT_DIR, "plots", chr_name)
  dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

  save_plot <- function(filename, plot = last_plot(), width = 10, height = 6) {
    path <- file.path(PLOT_DIR, filename)
    ggsave(path, plot, width = width, height = height, dpi = 300)
    message("  Saved: ", filename)
  }

  pfx <- function(n) paste0(chr_name, "_", n)

  # ── Plot 1 — Rank-abundance ────────────────────────────────────────────────
  cluster_sizes <- summary_df %>%
    distinct(block, cluster_id, cluster_size) %>%
    arrange(desc(cluster_size)) %>%
    mutate(rank = row_number())

  p1 <- ggplot(cluster_sizes, aes(x = rank, y = cluster_size)) +
    geom_point(alpha = 0.3, size = 1, colour = "#4472C4") +
    geom_line(alpha = 0.4, colour = "#4472C4") +
    scale_y_log10() +
    annotation_logticks(sides = "l") +
    theme_minimal(base_size = 13) +
    labs(title = paste("Rank-abundance of cluster sizes —", chr_name),
         x = "Rank", y = "Cluster size (log₁₀)")
  save_plot(pfx("01_rank_abundance.png"), p1)

  # ── Plot 2 — Cluster size histogram ───────────────────────────────────────
  p2 <- ggplot(cluster_sizes, aes(x = cluster_size)) +
    geom_histogram(bins = 40, fill = "#4472C4", colour = "white", linewidth = 0.2) +
    scale_x_log10() +
    annotation_logticks(sides = "b") +
    theme_minimal(base_size = 13) +
    labs(title = paste("Cluster size distribution —", chr_name),
         x = "Cluster size (log₁₀)", y = "Count")
  save_plot(pfx("02_cluster_size_histogram.png"), p2)

  # ── Plot 3 — Shannon entropy along chromosome ──────────────────────────────
  p3 <- ggplot(stats_df, aes(x = midpoint / 1e6, y = shannon_entropy)) +
    geom_point(alpha = 0.5, size = 1.5, colour = "#4472C4") +
    geom_smooth(method = "loess", span = 0.15, colour = "#C0392B",
                se = TRUE, linewidth = 0.8) +
    theme_minimal(base_size = 13) +
    labs(title = paste("Shannon entropy —", chr_name),
         x = "Position (Mb)", y = "Shannon entropy (bits)")
  save_plot(pfx("03_shannon_entropy.png"), p3, width = 12)

  # ── Plot 4 — Dominance along chromosome ───────────────────────────────────
  p4 <- ggplot(stats_df, aes(x = midpoint / 1e6, y = dominance)) +
    geom_point(alpha = 0.5, size = 1.5, colour = "#E67E22") +
    geom_smooth(method = "loess", span = 0.15, colour = "#C0392B",
                se = TRUE, linewidth = 0.8) +
    theme_minimal(base_size = 13) +
    labs(title = paste("Cluster dominance —", chr_name),
         subtitle = "Largest cluster / total haplotypes per block",
         x = "Position (Mb)", y = "Dominance index")
  save_plot(pfx("04_dominance.png"), p4, width = 12)

  # ── Plot 5 — Number of clusters per block (points + smoother) ─────────────
  p5 <- ggplot(stats_df, aes(x = midpoint / 1e6, y = n_clusters)) +
    geom_point(alpha = 0.4, size = 1.2, colour = "#4472C4") +
    geom_smooth(method = "loess", span = 0.15, colour = "#C0392B",
                se = TRUE, linewidth = 0.8) +
    theme_minimal(base_size = 13) +
    labs(title = paste("Number of clusters per block —", chr_name),
         x = "Position (Mb)", y = "Number of clusters")
  save_plot(pfx("05_n_clusters_along_chr.png"), p5, width = 12)

  # ── Plot 6 — Block length vs number of clusters (loess, log x) ────────────
  p6 <- ggplot(stats_df, aes(x = block_length / 1e3, y = n_clusters)) +
    geom_point(alpha = 0.35, size = 1.2, colour = "#4472C4") +
    geom_smooth(method = "loess", colour = "#C0392B", se = TRUE, linewidth = 0.8) +
    scale_x_log10() +
    annotation_logticks(sides = "b") +
    theme_minimal(base_size = 13) +
    labs(title = paste("Block length vs number of clusters —", chr_name),
         x = "Block length (kb, log₁₀)", y = "Number of clusters")
  save_plot(pfx("06_block_length_vs_n_clusters.png"), p6)

  # ── Plot 7 — Singleton rate along chromosome ───────────────────────────────
  p7 <- ggplot(stats_df, aes(x = midpoint / 1e6, y = singleton_rate)) +
    geom_point(alpha = 0.5, size = 1.5, colour = "#8E44AD") +
    geom_smooth(method = "loess", span = 0.15, colour = "#C0392B",
                se = TRUE, linewidth = 0.8) +
    theme_minimal(base_size = 13) +
    labs(title = paste("Singleton cluster rate —", chr_name),
         x = "Position (Mb)", y = "Singleton rate")
  save_plot(pfx("07_singleton_rate.png"), p7, width = 12)

  # ── Plot 8 — hap0/hap1 balance ────────────────────────────────────────────
  hap_balance <- summary_df %>%
    group_by(block, start, hap) %>%
    summarise(n = n_distinct(member), .groups = "drop") %>%
    pivot_wider(names_from = hap, values_from = n, values_fill = 0) %>%
    mutate(
      total     = `0` + `1`,
      hap1_frac = `1` / total,
      midpoint  = as.numeric(start) / 1e6
    )

  p8 <- ggplot(hap_balance, aes(x = midpoint, y = hap1_frac)) +
    geom_point(alpha = 0.4, size = 1.5, colour = "#27AE60") +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "gray50") +
    geom_smooth(method = "loess", span = 0.15, colour = "#C0392B",
                se = TRUE, linewidth = 0.8) +
    theme_minimal(base_size = 13) +
    labs(title = paste("hap1 fraction per block —", chr_name),
         subtitle = "Dashed line = balanced hap0/hap1",
         x = "Position (Mb)", y = "Fraction of hap1 members")
  save_plot(pfx("08_hap_balance.png"), p8, width = 12)

  # ── Summary ───────────────────────────────────────────────────────────────
  cat(sprintf("\n── %s ───────────────────────────────────────────────────\n", chr_name))
  cat(sprintf("  Blocks          : %d\n",        nrow(stats_df)))
  cat(sprintf("  Total clusters  : %d\n",        sum(stats_df$n_clusters)))
  cat(sprintf("  Mean clusters/block : %.1f\n",  mean(stats_df$n_clusters)))
  cat(sprintf("  Mean entropy    : %.2f bits\n", mean(stats_df$shannon_entropy)))
  cat(sprintf("  Mean dominance  : %.2f\n",      mean(stats_df$dominance)))
  cat(sprintf("  Mean singleton rate: %.2f\n",   mean(stats_df$singleton_rate)))
}

# ── Iterate over all chromosomes present in the data ─────────────────────────
chromosomes <- levels(droplevels(stats_all$chr))
message("\nProcessing ", length(chromosomes), " chromosomes...")

for (chr_name in chromosomes) {
  message("\n[", chr_name, "]")

  stats_chr   <- stats_all   %>% filter(chr == chr_name)
  summary_chr <- summary_all %>% filter(chr == chr_name)

  if (nrow(stats_chr) == 0) {
    message("  No data — skipping")
    next
  }

  plot_chromosome(chr_name, summary_chr, stats_chr)
}

message("\nAll done. Plots saved under: ", file.path(OUT_DIR, "plots"))