# =========================

# Haploblock QC Analysis

# =========================

library(ggplot2)
library(dplyr)
library(readr)
library(scales)

# -------------------------

# Load data

# -------------------------

df <- read.table("combined_haploblocks_clean.tsv",
header = FALSE,
sep = "\t",
stringsAsFactors = FALSE)

colnames(df) <- c("CHR", "START", "END", "SIZE")


# -------------------------
# Basic transforms
# -------------------------

df$CHR <- factor(df$CHR,
                 levels = unique(df$CHR[order(as.numeric(gsub("chr", "", df$CHR)))]))

df <- df %>%
  mutate(
    LOG_SIZE = log10(SIZE)
  )


dir.create("figures", showWarnings = FALSE)

# =========================

# 1. Global size distribution

# =========================

p1 <- ggplot(df, aes(LOG_SIZE)) +
geom_histogram(bins = 80, fill = "steelblue") +
theme_minimal() +
labs(title = "Haploblock Size Distribution",
x = "log10(Haploblock size bp)",
y = "Count")

ggsave("figures/01_size_distribution.pdf", p1, width = 8, height = 5)

# =========================

# 2. Chromosome boxplot

# =========================

p2 <- ggplot(df, aes(x = CHR, y = LOG_SIZE)) +
geom_boxplot(outlier.size = 0.5) +
theme_minimal() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
labs(title = "Haploblock Size per Chromosome",
x = "Chromosome",
y = "log10(Size bp)")

ggsave("figures/02_chr_boxplot.pdf", p2, width = 10, height = 5)

# =========================

# 3. ECDF (cumulative distribution)

# =========================

p3 <- ggplot(df, aes(SIZE)) +
stat_ecdf(geom = "step") +
scale_x_log10() +
theme_minimal() +
labs(title = "Cumulative Distribution of Haploblock Sizes",
x = "Haploblock size (bp, log scale)",
y = "Proportion")

ggsave("figures/03_ecdf.pdf", p3, width = 8, height = 5)

# =========================

# 4. Blocks per chromosome (density proxy)

# =========================

chr_counts <- df %>%
group_by(CHR) %>%
summarise(N_BLOCKS = n()) %>%
arrange(CHR)

p4 <- ggplot(chr_counts, aes(x = CHR, y = N_BLOCKS)) +
geom_bar(stat = "identity", fill = "darkgreen") +
theme_minimal() +
theme(axis.text.x = element_text(angle = 90)) +
labs(title = "Number of Haploblocks per Chromosome",
x = "Chromosome",
y = "Block count")

ggsave("figures/04_blocks_per_chr.pdf", p4, width = 10, height = 5)

# =========================

# 5. Mean size per chromosome

# =========================

chr_mean <- df %>%
group_by(CHR) %>%
summarise(MEAN_SIZE = mean(SIZE))

p5 <- ggplot(chr_mean, aes(x = CHR, y = MEAN_SIZE)) +
geom_bar(stat = "identity", fill = "purple") +
theme_minimal() +
theme(axis.text.x = element_text(angle = 90)) +
scale_y_log10() +
labs(title = "Mean Haploblock Size per Chromosome",
x = "Chromosome",
y = "Mean size (bp, log scale)")

ggsave("figures/05_mean_size_chr.pdf", p5, width = 10, height = 5)

# =========================

# 6. Size vs position (genome structure)

# =========================

p6 <- ggplot(df, aes(x = START, y = SIZE)) +
geom_point(alpha = 0.3, size = 0.5) +
facet_wrap(~CHR, scales = "free_x") +
scale_y_log10() +
theme_minimal() +
labs(title = "Haploblock Size vs Genomic Position",
x = "Position",
y = "Size (bp, log scale)")

ggsave("figures/06_size_vs_position.pdf", p6, width = 12, height = 8)

# =========================

# 7. Density (blocks per Mb)

# =========================

chr_density <- df %>%
group_by(CHR) %>%
summarise(DENSITY = n() / (max(END) / 1e6))

p7 <- ggplot(chr_density, aes(x = CHR, y = DENSITY)) +
geom_bar(stat = "identity", fill = "orange") +
theme_minimal() +
theme(axis.text.x = element_text(angle = 90)) +
labs(title = "Haploblock Density per Chromosome",
x = "Chromosome",
y = "Blocks per Mb")

ggsave("figures/07_density.pdf", p7, width = 10, height = 5)

# =========================

# Done

# =========================

print("All figures saved in /figures/")
