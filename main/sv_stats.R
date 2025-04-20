rm(list=ls())
ls()

library(tidyverse)
library(gtools)

setwd("G:/My Drive/PhD/project/structural_variant_genotyping_with_pangenome_graph/data_analysis/data_and_result/")

# =============================================================================
# Barplot of indels per chromosome
# =============================================================================

df = read.table("main/sv/sv_stats_indels_per_chrom.txt", header=TRUE, sep="\t")

# Convert to long format for ggplot
df_long <- df %>%
  pivot_longer(cols = c(Insertions, Deletions),
               names_to = "VariantType",
               values_to = "Count") 

# Order table using the Chromosome column
df_long <- df_long[mixedorder(df_long$Chromosome), ]

df_long$Chromosome <- factor(df_long$Chromosome,
                             levels = paste0("chr", 1:10))

png(file="main/sv/sv_stats_indels_per_chrom.png", width=14, height=8, units="in", res=500)

ggplot(df_long, aes(x = Chromosome, y = Count, fill = VariantType)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Insertions and Deletions Count per Chromosome",
       x = "Chromosome", y = "Count") +
  scale_fill_manual(values = c("Insertions" = "#1f77b4", "Deletions" = "#ff7f0e")) +
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 5000)) +
  theme_bw()

dev.off()


# =============================================================================
# Violin plot of indels per chromosome
# =============================================================================

df = read.table("main/sv/sv_stats_indels_len_per_chrom.txt", header=TRUE, sep="")

library(ggplot2)

set.seed(42)
df <- data.frame(
  Chromosome = rep(c("chr1", "chr2", "chr3"), each = 20),
  Position = sample(100:10000, 60, replace = TRUE),
  Type_of_Variant = rep(c("Insertion", "Deletion"), 30),
  Length = c(
    sample(1:10, 30, replace = TRUE),  # Insertions
    sample(3:15, 30, replace = TRUE)   # Deletions
  )
)

# Violin plot
ggplot(df, aes(x = Chromosome, y = Length, fill = Type_of_Variant)) +
  geom_violin(position = position_dodge(width = 0.9), trim = FALSE, alpha = 0.7) +
  geom_boxplot(position = position_dodge(width = 0.9), width = 0.2, outlier.shape = NA) +
  labs(title = "Distribution of Variant Lengths by Chromosome",
       x = "Chromosome", y = "Length") +
  scale_fill_manual(values = c("Insertion" = "#1f77b4", "Deletion" = "#ff7f0e")) +
  theme_bw()
