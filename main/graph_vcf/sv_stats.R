rm(list=ls())
ls()

library(tidyverse)
library(gtools)
library(ggpubr)

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

png(file="main/graph_vcf/sv/sv_stats_indels_per_chrom.png", width=14, height=8, units="in", res=500)

ggplot(df_long, aes(x = Chromosome, y = Count, fill = VariantType)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Insertions and Deletions Count per Chromosome",
       x = "Chromosome", y = "Count") +
  scale_fill_manual(values = c("Insertions" = "#1f77b4", "Deletions" = "#ff7f0e")) +
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 5000)) +
  theme_bw()

dev.off()


# =============================================================================
# Violin plot of indels length per chromosome
# =============================================================================

df = read.table("main/graph_vcf/sv/sv_stats_indels_len_per_chrom.txt", header=TRUE, sep="")

library(ggplot2)

df$Chromosome <- factor(df$Chromosome,
                             levels = paste0("chr", 1:10))

png(file="main/sv/sv_stats_indels_len_per_chrom.png", width=14, height=8, units="in", res=500)

ggplot(df, aes(x = Chromosome, y = log10(Length), fill = Type_of_Variant)) +
  geom_violin(position = position_dodge(width = 0.9), trim = FALSE, alpha = 0.7) +
  geom_boxplot(position = position_dodge(width = 0.9), width = 0.2, outlier.shape = NA) +
  labs(title = "Distribution of Variant Lengths by Chromosome",
       x = "Chromosome", y = "Length") +
  scale_fill_manual(values = c("Insertion" = "#1f77b4", "Deletion" = "#ff7f0e")) +
  theme_bw()

dev.off()

# =============================================================================

df_del <- subset(df, Type_of_Variant == "Deletion")
df_ins <- subset(df, Type_of_Variant == "Insertion")

del = ggplot(df_del, aes(x=Chromosome, y=log10(Length), fill=Type_of_Variant)) +
        geom_violin(trim=F) +
        geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
        labs(title = "Distribution of Deletion Lengths by Chromosome",
             x = "Chromosome", y = "log10(Length)") +
        scale_fill_manual(values = c("Deletion" = "#ff7f0e")) +
        theme_bw() +
        theme(legend.position = "none",
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
        scale_y_continuous(limits = c(0, 7), breaks = seq(0, 7, by = 1))

ins = ggplot(df_ins, aes(x=Chromosome, y=log10(Length), fill=Type_of_Variant)) +
        geom_violin(trim=F) +
        geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
        labs(title = "Distribution of Insertion Lengths by Chromosome",
             x = "Chromosome", y = "log10(Length)") +
        scale_fill_manual(values = c("Insertion" = "#1f77b4")) +
        theme_bw() +
        theme(legend.position = "none") +
        scale_y_continuous(limits = c(0, 7), breaks = seq(0, 7, by = 1))


png(file="main/sv/sv_stats_indels_len_per_chrom_2.png", width=14, height=8, units="in", res=500)

ggarrange(del, ins, 
          nrow=2,
          labels = "AUTO")

dev.off()