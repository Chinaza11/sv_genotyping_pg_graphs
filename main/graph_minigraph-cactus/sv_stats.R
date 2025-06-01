rm(list=ls())
ls()

library(tidyverse)
library(gtools)
library(ggpubr)

setwd("G:/My Drive/PhD/project/structural_variant_genotyping_with_pangenome_graph/data_analysis/data_and_result/main/graph_minigraph-cactus")

wild_sorghum_wgs <- c("Grif16309","pi302118","pi302267","pi329250","pi329251","pi329252","pi369484","pi369487","pi524718","pi532564","pi532565","pi532566","pi532568","pi535995")

# =============================================================================
# Barplot of indels per chromosome
# =============================================================================

df_indels_per_chrom <- read.table("indels_per_chrom.txt", header=TRUE, sep="\t")
species <- wild_sorghum_wgs[1]

indels_count_barplot <- function(species){

  # Subset dataframe and keep only rows for species of interest
  df_species <- subset(df_indels_per_chrom, Species == species)[, c("Chromosome", "Insertions", "Deletions")]
  
  # Convert to long format for ggplot
  df_long <- df_species %>%
    pivot_longer(cols = c(Insertions, Deletions),
                 names_to = "VariantType",
                 values_to = "Count") 
  
  # Order table using the Chromosome column
  df_long <- df_long[mixedorder(df_long$Chromosome), ]
  
  df_long$Chromosome <- factor(df_long$Chromosome,
                               levels = paste0("chr", 1:10))
  
  plot = ggplot(df_long, aes(x = Chromosome, y = Count, fill = VariantType)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Insertions and Deletions Count per Chromosome",
         x = "Chromosome", y = "Count") +
    scale_fill_manual(values = c("Insertions" = "#1f77b4", "Deletions" = "#ff7f0e")) +
    scale_y_continuous(limits = c(0, 3000), breaks = seq(0, 3000, by = 500)) +
    theme_bw()
  
  return(plot)
}


for (i in wild_sorghum_wgs){
  print(i)
  plot <- indels_count_barplot(i)
  print(plot)
}

# =============================================================================
# Violin plot of indels length per chromosome
# =============================================================================

species <- wild_sorghum_wgs[1]

indels_len_violinplot <- function(species){
  
  file = paste0(species,".indels_len_per_chrom.txt")
  df = read.table(file, header=TRUE, sep="")
  
  df$Chromosome <- factor(df$Chromosome,
                               levels = paste0("chr", 1:10))
  
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
          scale_y_continuous(limits = c(0, 7.5), breaks = seq(0, 8, by = 1))
  
  ins = ggplot(df_ins, aes(x=Chromosome, y=log10(Length), fill=Type_of_Variant)) +
          geom_violin(trim=F) +
          geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
          labs(title = "Distribution of Insertion Lengths by Chromosome",
               x = "Chromosome", y = "log10(Length)") +
          scale_fill_manual(values = c("Insertion" = "#1f77b4")) +
          theme_bw() +
          theme(legend.position = "none") +
          scale_y_continuous(limits = c(0, 7.5), breaks = seq(0, 8, by = 1))
  
  plot = ggarrange(del, ins, 
            nrow=2,
            labels = "AUTO")
  
  return(plot)
}

indels_len_violinplot(species)
