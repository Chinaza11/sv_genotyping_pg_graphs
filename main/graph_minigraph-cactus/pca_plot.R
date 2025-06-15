rm(list=ls())
ls()

library(ggplot2)
library(readxl)

setwd("G:/My Drive/PhD/project/structural_variant_genotyping_with_pangenome_graph/data_analysis/data_and_result/main/graph_minigraph-cactus/terra_raw")

# =============================================================================
# 
# =============================================================================

plot_pcs <- function(eigenvec_file, eigenval_file, pc_x = "PC1", pc_y = "PC2", pcs_file_path, save_percent_var_plot = F, percent_var_file_path) {
  
  # read eigenvalues and compute percent variance explained
  eigenval <- read.table(eigenval_file, header=F)
  colnames(eigenval) <- "eigenval"                              
  eigenval$percent_var <- round((eigenval$eigenval / sum(eigenval$eigenval)) * 100, 1)
  eigenval$PC <- factor(paste0("PC", 1:nrow(eigenval)), levels = paste0("PC", 1:nrow(eigenval)))
  eigenval <- eigenval[, c("PC", "eigenval", "percent_var")]
  
  # read eigenvectors
  header_line <- readLines(eigenvec_file, n = 1)
  header <- strsplit(sub("^#", "", header_line), "\\s+")[[1]]
  eigenvec <- read.table(eigenvec_file, header = FALSE, skip = 1)
  colnames(eigenvec) <- header
  
  # merge eigenvectors file with genotype info file
  df_ids <- read.csv("df_ids.csv")
  eigenvec <- merge(eigenvec, df_ids, by.x = "IID", by.y = "X.IID")
  genotype_info <- read_excel("Songsomboon_table_revised1.xlsx", sheet = "Supplement table 1", skip = 2)
  eigenvec <- merge(eigenvec, genotype_info, by= "Genotype")
  
  # extract percent variance for selected PCs
  x_label <- paste0(pc_x, " (", eigenval$percent_var[eigenval$PC == pc_x], "%)")
  y_label <- paste0(pc_y, " (", eigenval$percent_var[eigenval$PC == pc_y], "%)")
  
  custom_colors <- c(
    "#FF0000",  # Red
    "#00CFFF",  # Cyan
    "#1D6F2B",  # Dark Green
    "#FF00FF",  # Magenta
    "#FFA500",  # Orange
    "#00FF00",  # Bright Green
    "#FFB6C1",  # Light Pink
    "#000080",  # Navy blue
    "#FFFF00"   # Yellow
  )
  
  categories <- c("Region", "Cluster")
  
  for(i in categories){
    pc_plot = ggplot(data = eigenvec, aes(x = !!sym(pc_x), y = !!sym(pc_y), color = !!sym(i), shape = Type)) +
      geom_point(size=2) +
      xlab(x_label) +
      ylab(y_label) +
      theme_minimal() +
      scale_color_manual(values = custom_colors)
      
    print(pc_plot)
    
    # save plot
    png(file=paste0(pcs_file_path, "_", i, ".png"), width=14, height=8, units="in", res=500)
    print(pc_plot)
    dev.off()
  }
  
  # percent variance explained plot
  if(save_percent_var_plot == T){
    percent_var_plot = ggplot(eigenval, aes(PC, percent_var)) + geom_bar(stat = "identity") + 
      ylab("Percent variance explained") + 
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
      theme_minimal()
    
    percent_var_plot
    
    png(file=percent_var_file_path, width=14, height=8, units="in", res=500)
    print(percent_var_plot)
    dev.off()
  }
}

# =============================================================================
# 
# =============================================================================

plot_pcs(eigenvec_file = "pca_no_filter/plink.eigenvec", 
         eigenval_file = "pca_no_filter/plink.eigenval",
         pcs_file_path = "pca_no_filter/no_filter_pca1_n_pca2",
         save_percent_var_plot = T,
         percent_var_file_path = "pca_no_filter/no_filter_percent_var.png")

plot_pcs(eigenvec_file = "pca_no_filter/plink.eigenvec", 
         eigenval_file = "pca_no_filter/plink.eigenval", 
         pc_x = "PC3", 
         pc_y = "PC4", 
         pcs_file_path = "pca_no_filter/no_filter_pca3_n_pca4")

# =============================================================================
# 
# =============================================================================

plot_pcs(eigenvec_file = "pca_minus_maf_0.05/plink.eigenvec", 
         eigenval_file = "pca_minus_maf_0.05/plink.eigenval",
         pcs_file_path = "pca_minus_maf_0.05/filtered_pca1_n_pca2",
         save_percent_var_plot = T,
         percent_var_file_path = "pca_minus_maf_0.05/filtered_percent_var.png")

plot_pcs(eigenvec_file = "pca_minus_maf_0.05/plink.eigenvec", 
         eigenval_file = "pca_minus_maf_0.05/plink.eigenval", 
         pc_x = "PC3", 
         pc_y = "PC4", 
         pcs_file_path = "pca_minus_maf_0.05/filtered_pca3_n_pca4")

# =============================================================================
# 
# =============================================================================


plot_pcs(eigenvec_file = "pca_minus_maf_0.05_minus_ld/plink.eigenvec", 
         eigenval_file = "pca_minus_maf_0.05_minus_ld/plink.eigenval",
         pcs_file_path = "pca_minus_maf_0.05_minus_ld/filtered_pca1_n_pca2",
         save_percent_var_plot = T,
         percent_var_file_path = "pca_minus_maf_0.05_minus_ld/filtered_percent_var.png")

plot_pcs(eigenvec_file = "pca_minus_maf_0.05_minus_ld/plink.eigenvec", 
         eigenval_file = "pca_minus_maf_0.05_minus_ld/plink.eigenval", 
         pc_x = "PC3", 
         pc_y = "PC4", 
         pcs_file_path = "pca_minus_maf_0.05_minus_ld/filtered_pca3_n_pca4")

# =============================================================================
# 
# =============================================================================