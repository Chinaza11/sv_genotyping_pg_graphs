
# =============================================================================
# 
# =============================================================================

setwd("C:/Users/nnamd/Downloads")

df <- read.table("simulated.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(df) <- c("chrom1", "start", "chrom2", "end", "SV")  # Add more columns if your BED has more fields

df_sorted <- df[order(df$chrom1, df$start), ]

df_sorted$diff <- (df_sorted$end - df_sorted$start)

