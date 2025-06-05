rm(list=ls())
ls()

setwd("C:/Users/nnamd/Downloads")
options(stringsAsFactors = FALSE)

# =============================================================================
# Function to read vcf files
# =============================================================================

read.vcf <- function(file, special.char="##", ...) {
  my.search.term=paste0(special.char, ".*")
  all.lines=readLines(file)
  clean.lines=gsub(my.search.term, "",  all.lines)
  clean.lines=gsub("#CHROM", "CHROM", clean.lines)
  read.table(..., text=paste(clean.lines, collapse="\n"))
}

# =============================================================================
# Function to get derived counts
# =============================================================================

get.derived.counts <- function(my.vcf){
  ### Make an empty vector
  derived.counts = rep(0, nrow(my.vcf))
  
  ### Loop through the file and get the ALT allele
  ### count for every row (SNP site)
  for (i in 1:nrow(my.vcf)){
    test.line = as.vector(my.vcf[i,])
    x = unname(unlist((test.line[10:length(test.line)])))
    genotypes = unlist(lapply(strsplit(x, split=":"), `[[`, 1))
    alleles = unlist(strsplit(genotypes, "/"))
    derived.counts[i] = unname(table(alleles)["1"])
  }
  return(derived.counts)
}

# =============================================================================
# Function to get plots
# =============================================================================

get.plots = function(derived.counts, my.vcf, path_obs_afs, path_obs_n_exp_afs){
  ### count the total number of individuals in the pop
  N = ncol(my.vcf) - 9
  
  ### get the number of chromosomes
  n.chr = 2*N
  
  ### get a vector of empty bins
  possible.counts = seq(1, n.chr, by=1)
  
  ### calculate Waterson's theta
  Sn = nrow(my.vcf)
  theta.w = Sn/sum(1/(possible.counts[-length(possible.counts)]))
  
  ### Get the expected counts in each bin
  exp.counts = theta.w  * (1/possible.counts[-length(possible.counts)])
  
  ### Save results of observed AFS histogram
  png(file=path_obs_afs, width=14, height=8, units="in", res=500)
  histinfo = hist(derived.counts, breaks=possible.counts, plot=F)
  plot(histinfo, 
       col = "lightblue", 
       main = "Observed Allele Frequency Spectrum", 
       xlab = "Derived Allele Count", 
       ylab = "Frequency")
  dev.off()
  
  ### Plot the expected and the observed AFS
  ### side by side
  compare.counts = matrix(c(histinfo$counts, exp.counts), ncol=2, byrow=FALSE)
  colnames(compare.counts) = c("Obs", "Exp")
  
  png(file=path_obs_n_exp_afs, width=14, height=8, units="in", res=500)
  compare.plot = barplot(t(compare.counts), 
                         beside=T, 
                         col=c("lightblue", "black"),
                         main = "Observed vs Expected Allele Frequency Spectrum", 
                         xlab = "Allele Counts", 
                         ylab = "Frequency",
                         ylim = c(0, max(compare.counts) * 1.2))
  legend("topright", c("Obs", "Exp"), pch=c(15,15), col=c("lightblue", "black"))
  dev.off()
}


# =============================================================================
# Wild Sorghum WGS
# =============================================================================

ws.wgs.my.vcf = read.vcf("ws.wgs.merged.filtered.vcf", header=TRUE)

ws.wgs.derived.counts = get.derived.counts(ws.wgs.my.vcf)

get.plots(ws.wgs.derived.counts, 
          ws.wgs.my.vcf, 
          "G:/My Drive/PhD/project/structural_variant_genotyping_with_pangenome_graph/data_analysis/data_and_result/main/graph_minigraph-cactus/ws_wgs_obs_afs.png",
          "G:/My Drive/PhD/project/structural_variant_genotyping_with_pangenome_graph/data_analysis/data_and_result/main/graph_minigraph-cactus/ws_wgs_obs_n_exp_afs.png")

# =============================================================================
# TERRA RAW
# =============================================================================

tr.my.vcf = read.vcf("tr.merged.filtered.vcf", header=TRUE)

tr.derived.counts = get.derived.counts(tr.my.vcf)

get.plots(tr.derived.counts, 
          tr.my.vcf, 
          "G:/My Drive/PhD/project/structural_variant_genotyping_with_pangenome_graph/data_analysis/data_and_result/main/graph_minigraph-cactus/tr_obs_afs.png",
          "G:/My Drive/PhD/project/structural_variant_genotyping_with_pangenome_graph/data_analysis/data_and_result/main/graph_minigraph-cactus/tr_obs_n_exp_afs.png")

