setwd("C:/Users/nnamd/Downloads")
options(stringsAsFactors = FALSE)

### Function to read vcf files
read.vcf <- function(file, special.char="##", ...) {
  my.search.term=paste0(special.char, ".*")
  all.lines=readLines(file)
  clean.lines=gsub(my.search.term, "",  all.lines)
  clean.lines=gsub("#CHROM", "CHROM", clean.lines)
  read.table(..., text=paste(clean.lines, collapse="\n"))
}

### Read in sample data
my.vcf = read.vcf("simonii-chr1.vcf", header=TRUE)

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
histinfo = hist(derived.counts, breaks=possible.counts)

### Plot the expected and the observed AFS
### side by side
compare.counts = matrix(c(histinfo$counts, exp.counts), ncol=2, byrow=FALSE)
colnames(compare.counts) = c("Obs", "Exp")
barplot(t(compare.counts), beside=T, col=c("lightblue", "black"))
legend("topright", c("Obs", "Exp"), pch=c(15,15), col=c("lightblue", "black"))
