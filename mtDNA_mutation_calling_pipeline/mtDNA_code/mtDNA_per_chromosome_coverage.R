setwd("/Users/al28/Mounts/cancer/mtDNA/Analysis/bedtools_coverage/blood/")

test <- read.table("BMb.genomecov.txt", sep = "\t", stringsAsFactors = F, header = F)
test[,6] <- as.numeric(test[,2]) * as.numeric(test[,3])

chrom_sizes <- read.table("/Users/al28/Mounts/cancer/bed_ref/hg19_chrom_sizes.csv", sep = ",", stringsAsFactors = F, header = F)

per_chromosome_coverage <- array(data = NA, dim = c(25,2))
per_chromosome_coverage[,1] <- c(1:22,"X","Y","MT")

for(i in 1:nrow(per_chromosome_coverage)){
  if(i %in% 1:24){
    per_chromosome_coverage[i,2] <- sum(test[which(test[,1] == per_chromosome_coverage[i,1]),6]) / chrom_sizes[i,2]
  } else{
    per_chromosome_coverage[i,2] <- sum(test[which(test[,1] == per_chromosome_coverage[i,1]),6]) / 16569
  }
}