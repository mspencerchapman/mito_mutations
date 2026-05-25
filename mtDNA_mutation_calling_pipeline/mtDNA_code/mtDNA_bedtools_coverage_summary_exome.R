all_files <- list.files()
files <- all_files[grepl(pattern = "final",x = all_files)]

chrom_sizes <- read.table(file = "/lustre/scratch125/casm/team268im/al28/bed_ref/exome_chrom_sizes.csv", sep = ",", stringsAsFactors = F)
chrom_sizes[25,] <- c("MT",16569)

for(i in 1:length(files)){
  if(!(gsub(pattern = ".final", replacement = ".summary", x = files[i]) %in% all_files)){
    bedtools_cov <- read.table(files[i], sep = "\t", header = F)
    bedtools_cov$depth <- as.numeric(bedtools_cov$V2) * as.numeric(bedtools_cov$V3)
    summary_vec <- vector(length = 26)
    summary_vec[1] <- sum(bedtools_cov$depth[which(bedtools_cov$V1 %in% 1:22)]) / sum(as.numeric(chrom_sizes$V2[1:22]))
    for(j in 1:22){
      summary_vec[j + 1] <- sum(bedtools_cov$depth[which(bedtools_cov$V1 == j)]) / sum(as.numeric(chrom_sizes$V2[j]))
    }
    summary_vec[24] <- sum(bedtools_cov$depth[which(bedtools_cov$V1 == "X")]) / sum(as.numeric(chrom_sizes$V2[23]))
    summary_vec[25] <- sum(bedtools_cov$depth[which(bedtools_cov$V1 == "Y")]) / sum(as.numeric(chrom_sizes$V2[24]))
    summary_vec[26] <- sum(bedtools_cov$depth[which(bedtools_cov$V1 == "MT")]) / sum(as.numeric(chrom_sizes$V2[25]))
    write.table(summary_vec, file = gsub(pattern = ".final", replacement = ".summary", x = files[i]), sep = "\t", quote = F, row.names = F, col.names = F)
  }
}