setwd("/Users/al28/Mounts/cancer/mtDNA/")

exome_cov_table <- read.csv(file = "exome_coverage_pileup_and_bedtools_annotated.csv", header = T)

new_study_sample_list <- read.table("sample_lists/psoriasis_exome.tsv", sep = "\t", stringsAsFactors = F, header = F)
new_study_sample_list <- new_study_sample_list[which(!(new_study_sample_list$V2 %in% exome_cov_table$Sample)),]
new_study_name <- "psoriasis_exome"

new_study_array <- as.data.frame(array(data = NA, dim = c(nrow(new_study_sample_list),ncol(exome_cov_table))))
colnames(new_study_array) <- colnames(exome_cov_table)

new_study_array$Tissue <- new_study_name
new_study_array$Study <- new_study_sample_list$V1
new_study_array$Sample <- new_study_sample_list$V2

seqx <- read.csv("/Users/al28/Desktop/canapps_2545.csv", header = T)
seqx <- rbind(seqx, read.csv("/Users/al28/Desktop/canapps_2636.csv", header = T))
seqx <- rbind(seqx, read.csv("/Users/al28/Desktop/canapps_2690.csv", header = T))

for(i in 1:nrow(new_study_array)){
  new_study_array$SeqX[i] <- seqx$Seq.X[which(seqx$Sample == new_study_array$Sample[i])] 
}

new_study_array$Normal <- ""

new_exome_cov_table <- rbind(exome_cov_table, new_study_array)

write.table(x = new_exome_cov_table, file = "exome_coverage_pileup_and_bedtools_annotated.csv", sep = ",", quote = F, col.names = T, row.names = F)
