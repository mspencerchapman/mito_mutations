args <- commandArgs(TRUE)
new_study_name <- args[1]
canapps_dir<-args[2] #"/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/canapps_info/blood_CML"

setwd("/lustre/scratch125/casm/team268im/al28/mtDNA/")

genome_cov_table <- read.csv(file = "whole_genome_coverage_pileup_and_bedtools_annotated.csv", header = T)

new_study_sample_list <- read.table(paste0("sample_lists/",new_study_name,".tsv"), sep = "\t", stringsAsFactors = F, header = F)
new_study_sample_list <- new_study_sample_list[which(!(new_study_sample_list$V2 %in% genome_cov_table$Sample)),]

new_study_array <- as.data.frame(array(data = NA, dim = c(nrow(new_study_sample_list),ncol(genome_cov_table))))
colnames(new_study_array) <- colnames(genome_cov_table)

new_study_array$Tissue <- new_study_name
new_study_array$Study <- new_study_sample_list$V1
new_study_array$Sample <- new_study_sample_list$V2

###This is a manual step - the CanApps info must be downloaded and placed as a csv file in the 'canapps_dir'
canapps_files<-list.files(path=canapps_dir,pattern=".csv",full.names = T)
seqx<-Reduce(rbind,lapply(canapps_files,function(file) {read.csv(file,header=T)}))
for(i in 1:nrow(new_study_array)){
  new_study_array$SeqX[i] <- seqx$Seq.X[which(seqx$Sample == new_study_array$Sample[i])] 
}

new_study_array$Normal <- ""

new_genome_cov_table <- rbind(genome_cov_table, new_study_array)

write.table(x = new_genome_cov_table, file = "whole_genome_coverage_pileup_and_bedtools_annotated.csv", sep = ",", quote = F, col.names = T, row.names = F)
