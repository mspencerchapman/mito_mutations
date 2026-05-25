setwd("/Users/al28/Mounts/cancer/mtDNA/")

metanorm_list <- read.table("metanorm_masterlist_2020-06-18.txt", sep = "\t", stringsAsFactors = F, header = T)

sample_lists <- list.files("sample_lists/")
for(i in 1:length(sample_lists)){
  temp_list <- read.table(paste0("sample_lists/",sample_lists[i]), sep = "\t")
  colnames(temp_list) <- c("project","sample")
  temp_list$tissue <- substr(sample_lists[i],1,nchar(sample_lists[i])-4)
  if(i == 1){
    current_sample_list <- temp_list
  } else{
    current_sample_list <- rbind(current_sample_list,temp_list)
  }
}

new_samples <- metanorm_list[which(!(paste(metanorm_list$ID_INT_PROJECT,metanorm_list$SAMPLE_NAME,sep = "_") %in% paste(current_sample_list$project,current_sample_list$sample,sep = "_"))),]
