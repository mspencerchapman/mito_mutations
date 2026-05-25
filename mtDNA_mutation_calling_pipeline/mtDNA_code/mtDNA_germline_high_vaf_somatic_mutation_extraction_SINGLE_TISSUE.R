args <- commandArgs(TRUE)
tissue <- as.character(args[1])

setwd("/lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/global_vaf_heatmaps/")
excluded_samples <- read.table("/lustre/scratch125/casm/team268im/al28/mtDNA/excluded_samples.csv", sep = ",", stringsAsFactors = F, header = T)

tissue_list <- list.files(".")
tissue_list <- tissue_list[grep(pattern = "PCAWG",x = tissue_list, invert = T)]

z<-which(tissue_list==tissue)

library(tidyverse)

##NOW RUN  IT USING THIS INDEX
if(tissue_list[z] == "blood"){
  donors <- "TG01"
} else if(tissue_list[z] == "colon_ibd"){
  donors <- unique(gsub("_.*","",list.files(tissue_list[z]))[which(gsub("_.*","",list.files(tissue_list[z])) != "colon")])
} else if(tissue_list[z] == "colon"){
  donors <- unique(gsub("_.*","",list.files(tissue_list[z]))[which(gsub("_.*","",list.files(tissue_list[z])) != "colon")])
  donors <- donors[which(!(donors %in% c("PD36813","PD36814","PD37226","PD37265","PD37457")))]
} else if(tissue_list[z] == "immune"){
  donors <- unique(gsub("_.*","",list.files(tissue_list[z]))[which(gsub("_.*","",list.files(tissue_list[z])) != "immune")])
} else if(tissue_list[z] == "blood_MPN") {
  library(tidyverse)
  donors<-data.frame(ID=tissue_samples) %>%
    separate(ID,into = c("text", "num"),sep = "(?<=[A-Za-z])(?=[0-9])")%>%
    separate(num,into=c("num","text2"),sep= "(?<=[0-9])(?=[A-Za-z])")%>%
    unite(col="new_PD",text,num,sep = "")%>%
    pull(new_PD)%>%
    unique()
  
} else {
  donors <- unique(substr(list.files(tissue_list[z]),1,7)[which(substr(list.files(tissue_list[z]),1,7) != substr(paste(tissue_list[z],"combined_",sep = "_"),1,7))])
}

for(y in 1:length(donors)){
  germline_variants <- vector()
  somatic_variants <- vector()
  
  high_vaf_calls <- read.table(paste0(tissue_list[z],"/",donors[y],"_donor_vaf_table.tsv"))
  high_vaf_calls <- high_vaf_calls[which(!(rownames(high_vaf_calls) %in% excluded_samples$sample[which(excluded_samples$tissue == tissue_list[z] & excluded_samples$donor == donors[y])])),]
  
  if(nrow(high_vaf_calls) > 0){
    
    for(j in 1:ncol(high_vaf_calls)){
      if((length(which(high_vaf_calls[,j] > 0.9)) >= (floor(0.9 * nrow(high_vaf_calls))) & nrow(high_vaf_calls) > 2) | (length(which(high_vaf_calls[,j] > 0.9)) == nrow(high_vaf_calls) & nrow(high_vaf_calls) <= 2)){
        germline_variants <- c(germline_variants, j)
      } else if(length(which(high_vaf_calls[,j] > 0.1)) <= (floor(0.9 * nrow(high_vaf_calls))) & length(which(high_vaf_calls[,j] > 0.1)) >= 1 & length(which(high_vaf_calls[,j] < 0.01)) >= (nrow(high_vaf_calls) - (floor(0.9 * nrow(high_vaf_calls)))) & length(which(high_vaf_calls[,j] < 0.01)) >= 1){
        somatic_variants <- c(somatic_variants, j)
      }
    }
    
    if(length(germline_variants) > 0){
      germline_list <- strsplit(colnames(high_vaf_calls)[germline_variants], split = "_")
      germline_array <- data.frame(matrix(unlist(germline_list), nrow=length(germline_list), byrow=T), stringsAsFactors = F)
      colnames(germline_array) <- c("pos","ref","mut")
      
      for(k in 1:nrow(germline_array)){
        germline_array$pos[k] <- substr(germline_array$pos[k],2,nchar(germline_array$pos[k]))
      }
      
      germline_array$chr <- "MT"
      germline_array$donor <- donors[y]
      germline_array$pos <- as.numeric(germline_array$pos)
      
      if(!(exists("combined_germline_array"))){
        combined_germline_array <- germline_array
      } else{
        combined_germline_array <- rbind(combined_germline_array, germline_array)
      }
    }
    
    if(length(somatic_variants) > 0){
      for(k in 1:length(somatic_variants)){
        somatic_sample_list <- rownames(high_vaf_calls)[which(high_vaf_calls[,somatic_variants[k]] > 0.1)]
        somatic_list <- strsplit(colnames(high_vaf_calls)[somatic_variants[k]], split = "_")
        somatic_temp_array <- data.frame(matrix(unlist(somatic_list), nrow=length(somatic_list), byrow=T), stringsAsFactors = F)
        if(length(somatic_sample_list) > 1){
          for(x in 2:length(somatic_sample_list)){
            somatic_temp_array <- rbind(somatic_temp_array,somatic_temp_array[1,])
          }
        }
        colnames(somatic_temp_array) <- c("pos","ref","mut")
        somatic_temp_array$chr <- "MT"
        somatic_temp_array$donor <- donors[y]
        somatic_temp_array$sample <- somatic_sample_list
        somatic_temp_array$vaf <- high_vaf_calls[which(high_vaf_calls[,somatic_variants[k]] > 0.1),somatic_variants[k]]
        
        if(k == 1){
          somatic_array <- somatic_temp_array
        } else{
          somatic_array <- rbind(somatic_array, somatic_temp_array)
        }
      }
      
      for(k in 1:nrow(somatic_array)){
        somatic_array$pos[k] <- substr(somatic_array$pos[k],2,nchar(somatic_array$pos[k]))
      }
      somatic_array$pos <- as.numeric(somatic_array$pos)
      
      if(!(exists("combined_somatic_array"))){
        combined_somatic_array <- somatic_array
      } else{
        combined_somatic_array <- rbind(combined_somatic_array, somatic_array)
      }
    }
  }
}

if(exists("combined_germline_array")){
  combined_germline_array <- combined_germline_array[c("donor","chr","pos","ref","mut")]
  combined_germline_array$mut_id <- paste(combined_germline_array$donor,combined_germline_array$chr,combined_germline_array$pos,combined_germline_array$ref,combined_germline_array$mut, sep = "_")
  combined_germline_array$tissue <- tissue_list[z]
  combined_germline_array <- combined_germline_array[order(combined_germline_array$donor,combined_germline_array$pos),]
  
  write.table(x = combined_germline_array, file = paste0(tissue_list[z],"/",tissue_list[z],"_combined_germline_indels_not_merged.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  for(i in 2:nrow(combined_germline_array)){
    if(combined_germline_array$donor[i] == combined_germline_array$donor[i - 1] &
       combined_germline_array$pos[i] == combined_germline_array$pos[i - 1] + 1){
      if(!(exists("germline_merger"))){
        germline_merger <- as.data.frame(t(c(1, i)))
      }else if((i - 1) %in% germline_merger$V2){
        germline_merger <- rbind(germline_merger, t(c(max(germline_merger$V1), i)))
      }else{
        germline_merger <- rbind(germline_merger, t(c(max(germline_merger$V1) + 1, i)))
      }
    }
  }
  
  if(exists("germline_merger")){
    for(i in 1:max(germline_merger$V1)){
      min_row <- min(germline_merger$V2[which(germline_merger$V1 == i)]) - 1
      max_row <- max(germline_merger$V2[which(germline_merger$V1 == i)])
      merged_ref_str <- paste0(combined_germline_array$ref[min_row:max_row], collapse = "")
      if(all(combined_germline_array$mut[min_row:max_row] == ".")){
        merged_mut_str <- "."
      } else if(any(combined_germline_array$mut[min_row:max_row] %in% c(".","INS"))){
        merged_mut_str <- "COMPLEX"
      } else{
        merged_mut_str <- paste0(combined_germline_array$mut[min_row:max_row], collapse = "")
      }
      merged_mut_id <- paste(combined_germline_array$donor[min_row],combined_germline_array$chr[min_row],combined_germline_array$pos[min_row],merged_ref_str,merged_mut_str,sep = "_")
      combined_germline_array[min_row,] <- c(combined_germline_array$donor[min_row],combined_germline_array$chr[min_row],combined_germline_array$pos[min_row],merged_ref_str,merged_mut_str,merged_mut_id,combined_germline_array$tissue[min_row])
      combined_germline_array <- combined_germline_array[-((min_row + 1):max_row),]
      germline_merger$V2 <- germline_merger$V2 - (max_row - min_row)
      
      combined_germline_array$pos <- as.numeric(combined_germline_array$pos)
    }
  }
  
  write.table(x = combined_germline_array, file = paste0(tissue_list[z],"/",tissue_list[z],"_combined_germline_indels_merged.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  rm(combined_germline_array)
}

if(exists("combined_somatic_array")){
  combined_somatic_array <- combined_somatic_array[c("sample","chr","pos","ref","mut","vaf","donor")]
  combined_somatic_array$mut_id <- paste(combined_somatic_array$donor,combined_somatic_array$chr,combined_somatic_array$pos,combined_somatic_array$ref,combined_somatic_array$mut, sep = "_")
  combined_somatic_array$tissue <- tissue_list[z]
  combined_somatic_array <- combined_somatic_array[order(combined_somatic_array$donor,combined_somatic_array$sample,combined_somatic_array$pos),]
  
  write.table(x = combined_somatic_array, file = paste0(tissue_list[z],"/",tissue_list[z],"_combined_somatic_indels_not_merged.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  for(i in 2:nrow(combined_somatic_array)){
    if(combined_somatic_array$donor[i] == combined_somatic_array$donor[i - 1] &
       combined_somatic_array$sample[i] == combined_somatic_array$sample[i - 1] &
       combined_somatic_array$pos[i] == combined_somatic_array$pos[i - 1] + 1){
      if(!(exists("somatic_merger"))){
        somatic_merger <- as.data.frame(t(c(1, i)))
      }else if((i - 1) %in% somatic_merger$V2){
        somatic_merger <- rbind(somatic_merger, t(c(max(somatic_merger$V1), i)))
      }else{
        somatic_merger <- rbind(somatic_merger, t(c(max(somatic_merger$V1) + 1, i)))
      }
    }
  }
  
  if(exists("somatic_merger")){
    for(i in 1:max(somatic_merger$V1)){
      min_row <- min(somatic_merger$V2[which(somatic_merger$V1 == i)]) - 1
      max_row <- max(somatic_merger$V2[which(somatic_merger$V1 == i)])
      merged_ref_str <- paste0(combined_somatic_array$ref[min_row:max_row], collapse = "")
      if(all(combined_somatic_array$mut[min_row:max_row] == ".")){
        merged_mut_str <- "."
      } else if(any(combined_somatic_array$mut[min_row:max_row] %in% c(".","INS"))){
        merged_mut_str <- "COMPLEX"
      } else{
        merged_mut_str <- paste0(combined_somatic_array$mut[min_row:max_row], collapse = "")
      }
      merged_vaf <- mean(combined_somatic_array$vaf[min_row:max_row][which(combined_somatic_array$mut[min_row:max_row] != "INS")])
      merged_mut_id <- paste(combined_somatic_array$donor[min_row],combined_somatic_array$chr[min_row],combined_somatic_array$pos[min_row],merged_ref_str,merged_mut_str,sep = "_")
      combined_somatic_array[min_row,] <- c(combined_somatic_array$sample[min_row],combined_somatic_array$chr[min_row],combined_somatic_array$pos[min_row],merged_ref_str,merged_mut_str,merged_vaf,combined_somatic_array$donor[min_row],merged_mut_id,combined_somatic_array$tissue[min_row])
      combined_somatic_array <- combined_somatic_array[-((min_row + 1):max_row),]
      somatic_merger$V2 <- somatic_merger$V2 - (max_row - min_row)
      
      combined_somatic_array$pos <- as.numeric(combined_somatic_array$pos)
      combined_somatic_array$vaf <- as.numeric(combined_somatic_array$vaf)
    }
  }
  
  write.table(x = combined_somatic_array, file = paste0(tissue_list[z],"/",tissue_list[z],"_combined_somatic_indels_merged.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  rm(combined_somatic_array)
}

if(exists("somatic_merger")){
  rm(somatic_merger)
}
if(exists("germline_merger")){
  rm(germline_merger)
}






