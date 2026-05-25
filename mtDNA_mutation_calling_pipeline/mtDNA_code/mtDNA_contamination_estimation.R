setwd("/Users/al28/Mounts/cancer/mtDNA/Analysis/global_vaf_heatmaps/")

tissues <- list.files(path = ".")
tissues <- tissues[!(grepl(tissues, pattern = "PCAWG"))]

excluded_samples <- read.table("/Users/al28/Mounts/cancer/mtDNA/excluded_samples.csv", sep = ",", stringsAsFactors = F, header = T)

for(i in 1:length(tissues)){
  setwd(tissues[i])
  
  tissue_germline_calls <- read.table(paste0(tissues[i],"_combined_germline_indels_not_merged.tsv"), sep = "\t", header = T, stringsAsFactors = F)
  
  donor_vafs <- list.files(path = ".")
  donor_vafs <- donor_vafs[grepl(donor_vafs, pattern = "_donor_vaf_table.tsv")]
  for(j in 1:length(donor_vafs)){
    donor_id <- gsub(pattern = "_donor_vaf_table.tsv", replacement = "", x = donor_vafs[j])
    donor_germline_muts <- tissue_germline_calls[which(tissue_germline_calls$donor == donor_id),]
    
    if(nrow(donor_germline_muts) > 0){
      temp_vaf_table <- read.table(donor_vafs[j], sep = "\t")
      
      temp_excluded_samples <- which(row.names(temp_vaf_table) %in% excluded_samples$sample)
      
      temp_germline_table <- temp_vaf_table[,which(colnames(temp_vaf_table) %in% paste0("X",donor_germline_muts$pos,"_",donor_germline_muts$ref,"_",donor_germline_muts$mut))]
      
      temp_ave_diff_table <- array(data = NA, dim = c(nrow(temp_germline_table),ncol(temp_germline_table)))
      temp_contamination <- as.data.frame(array(data = NA, dim = c(nrow(temp_germline_table),9)))
      
      colnames(temp_contamination) <- c("tissue","donor","sampleID","sum_all_germline_diffs","sum_exc_top_bottom","min_diff","second_min_diff","contamination","reversion_mut")
      temp_contamination$tissue <- tissues[i]
      temp_contamination$donor <- donor_id
      temp_contamination$sampleID <- rownames(temp_germline_table)
      
      if(nrow(temp_contamination) > 1 & (length(temp_excluded_samples) + 1) < nrow(temp_germline_table)){
        for(x in 1:nrow(temp_germline_table)){
          for(y in 1:ncol(temp_germline_table)){
            if(length(excluded_samples) < 1){
              temp_ave_diff_table[x,y] <- temp_germline_table[x,y] - ((sum(temp_germline_table[,y]) - temp_germline_table[x,y]) / (nrow(temp_germline_table) - 1))
            } else if(!(x %in% temp_excluded_samples)){
              temp_ave_diff_table[x,y] <- temp_germline_table[x,y] - ((sum(temp_germline_table[,y]) - temp_germline_table[x,y] - sum(temp_germline_table[temp_excluded_samples,y])) / (nrow(temp_germline_table) - 1 - length(temp_excluded_samples)))
            } else{
              temp_ave_diff_table[x,y] <- temp_germline_table[x,y] - ((sum(temp_germline_table[,y]) - sum(temp_germline_table[temp_excluded_samples,y])) / (nrow(temp_germline_table) - length(temp_excluded_samples)))
            }
          }
          
          temp_contamination$sum_all_germline_diffs[x] <- sum(temp_ave_diff_table[x,])
          temp_contamination$sum_exc_top_bottom[x] <- sum(temp_ave_diff_table[x,]) - min(temp_ave_diff_table[x,]) - max(temp_ave_diff_table[x,])
          temp_contamination$min_diff[x] <- min(temp_ave_diff_table[x,])
          temp_contamination$second_min_diff[x] <- sort(temp_ave_diff_table[x,])[2]
          
          if(temp_contamination$sum_exc_top_bottom[x] <= -0.05 & temp_contamination$second_min_diff[x] <= -0.02){
            temp_contamination$contamination[x] <- "YY"
          } else if(temp_contamination$sum_exc_top_bottom[x] <= -0.05 & temp_contamination$second_min_diff[x] <= -0.01){
            temp_contamination$contamination[x] <- "YX"
          } 
          else if(temp_contamination$sum_exc_top_bottom[x] <= -0.05){
            temp_contamination$contamination[x] <- "YN"
          } else if(temp_contamination$second_min_diff[x] <= -0.02){
            temp_contamination$contamination[x] <- "NY"
          } else if(temp_contamination$second_min_diff[x] <= -0.01){
            temp_contamination$contamination[x] <- "NX"
          } else{
            temp_contamination$contamination[x] <- "NN"
          }
          
          if((temp_contamination$sum_exc_top_bottom[x] - temp_contamination$sum_all_germline_diffs[x]) >= 0.05){
            temp_contamination$reversion_mut[x] <- "Y"
          } else{
            temp_contamination$reversion_mut[x] <- ""
          }
        }
        
        if(i == 1 & j == 1){
          contamination_results  <- temp_contamination 
        } else{
          contamination_results <- rbind(contamination_results, temp_contamination)
        }
      }
    }
  }
  setwd("..")
}

write.table(x = contamination_results, file = "/Users/al28/Mounts/cancer/mtDNA/mtDNA_contamination_estimates.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
