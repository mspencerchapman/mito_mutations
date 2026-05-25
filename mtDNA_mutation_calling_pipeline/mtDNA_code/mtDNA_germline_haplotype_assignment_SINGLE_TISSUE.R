args <- commandArgs(TRUE)
tissue <- as.character(args[1])

setwd("/lustre/scratch125/casm/team268im/al28/mtDNA/")

library(tidyverse)

#Set up the reference haplotype array
haplotype_markers <- read.table("mtDNA_haplotype_markers_mitomap_2021-01-15.txt", sep = "\t", stringsAsFactors = F, header = T)
haplotype_markers$mut_ref <- paste(haplotype_markers$tpos,haplotype_markers$tnt,haplotype_markers$qnt,sep = "_")
haplotype_array <- as.data.frame(array(data = 0, dim = c(length(unique(haplotype_markers$haplo)),length(unique(haplotype_markers$mut_ref)))))
colnames(haplotype_array) <- unique(haplotype_markers$mut_ref)
rownames(haplotype_array) <- unique(haplotype_markers$haplo)
for(i in 1:nrow(haplotype_markers)){haplotype_array[haplotype_markers$haplo[i],haplotype_markers$mut_ref[i]] <- 1}

# uninformative_muts <- c("73_A_G","263_A_G","750_A_G","1438_A_G","2706_A_G","4769_A_G","7028_C_T","8860_A_G","11719_G_A","14766_C_T","15326_A_G","16519_T_C")
# haplotype_array <- haplotype_array[,-which(colnames(haplotype_array) %in% uninformative_muts)]

setwd("Analysis/global_vaf_heatmaps/")

##START THE ANALYSIS (IN THE COMPLETE SCRIPT THIS IS A OR LOOP)

germline_muts <- read.table(paste0(tissue,"/",tissue,"_combined_germline_indels_merged.tsv"), sep = "\t", stringsAsFactors = F, header = T)
germline_muts$mut[which(germline_muts$mut == ".")] <- "-"

split_rows <- NULL

for(i in 1:nrow(germline_muts)){
  if(nchar(germline_muts$ref[i]) > 1 & nchar(germline_muts$mut[i]) > 1){
    split_rows <- c(split_rows,i)
  }  
}

if(length(split_rows) > 1){
  for(i in 1:length(split_rows)){
    temp_array <- as.data.frame(array(data = NA, dim = c(nchar(germline_muts$ref[split_rows[i]]),ncol(germline_muts))))
    colnames(temp_array) <- colnames(germline_muts)
    temp_array$donor <- germline_muts$donor[split_rows[i]]
    temp_array$chr <- germline_muts$chr[split_rows[i]]
    temp_array$mut_id <- germline_muts$mut_id[split_rows[i]]
    temp_array$tissue <- germline_muts$tissue[split_rows[i]]
    for(j in 1:nchar(germline_muts$ref[split_rows[i]])){
      temp_array$pos[j] <- germline_muts$pos[split_rows[i]] + j - 1
      temp_array$ref[j] <- substr(germline_muts$ref[split_rows[i]],j,j)
      temp_array$mut[j] <- substr(germline_muts$mut[split_rows[i]],j,j)
    }
    germline_muts <- rbind(germline_muts, temp_array)
  }
  
  germline_muts <- germline_muts[-split_rows,]
}

germline_muts$mut_ref <- paste(germline_muts$pos,germline_muts$ref,germline_muts$mut,sep = "_")

germline_array <- as.data.frame(array(data = 0, dim = c(ncol(haplotype_array), length(unique(germline_muts$donor)))))

rownames(germline_array) <- unique(colnames(haplotype_array))
colnames(germline_array) <- unique(germline_muts$donor)

marker_muts <- germline_muts[which(germline_muts$mut_ref %in% colnames(haplotype_array)),]

for(i in 1:nrow(marker_muts)){
  germline_array[marker_muts$mut_ref[i],marker_muts$donor[i]] <- 1
}

cross_prod_array <- as.data.frame(array(data = NA, dim = c(nrow(haplotype_array),ncol(germline_array))))
colnames(cross_prod_array) <- colnames(germline_array)
rownames(cross_prod_array) <- rownames(haplotype_array)

haplotype_sum <- rowSums(haplotype_array)
donor_sum <- colSums(germline_array)

for(i in 1:nrow(cross_prod_array)){
  for(j in 1:ncol(cross_prod_array)){
    cross_prod_array[i,j] <- (3 * (as.numeric(haplotype_array[i,]) %*% as.numeric(germline_array[,j]))) - haplotype_sum[i] - donor_sum[j]
  }
}

temp_hap <- as.data.frame(array(data = NA, dim = c(ncol(cross_prod_array),11))) 
colnames(temp_hap) <- c("tissue","donor","total_germline_muts","total_marker_muts","top_hap","top_hap_collapsed","top_hap_score","top_hap_max_score","top_hap_present_muts","within_2","within_2_collapsed")

temp_hap$tissue <- tissue

for(i in 1:ncol(cross_prod_array)){
  temp_hap$donor[i] <- colnames(cross_prod_array)[i]
  temp_hap$total_germline_muts[i] <- length(which(germline_muts$donor == temp_hap$donor[i]))
  temp_hap$total_marker_muts[i] <- length(which(marker_muts$donor == temp_hap$donor[i]))
  
  cross_prod_slice <- cross_prod_array[,i]
  names(cross_prod_slice) <- rownames(cross_prod_array)
  cross_prod_slice <- sort(cross_prod_slice, decreasing = T)
  
  temp_top_haps <- names(cross_prod_slice)[which(cross_prod_slice == max(cross_prod_slice))]
  temp_within_2 <- names(cross_prod_slice)[which(cross_prod_slice >= (max(cross_prod_slice) - 2))]
  
  temp_hap$top_hap[i] <- paste(temp_top_haps, collapse = ",")
  if(length(temp_top_haps) == 1){
    temp_hap$top_hap_collapsed[i] <- temp_top_haps
  }else if(length(unique(str_replace_all(string = temp_top_haps, pattern = "[a-z]", replacement = ''))) == 1){
    temp_hap$top_hap_collapsed[i] <- unique(str_replace_all(string = temp_top_haps, pattern = "[a-z]", replacement = ''))
  } else if(length(unique(str_replace_all(string = unique(str_replace_all(string = temp_top_haps, pattern = "[a-z]", replacement = '')), pattern = "[0-9]", replacement = ''))) == 1){
    temp_hap$top_hap_collapsed[i] <- unique(str_replace_all(string = unique(str_replace_all(string = temp_top_haps, pattern = "[a-z]", replacement = '')), pattern = "[0-9]", replacement = ''))
  }
  
  temp_hap$top_hap_score[i] <- max(cross_prod_slice)
  
  temp_top_hap_maxs <- vector(length = length(temp_top_haps))
  for(j in 1:length(temp_top_hap_maxs)){
    temp_top_hap_maxs[j] <- length(which(haplotype_markers$haplo == temp_top_haps[j]))
  }
  
  temp_hap$top_hap_max_score[i] <- paste(temp_top_hap_maxs,collapse = ",")
  
  temp_top_hap_present_scores <- vector(length = length(temp_top_haps))
  for(j in 1:length(temp_top_hap_present_scores)){
    temp_top_hap_present_scores[j] <- as.numeric(haplotype_array[which(rownames(haplotype_array) == temp_top_haps[j]),]) %*% as.numeric(germline_array[,which(colnames(germline_array) == temp_hap$donor[i])])
  }
  
  temp_hap$top_hap_present_muts[i] <- paste(temp_top_hap_present_scores, collapse = ",")
  
  temp_hap$within_2[i] <- paste(temp_within_2, collapse = ",")
  if(length(temp_within_2) == 1){
    temp_hap$within_2_collapsed[i] <- temp_within_2
  }else if(length(unique(str_replace_all(string = temp_within_2, pattern = "[a-z]", replacement = ''))) == 1){
    temp_hap$within_2_collapsed[i] <- unique(str_replace_all(string = temp_within_2, pattern = "[a-z]", replacement = ''))
  } else if(length(unique(str_replace_all(string = unique(str_replace_all(string = temp_within_2, pattern = "[a-z]", replacement = '')), pattern = "[0-9]", replacement = ''))) == 1){
    temp_hap$within_2_collapsed[i] <- unique(str_replace_all(string = unique(str_replace_all(string = temp_within_2, pattern = "[a-z]", replacement = '')), pattern = "[0-9]", replacement = ''))
  }
}

write.table(x = temp_hap, paste0(tissue,"/",tissue,"_combined_haplotype_assignment.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

