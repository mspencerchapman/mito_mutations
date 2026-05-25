setwd("/Users/al28/Mounts/cancer/mtDNA/")

shearwater_panel_samples <- read.table(file = "shearwater_normal_panel_exome", sep = "\t", header = T, stringsAsFactors = F)

tissue_list <- list.files("Analysis/global_vaf_heatmaps_exome/")

for(z in 1:length(tissue_list)){
  if(file.exists(paste0("Analysis/global_vaf_heatmaps_exome/",tissue_list[z],"/",tissue_list[z],"_combined_germline_indels_merged.tsv"))){
    if(!(exists("all_tissue_germline"))){
      all_tissue_germline <- read.table(file = paste0("Analysis/global_vaf_heatmaps_exome/",tissue_list[z],"/",tissue_list[z],"_combined_germline_indels_merged.tsv"), header = T, stringsAsFactors = F, sep = "\t")
    } else{
      all_tissue_germline <- rbind(all_tissue_germline, read.table(file = paste0("Analysis/global_vaf_heatmaps_exome/",tissue_list[z],"/",tissue_list[z],"_combined_germline_indels_merged.tsv"), header = T, stringsAsFactors = F, sep = "\t"))
    }
  }
}

shearwater_background <- array(data = NA, dim = c(nrow(shearwater_panel_samples),16569,12))

for(i in 1:nrow(shearwater_panel_samples)){
  temp_counts <- read.table(file = paste0("Data/",shearwater_panel_samples[i,"tissue"],"/",shearwater_panel_samples[i,"sampleID"],"_MT_count.csv"), sep = ",", stringsAsFactors = F, header = T)
  germline_pos <- all_tissue_germline[which(all_tissue_germline$tissue == shearwater_panel_samples[i,"tissue"] & all_tissue_germline$donor == substr(shearwater_panel_samples[i,"sampleID"],1,7)),]
  germline_array <- rep(F, 16569)
  germline_array[germline_pos$pos] <- T
  
  alt_ref_sites <- c(263,8860,15326,750,4769,1438,16519)
  for(j in 1:length(alt_ref_sites)){
    germline_array[alt_ref_sites[j]] <- !(germline_array[alt_ref_sites[j]])
  }
  
  temp_counts[germline_array,] <- 0

    for(k in 1:12){
      shearwater_background[i,,k] <- temp_counts[,k]
    }
}

saveRDS(shearwater_background, file = "shearwater_background_exome.Rds")
