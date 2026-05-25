args <- commandArgs(TRUE)

tissue <- as.character(args[1])

library("tidyverse")

setwd("/lustre/scratch125/casm/team268im/al28/mtDNA/")

exome_coverage <- read.table(file = "exome_coverage_pileup_and_bedtools_annotated.csv", sep = ",", header = T)

mtdna_ref <- readLines(con = "/lustre/scratch125/casm/team268im/al28/bed_ref/rCRS.fasta", warn = F)
mtdna_ref <- paste0(mtdna_ref[2:278], collapse = "")

tissue_list <- tissue

for(i in 1:length(tissue_list)){
  tissue_samples <- as.character(exome_coverage$Sample[which(exome_coverage$Tissue == tissue_list[i])])
  
  donors <- unique(substr(tissue_samples,1,7))
  
  for(j in 1:length(donors)){
    
    donor_samples <- tissue_samples[which(substr(tissue_samples,1,7) == donors[j])]
    
    for(k in 1:length(donor_samples)){
      if(file.exists(paste0("Data/",tissue_list[i],"/",donor_samples[k],"_MT_count.csv"))){
        sample_cov <- read.table(paste0("Data/",tissue_list[i],"/",donor_samples[k],"_MT_count.csv"), sep = ",", header = T)
        if(!(exists("donor_cov"))){
          donor_cov <- sample_cov
        } else{
          donor_cov <- donor_cov + sample_cov
        }
      }
    }
    
    donor_global_vaf <- (donor_cov[,1:6] + donor_cov[,7:12]) / rowSums(donor_cov)
    
    ref_base_array <- c("A","T","C","G","N")
    alt_base_array <- c("A","T","C","G","-","INS")
    
    for(x in 1:nrow(donor_cov)){
      ref_base <- substr(mtdna_ref,x,x)
      ref_base_pos <- which(ref_base_array == ref_base)
      for(y in 1:6){
        if(y != ref_base_pos){
          if(!(exists("germline_sites"))){
            germline_sites <- as.matrix(t(c(x,ref_base,alt_base_array[y],donor_global_vaf[x,y])))
          } else{
            germline_sites <- rbind(germline_sites,as.matrix(t(c(x,ref_base,alt_base_array[y],donor_global_vaf[x,y]))))
          }
        }
      }
    }
    
    germline_data <- as.data.frame(germline_sites)
    colnames(germline_data) <- c("pos","ref","mut","global_vaf")
    germline_data$global_vaf <- as.numeric(as.character(germline_data$global_vaf))
    germline_data$pos <- as.numeric(as.character(germline_data$pos))
    
    threshold <- min(max(c((mean(rowSums(donor_cov)) * 0.25) / (length(donor_samples)), 100)) / mean(rowSums(donor_cov)),0.25)
    
    ggplot(germline_data) +
      geom_histogram(aes(x = global_vaf),bins = 1000) +
      scale_x_log10() +
      geom_vline(xintercept = threshold, color = "red") +
      theme_bw()
    ggsave(paste0("/lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/global_vaf_heatmaps_exome/",tissue_list[i],"/",donors[j],"_global_vaf_histogram.pdf"), width = 15, height = 5)
    
    threshold_germline <- germline_data[which(germline_data$global_vaf > threshold),]
    
    sample_vaf_array <- array(data = NA, dim = c(length(donor_samples),nrow(threshold_germline)))
    excluded_samples <- vector()
    
    
    for(k in 1:length(donor_samples)){
      if(file.exists(paste0("Data/",tissue_list[i],"/",donor_samples[k],"_MT_count.csv"))){
        sample_cov <- read.table(paste0("Data/",tissue_list[i],"/",donor_samples[k],"_MT_count.csv"), sep = ",", header = T)
        for(x in 1:nrow(threshold_germline)){
          sample_vaf_array[k,x] <- sum(sample_cov[threshold_germline$pos[x],c(which(alt_base_array == threshold_germline$mut[x]),(which(alt_base_array == threshold_germline$mut[x]) + 6))]) / sum(sample_cov[threshold_germline$pos[x],]) 
        }
      } else{
        excluded_samples <- c(excluded_samples,donor_samples[k])
      }
    }
    rownames(sample_vaf_array) <- donor_samples
    colnames(sample_vaf_array) <- paste(threshold_germline$pos,threshold_germline$ref,threshold_germline$mut,sep = "_")
    
    sample_vaf_array <- sample_vaf_array[which(!(rownames(sample_vaf_array) %in% excluded_samples)), ,drop = F]
    
    if(nrow(sample_vaf_array) > 1){
      sample_vaf_array <- sample_vaf_array[hclust(dist(sample_vaf_array, method = "euclidean"))$order,]
      sample_vaf_array <- sample_vaf_array[,hclust(dist(t(sample_vaf_array), method = "euclidean"))$order]
    }
    
    write.table(x = sample_vaf_array, file = paste0("/lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/global_vaf_heatmaps_exome/",tissue_list[i],"/",donors[j],"_donor_vaf_table.tsv"), sep = "\t", row.names = T, col.names = T, quote = F)
    
    tidy_sample_vaf_array <- reshape2::melt(sample_vaf_array)
    colnames(tidy_sample_vaf_array) <- c("sample","variant","vaf")
    
    ggplot(tidy_sample_vaf_array) +
      geom_tile(aes(x = sample,y = variant, fill = vaf)) +
      theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1))
    ggsave(paste0("/lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/global_vaf_heatmaps_exome/",tissue_list[i],"/",donors[j],"_vaf_heatmap_global_vaf_threshold.pdf"), width = 40, height = 40)
    
    rm(germline_sites)
    rm(donor_cov)
  }
}