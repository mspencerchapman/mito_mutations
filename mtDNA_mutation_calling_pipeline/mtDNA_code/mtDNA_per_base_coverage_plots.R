library("tidyverse")

setwd("/Users/al28/Mounts/cancer/mtDNA/")

wgs_coverage <- read.table(file = "whole_genome_coverage.csv", sep = ",", header = T)
wgs_coverage$mtDNAcoverage <- NA

tissue_list <- as.character(unique(wgs_coverage$Tissue))

for(i in 1:length(tissue_list)){
  tissue_samples <- as.character(wgs_coverage$Sample[which(wgs_coverage$Tissue == tissue_list[i])])
  per_base_coverage <- as.data.frame(array(data = NA, dim = c(16569,length(tissue_samples))))
  colnames(per_base_coverage) <- tissue_samples
  
  average_coverage <- vector(length = length(tissue_samples))
  names(average_coverage) <- tissue_samples
  
  for(j in 1:length(tissue_samples)){
    if(file.exists(paste0("Data/",tissue_list[i],"/",tissue_samples[j],"_MT_count.csv"))){
      sample_cov <- read.table(paste0("Data/",tissue_list[i],"/",tissue_samples[j],"_MT_count.csv"), sep = ",", header = T)
      per_base_coverage[,j] <- rowSums(sample_cov)
      average_coverage[j] <- mean(per_base_coverage[,j])
    }
  }
  
  per_base_cov_norm <- as.data.frame(array(data = NA, dim = c(16569,length(tissue_samples))))
  colnames(per_base_cov_norm) <- tissue_samples
  for(j in 1:length(tissue_samples)){
    per_base_cov_norm[,j] <- per_base_coverage[,j] / average_coverage[j]
  }
  
  if(!(dir.exists(paste0("Analysis/per_sample_norm_cov_plots/",tissue_list[i])))){
    dir.create(paste0("Analysis/per_sample_norm_cov_plots/",tissue_list[i]))
  }
  
  for(j in 1:length(tissue_samples)){
    per_base_cov_norm_tidy <- as.data.frame(array(data = NA, dim = c(16569,2)))
    colnames(per_base_cov_norm_tidy) <- c("pos","norm_cov")
    per_base_cov_norm_tidy$pos <- 1:16569
    per_base_cov_norm_tidy$norm_cov <- per_base_cov_norm[,j]
    
    norm_cov_plot <- ggplot(data = per_base_cov_norm_tidy) +
      geom_point(aes(x = pos, y = norm_cov), size = 0.2) +
      scale_y_continuous(limits = c(0, 1.1 * max(per_base_cov_norm_tidy$norm_cov))) +
      theme_bw()
    ggsave(paste0("Analysis/per_sample_norm_cov_plots/",tissue_list[i],"/",tissue_samples[j],"_norm_cov_plot.pdf"), height = 8, width = 15)
  }
}

