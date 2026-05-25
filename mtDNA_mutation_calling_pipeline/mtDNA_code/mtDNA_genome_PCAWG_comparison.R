library("tidyverse")

setwd("/Users/al28/Mounts/cancer/mtDNA/")

pcawg_coverage <- read.table(file = "pcawg_mtdna_copy_number.csv", sep = ",", header = T)
wgs_coverage <- read.table(file = "whole_genome_coverage.csv", sep = ",", header = T)
wgs_coverage$mtDNAcoverage <- NA

for(i in 1:nrow(wgs_coverage)){
  if(file.exists(paste0("Data/",wgs_coverage$Tissue[i],"/",wgs_coverage$Sample[i],"_MT_count.csv"))){
    sample_cov <- read.table(paste0("Data/",wgs_coverage$Tissue[i],"/",wgs_coverage$Sample[i],"_MT_count.csv"), sep = ",", header = T)
    wgs_coverage$mtDNAcoverage[i] <- sum(sample_cov) / 16569
   }
}

wgs_coverage$mtDNAgenomes <- 2 * (wgs_coverage$mtDNAcoverage / wgs_coverage$SeqX)

coverage_calculated <- wgs_coverage[which(!(is.na(wgs_coverage$mtDNAgenomes)) & wgs_coverage$SeqX >= 5),]

tissues <- as.character(unique(coverage_calculated$Tissue))

tissue_data <- as.data.frame(array(data = NA, dim = c(length(tissues),6)))
colnames(tissue_data) <- c("tissue","median_genomes","sample_count","cum_sample_count","cum_sample_start","ranked_median")

for(i in 1:nrow(tissue_data)){
 tissue_data$tissue[i] <- tissues[i]
 tissue_data$median_genomes[i] <- median(coverage_calculated$mtDNAgenomes[which(coverage_calculated$Tissue == tissue_data$tissue[i])])
 tissue_data$sample_count[i] <- length(which(coverage_calculated$Tissue == tissue_data$tissue[i]))
}

tissue_data$ranked_median[order(tissue_data$median_genomes)] <- 1:length(tissues)
tissue_data <- tissue_data[order(tissue_data$ranked_median, decreasing = F),]

for(i in 1:nrow(tissue_data)){
  tissue_data$cum_sample_count[i] <- sum(tissue_data$sample_count[1:i])
  if(i == 1){
    tissue_data$cum_sample_start[i] <- 1
  } else{
    tissue_data$cum_sample_start[i] <- tissue_data$cum_sample_count[i - 1] + 1
  }
}

coverage_calculated$tissue_ranked <- NA

for(i in 1:nrow(coverage_calculated)){
  coverage_calculated$tissue_ranked[i] <- tissue_data$ranked_median[which(tissue_data$tissue == coverage_calculated$Tissue[i])] 
}

coverage_calculated <- coverage_calculated[order(coverage_calculated$tissue_ranked,coverage_calculated$mtDNAgenomes,decreasing = F),]
coverage_calculated$index <- 1:nrow(coverage_calculated)

pcawg_coverage$tissue_ranked <- NA

for(i in 1:nrow(pcawg_coverage)){
  pcawg_coverage$tissue_ranked[i] <- tissue_data$ranked_median[which(tissue_data$tissue == pcawg_coverage$combined_category[i])]
}

pcawg_coverage <- pcawg_coverage[order(pcawg_coverage$tissue_ranked, pcawg_coverage$mt_copy_number, decreasing = F),]

tissue_data$pcawg_index_increment <- NA

for(i in 1:nrow(tissue_data)){
  tissue_data$median_pcawg[i] <- median(pcawg_coverage$mt_copy_number[which(pcawg_coverage$combined_category == tissue_data$tissue[i])])
  tissue_data$pcawg_count[i] <- length(which(pcawg_coverage$combined_category == tissue_data$tissue[i]))
  if(tissue_data$pcawg_count[i] != 0){
    tissue_data$pcawg_index_increment[i] <- (tissue_data$sample_count[i] - 1) / (tissue_data$pcawg_count[i] - 1)
  }
}

pcawg_coverage$index <- NA

for(i in 1:nrow(pcawg_coverage)){
  if(i == 1){
    pcawg_coverage$index[i] <- tissue_data$cum_sample_start[which(tissue_data$tissue == pcawg_coverage$combined_category[i])]
  } else if(pcawg_coverage$combined_category[i] == pcawg_coverage$combined_category[i - 1]){
    pcawg_coverage$index[i] <- pcawg_coverage$index[i - 1] + tissue_data$pcawg_index_increment[which(tissue_data$tissue == pcawg_coverage$combined_category[i])]
  } else{
    pcawg_coverage$index[i] <- tissue_data$cum_sample_start[which(tissue_data$tissue == pcawg_coverage$combined_category[i])]
  }
}

ggplot() +
  geom_point(data = pcawg_coverage, aes(x = index, y = mt_copy_number), color = "red", alpha = 0.5, size = 0.9) +
  geom_point(data = coverage_calculated, aes(x = index, y = mtDNAgenomes), color = "blue", alpha = 0.5, size = 0.9) +
  geom_vline(xintercept = c(0.5,(tissue_data$cum_sample_count + 0.5))) +
  geom_segment(aes(x = tissue_data$cum_sample_start,y = tissue_data$median_genomes,xend = tissue_data$cum_sample_count, yend = tissue_data$median_genomes), color = "blue") +
  geom_segment(aes(x = tissue_data$cum_sample_start,y = tissue_data$median_pcawg,xend=tissue_data$cum_sample_count,yend=tissue_data$median_pcawg), color = "red") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(panel.grid = element_blank())
