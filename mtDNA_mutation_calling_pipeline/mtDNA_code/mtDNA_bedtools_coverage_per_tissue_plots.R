library(tidyverse)
library(readxl)

mt_cov <- read.table(file = "/Users/al28/Mounts/cancer/mtDNA/whole_genome_coverage_pileup_and_bedtools_annotated.csv", sep = ",", header = T)

## All tissues
median_tissue <- vector(length = length(unique(mt_cov$Tissue)))
names(median_tissue) <- unique(mt_cov$Tissue)

for(i in 1:length(median_tissue)){
  median_tissue[i] <- median(na.rm = T, mt_cov$bedtools_mtDNA_genomes[which(mt_cov$Tissue == names(median_tissue)[i])]) 
}
median_tissue <- sort(median_tissue, decreasing = F)

median_tissue_cold <- median_tissue[which(!(grepl("warm_autopsy", names(median_tissue))))]

mt_cov_cold <- mt_cov[which(!(grepl("warm_autopsy", mt_cov$Tissue))),]
mt_cov_cold <- mt_cov_cold[which(!(is.na(mt_cov_cold$bedtools_autosome_coverage))),]

for(i in 1:nrow(mt_cov_cold)){
  mt_cov_cold$tissue_med_rank[i] <- which(names(median_tissue_cold) == mt_cov_cold$Tissue[i])
}

mt_cov_cold <- mt_cov_cold[order(mt_cov_cold$tissue_med_rank, mt_cov_cold$bedtools_mtDNA_genomes),]
mt_cov_cold$index <- 1:nrow(mt_cov_cold)

tissue_data <- as.data.frame(array(data = NA, dim = c(length(median_tissue_cold), 4)))
colnames(tissue_data) <- c("tissue","median","index_start","index_end")

tissue_data$tissue <- names(median_tissue_cold)
tissue_data$median <- median_tissue_cold
for(i in 1:nrow(tissue_data)){
  tissue_data$index_start[i] <- min(mt_cov_cold$index[which(mt_cov_cold$Tissue == tissue_data$tissue[i])])
  tissue_data$index_end[i] <- max(mt_cov_cold$index[which(mt_cov_cold$Tissue == tissue_data$tissue[i])])
}

ggplot() +
  geom_vline(xintercept = c(0.5,(tissue_data$index_end + 0.5)), color = "grey") +
  geom_point(data = mt_cov_cold, aes(x = index, y = bedtools_mtDNA_genomes), alpha = 0.5, size = 0.9) +
  geom_segment(aes(x = tissue_data$index_start,y = tissue_data$median,xend = tissue_data$index_end, yend = tissue_data$median), color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(expand = c(0,0), breaks = c(1,2500,5000,7500,10000,12500,nrow(mt_cov_cold))) +
  geom_text(data = tissue_data, aes(label = paste0(tissue,"\n",round(median)), x = (index_start + index_end) / 2, y = rep(x = c(42500,45000,40000), times = 8)), size = 2) +
  labs(x = "Index", y = "mtDNA genomes") +
  theme(plot.margin = unit(c(0.5,0.7,0.5,0.5),"cm"))
ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/bedtools_plots/per_tissue_mtDNA_genomes_exc_warm_autopsy.pdf", width = 15, height = 5)

ggplot() +
  geom_vline(xintercept = c(0.5,(tissue_data$index_end + 0.5)), color = "grey") +
  geom_point(data = mt_cov_cold, aes(x = index, y = bedtools_mtDNA_genomes), alpha = 0.5, size = 0.9) +
  geom_segment(aes(x = tissue_data$index_start,y = tissue_data$median,xend = tissue_data$index_end, yend = tissue_data$median), color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text(data = tissue_data, aes(label = paste0(tissue,"\n",round(median)), x = (index_start + index_end) / 2, y = rep(x = c(19000,20500,17500), times = 8)), size = 2) +
  scale_y_continuous(limits = c(0,21000)) +
  scale_x_continuous(expand = c(0,0), breaks = c(1,2500,5000,7500,10000,nrow(mt_cov_cold))) +
  labs(x = "Index", y = "mtDNA genomes") +
  theme(plot.margin = unit(c(0.5,0.7,0.5,0.5),"cm"))
ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/bedtools_plots/per_tissue_mtDNA_genomes_exc_warm_autopsy_capped_21000.pdf", width = 15, height = 5)

## Bladder
bladder_wgs_cov <- mt_cov_cold[which(mt_cov_cold$Tissue == "bladder_wgs"),]

median_bladder_wgs <- vector(length = length(unique(substr(bladder_wgs_cov$Sample,1,7))))
names(median_bladder_wgs) <- unique(substr(bladder_wgs_cov$Sample,1,7))

for(i in 1:length(median_bladder_wgs)){
  median_bladder_wgs[i] <- median(na.rm = T, bladder_wgs_cov$bedtools_mtDNA_genomes[which(substr(bladder_wgs_cov$Sample,1,7) == names(median_bladder_wgs)[i])]) 
}

median_bladder_wgs <- sort(median_bladder_wgs, decreasing = F)

for(i in 1:nrow(bladder_wgs_cov)){
  bladder_wgs_cov$donor_med_rank[i] <- which(names(median_bladder_wgs) == substr(bladder_wgs_cov$Sample[i],1,7))
}

bladder_wgs_cov_over_0_samples <- bladder_wgs_cov[which(bladder_wgs_cov$donor_med_rank %in% which(table(bladder_wgs_cov$donor_med_rank) > 0)),]
bladder_wgs_cov_over_0_samples <- bladder_wgs_cov_over_0_samples[order(bladder_wgs_cov_over_0_samples$donor_med_rank,bladder_wgs_cov_over_0_samples$bedtools_mtDNA_genomes),]
bladder_wgs_cov_over_0_samples$index <- 1:nrow(bladder_wgs_cov_over_0_samples)

bladder_wgs_donor_data <- as.data.frame(array(data = NA, dim = c(length(unique(substr(bladder_wgs_cov_over_0_samples$Sample,1,7))), 4)))
colnames(bladder_wgs_donor_data) <- c("donor","median","index_start","index_end")
bladder_wgs_donor_data$donor <- unique(substr(bladder_wgs_cov_over_0_samples$Sample,1,7))

for(i in 1:nrow(bladder_wgs_donor_data)){
  bladder_wgs_donor_data$median[i] <- median_bladder_wgs[which(names(median_bladder_wgs) == bladder_wgs_donor_data$donor[i])]
  bladder_wgs_donor_data$index_start[i] <- min(bladder_wgs_cov_over_0_samples$index[which(substr(bladder_wgs_cov_over_0_samples$Sample,1,7) == bladder_wgs_donor_data$donor[i])])
  bladder_wgs_donor_data$index_end[i] <- max(bladder_wgs_cov_over_0_samples$index[which(substr(bladder_wgs_cov_over_0_samples$Sample,1,7) == bladder_wgs_donor_data$donor[i])])
}

bladder_metadata <- read_xlsx("/Users/al28/Mounts/cancer/bladder_manuscript/supplementary_material/SuppTable2_MicrobiopsyInformation.xlsx")

for(i in 1:nrow(bladder_wgs_cov_over_0_samples)){
  bladder_wgs_cov_over_0_samples$histological_feature[i] <- bladder_metadata$histological_feature[which(bladder_metadata$data_access_id == bladder_wgs_cov_over_0_samples$Sample[i])]
}

ggplot() +
  geom_vline(xintercept = c(0.5,(bladder_wgs_donor_data$index_end + 0.5)), color = "grey") +
  geom_segment(aes(x = bladder_wgs_donor_data$index_start,y = bladder_wgs_donor_data$median,xend = bladder_wgs_donor_data$index_end, yend = bladder_wgs_donor_data$median), color = "red") +
  geom_point(data = bladder_wgs_cov_over_0_samples, aes(x = index, y = bedtools_mtDNA_genomes, color = histological_feature), alpha = 0.5, size = 1.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text(data = bladder_wgs_donor_data, aes(label = paste0(donor,"\n",round(median)), x = (index_start + index_end) / 2, y = rep(c(4000,4200), times = 10)), size = 2) +
  labs(x = "Index", y = "mtDNA genomes", color = "Histological feature") +
  scale_x_continuous(expand = c(0,0), breaks = c(1,25,50,75,100,nrow(bladder_wgs_cov_over_0_samples)))
ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/bedtools_plots/per_donor_bladder_wgs_samples_mtDNA_genomes.pdf", width = 15, height = 5)

## Blood Emily
blood_emily_cov <- mt_cov_cold[which(mt_cov_cold$Tissue == "blood_emily"),]

median_blood_emily <- vector(length = length(unique(substr(blood_emily_cov$Sample,1,7))))
names(median_blood_emily) <- unique(substr(blood_emily_cov$Sample,1,7))

for(i in 1:length(median_blood_emily)){
  median_blood_emily[i] <- median(na.rm = T, blood_emily_cov$bedtools_mtDNA_genomes[which(substr(blood_emily_cov$Sample,1,7) == names(median_blood_emily)[i])]) 
}

median_blood_emily <- sort(median_blood_emily, decreasing = F)

for(i in 1:nrow(blood_emily_cov)){
  blood_emily_cov$donor_med_rank[i] <- which(names(median_blood_emily) == substr(blood_emily_cov$Sample[i],1,7))
}

blood_emily_cov_over_100_samples <- blood_emily_cov[which(blood_emily_cov$donor_med_rank %in% which(table(blood_emily_cov$donor_med_rank) > 100)),]
blood_emily_cov_over_100_samples <- blood_emily_cov_over_100_samples[order(blood_emily_cov_over_100_samples$donor_med_rank,blood_emily_cov_over_100_samples$bedtools_mtDNA_genomes),]
blood_emily_cov_over_100_samples$index <- 1:nrow(blood_emily_cov_over_100_samples)

blood_emily_donor_data <- as.data.frame(array(data = NA, dim = c(length(unique(substr(blood_emily_cov_over_100_samples$Sample,1,7))), 4)))
colnames(blood_emily_donor_data) <- c("donor","median","index_start","index_end")
blood_emily_donor_data$donor <- unique(substr(blood_emily_cov_over_100_samples$Sample,1,7))

for(i in 1:nrow(blood_emily_donor_data)){
  blood_emily_donor_data$median[i] <- median_blood_emily[which(names(median_blood_emily) == blood_emily_donor_data$donor[i])]
  blood_emily_donor_data$index_start[i] <- min(blood_emily_cov_over_100_samples$index[which(substr(blood_emily_cov_over_100_samples$Sample,1,7) == blood_emily_donor_data$donor[i])])
  blood_emily_donor_data$index_end[i] <- max(blood_emily_cov_over_100_samples$index[which(substr(blood_emily_cov_over_100_samples$Sample,1,7) == blood_emily_donor_data$donor[i])])
}

ggplot() +
  geom_vline(xintercept = c(0.5,(blood_emily_donor_data$index_end + 0.5)), color = "grey") +
  geom_point(data = blood_emily_cov_over_100_samples, aes(x = index, y = bedtools_mtDNA_genomes), alpha = 0.5, size = 0.9) +
  geom_segment(aes(x = blood_emily_donor_data$index_start,y = blood_emily_donor_data$median,xend = blood_emily_donor_data$index_end, yend = blood_emily_donor_data$median), color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text(data = blood_emily_donor_data, aes(label = paste0(donor,"\n",round(median)), x = (index_start + index_end) / 2, y = 7900), size = 2) +
  labs(x = "Index", y = "mtDNA genomes") +
  scale_x_continuous(expand = c(0,0), breaks = c(1,500,1000,1500,2500,nrow(blood_emily_cov_over_100_samples))) +
  theme(plot.margin = unit(c(0.5,0.7,0.5,0.5),"cm"))
ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/bedtools_plots/per_donor_blood_emily_above_100_samples_mtDNA_genomes.pdf", width = 15, height = 5)

blood_emily_cov_under_100_samples <- blood_emily_cov[which(blood_emily_cov$donor_med_rank %in% which(table(blood_emily_cov$donor_med_rank) <= 100)),]
blood_emily_cov_under_100_samples <- blood_emily_cov_under_100_samples[order(blood_emily_cov_under_100_samples$donor_med_rank,blood_emily_cov_under_100_samples$bedtools_mtDNA_genomes),]
blood_emily_cov_under_100_samples$index <- 1:nrow(blood_emily_cov_under_100_samples)

blood_emily_donor_data <- as.data.frame(array(data = NA, dim = c(length(unique(substr(blood_emily_cov_under_100_samples$Sample,1,7))), 4)))
colnames(blood_emily_donor_data) <- c("donor","median","index_start","index_end")
blood_emily_donor_data$donor <- unique(substr(blood_emily_cov_under_100_samples$Sample,1,7))

for(i in 1:nrow(blood_emily_donor_data)){
  blood_emily_donor_data$median[i] <- median_blood_emily[which(names(median_blood_emily) == blood_emily_donor_data$donor[i])]
  blood_emily_donor_data$index_start[i] <- min(blood_emily_cov_under_100_samples$index[which(substr(blood_emily_cov_under_100_samples$Sample,1,7) == blood_emily_donor_data$donor[i])])
  blood_emily_donor_data$index_end[i] <- max(blood_emily_cov_under_100_samples$index[which(substr(blood_emily_cov_under_100_samples$Sample,1,7) == blood_emily_donor_data$donor[i])])
}

blood_emily_aa_metadata <- read.table("/Users/al28/Mounts/cancer/mtDNA/project_metadata/blood_emily/AA_samples_clonal.csv", sep = ",", stringsAsFactors = F, header = T)

blood_emily_cov_under_100_samples$clonal_status <- NA
for(i in 1:nrow(blood_emily_cov_under_100_samples)){
  if(blood_emily_cov_under_100_samples$Study[i] == 2178){
    blood_emily_cov_under_100_samples$clonal_status[i] <- "Alk"
  } else if(blood_emily_cov_under_100_samples$Study[i] == 2066){
    if(blood_emily_cov_under_100_samples$Sample[i] %in% blood_emily_aa_metadata$current_name){
      blood_emily_cov_under_100_samples$clonal_status[i] <- "AA_clonal"
    } else{
      blood_emily_cov_under_100_samples$clonal_status[i] <- "AA_non_clonal"
    }
  }
}

ggplot() +
  geom_vline(xintercept = c(0.5,(blood_emily_donor_data$index_end + 0.5)), color = "grey") +
  geom_point(data = blood_emily_cov_under_100_samples, aes(x = index, y = bedtools_mtDNA_genomes, color = clonal_status), alpha = 0.5, size = 0.9) +
  #  geom_segment(aes(x = blood_emily_donor_data$index_start,y = blood_emily_donor_data$median,xend = blood_emily_donor_data$index_end, yend = blood_emily_donor_data$median), color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text(data = blood_emily_donor_data, aes(label = paste0(donor,"\n",round(median)), x = (index_start + index_end) / 2, y = c(rep(c(20000,19000), times = 7), 20000)), size = 2) +
  labs(x = "Index", y = "mtDNA genomes", color = "Clonal colony") +
  scale_x_continuous(expand = c(0,0), breaks = c(1,25,50,75,97))
ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/bedtools_plots/per_donor_blood_emily_below_100_samples_mtDNA_genomes.pdf", width = 15, height = 5)

## Liver
liver_cov <- mt_cov_cold[which(mt_cov_cold$Tissue == "liver"),]

median_liver <- vector(length = length(unique(substr(liver_cov$Sample,1,7))))
names(median_liver) <- unique(substr(liver_cov$Sample,1,7))

for(i in 1:length(median_liver)){
  median_liver[i] <- median(na.rm = T, liver_cov$bedtools_mtDNA_genomes[which(substr(liver_cov$Sample,1,7) == names(median_liver)[i])]) 
}

median_liver <- sort(median_liver, decreasing = F)

for(i in 1:nrow(liver_cov)){
  liver_cov$donor_med_rank[i] <- which(names(median_liver) == substr(liver_cov$Sample[i],1,7))
}

liver_cov_over_20_samples <- liver_cov[which(liver_cov$donor_med_rank %in% which(table(liver_cov$donor_med_rank) > 20)),]
liver_cov_over_20_samples <- liver_cov_over_20_samples[order(liver_cov_over_20_samples$donor_med_rank,liver_cov_over_20_samples$bedtools_mtDNA_genomes),]
liver_cov_over_20_samples$index <- 1:nrow(liver_cov_over_20_samples)

liver_donor_data <- as.data.frame(array(data = NA, dim = c(length(unique(substr(liver_cov_over_20_samples$Sample,1,7))), 4)))
colnames(liver_donor_data) <- c("donor","median","index_start","index_end")
liver_donor_data$donor <- unique(substr(liver_cov_over_20_samples$Sample,1,7))

for(i in 1:nrow(liver_donor_data)){
  liver_donor_data$median[i] <- median_liver[which(names(median_liver) == liver_donor_data$donor[i])]
  liver_donor_data$index_start[i] <- min(liver_cov_over_20_samples$index[which(substr(liver_cov_over_20_samples$Sample,1,7) == liver_donor_data$donor[i])])
  liver_donor_data$index_end[i] <- max(liver_cov_over_20_samples$index[which(substr(liver_cov_over_20_samples$Sample,1,7) == liver_donor_data$donor[i])])
}

liver_metadata <- read.table("/Users/al28/Mounts/cancer/mtDNA/project_metadata/liver/liver.patients.info.for.andrew.csv", sep = ",", stringsAsFactors = F, header = T)

liver_cov_over_20_samples$Aetiology <- NA
for(i in 1:nrow(liver_cov_over_20_samples)){
  liver_cov_over_20_samples$Aetiology[i] <- liver_metadata$Aetiology[which(substr(liver_metadata$PD.ID,1,7) == substr(liver_cov_over_20_samples$Sample[i],1,7))]
}

ggplot() +
  geom_vline(xintercept = c(0.5,(liver_donor_data$index_end + 0.5)), color = "grey") +
  geom_point(data = liver_cov_over_20_samples, aes(x = index, y = bedtools_mtDNA_genomes, color = Aetiology), alpha = 0.5, size = 0.9) +
  geom_segment(aes(x = liver_donor_data$index_start,y = liver_donor_data$median,xend = liver_donor_data$index_end, yend = liver_donor_data$median), color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text(data = liver_donor_data, aes(label = paste0(donor,"\n",round(median)), x = (index_start + index_end) / 2, y = c(rep(c(12000,14000,10000), times = 10),12000,14000)), size = 2) +
  labs(x = "Index", y = "mtDNA genomes", color = "Aetiology") +
  scale_x_continuous(expand = c(0,0), breaks = c(1,250,500,750,1000,1250,1320))
ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/bedtools_plots/per_donor_liver_above_20_samples_mtDNA_genomes.pdf", width = 15, height = 5)

## Blood Foetal
blood_foetal_cov <- mt_cov_cold[which(mt_cov_cold$Tissue == "blood_foetal"),]

median_blood_foetal <- vector(length = length(unique(substr(blood_foetal_cov$Sample,1,7))))
names(median_blood_foetal) <- unique(substr(blood_foetal_cov$Sample,1,7))

for(i in 1:length(median_blood_foetal)){
  median_blood_foetal[i] <- median(na.rm = T, blood_foetal_cov$bedtools_mtDNA_genomes[which(substr(blood_foetal_cov$Sample,1,7) == names(median_blood_foetal)[i])]) 
}

median_blood_foetal <- sort(median_blood_foetal, decreasing = F)

for(i in 1:nrow(blood_foetal_cov)){
  blood_foetal_cov$donor_med_rank[i] <- which(names(median_blood_foetal) == substr(blood_foetal_cov$Sample[i],1,7))
}

blood_foetal_cov_over_0_samples <- blood_foetal_cov[which(blood_foetal_cov$donor_med_rank %in% which(table(blood_foetal_cov$donor_med_rank) > 0)),]
blood_foetal_cov_over_0_samples <- blood_foetal_cov_over_0_samples[order(blood_foetal_cov_over_0_samples$donor_med_rank,blood_foetal_cov_over_0_samples$bedtools_mtDNA_genomes),]
blood_foetal_cov_over_0_samples$index <- 1:nrow(blood_foetal_cov_over_0_samples)

blood_foetal_donor_data <- as.data.frame(array(data = NA, dim = c(length(unique(substr(blood_foetal_cov_over_0_samples$Sample,1,7))), 4)))
colnames(blood_foetal_donor_data) <- c("donor","median","index_start","index_end")
blood_foetal_donor_data$donor <- unique(substr(blood_foetal_cov_over_0_samples$Sample,1,7))

for(i in 1:nrow(blood_foetal_donor_data)){
  blood_foetal_donor_data$median[i] <- median_blood_foetal[which(names(median_blood_foetal) == blood_foetal_donor_data$donor[i])]
  blood_foetal_donor_data$index_start[i] <- min(blood_foetal_cov_over_0_samples$index[which(substr(blood_foetal_cov_over_0_samples$Sample,1,7) == blood_foetal_donor_data$donor[i])])
  blood_foetal_donor_data$index_end[i] <- max(blood_foetal_cov_over_0_samples$index[which(substr(blood_foetal_cov_over_0_samples$Sample,1,7) == blood_foetal_donor_data$donor[i])])
}

ggplot() +
  geom_vline(xintercept = c(0.5,(blood_foetal_donor_data$index_end + 0.5)), color = "grey") +
  geom_point(data = blood_foetal_cov_over_0_samples, aes(x = index, y = bedtools_mtDNA_genomes), alpha = 0.5, size = 0.9) +
  geom_segment(aes(x = blood_foetal_donor_data$index_start,y = blood_foetal_donor_data$median,xend = blood_foetal_donor_data$index_end, yend = blood_foetal_donor_data$median), color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text(data = blood_foetal_donor_data, aes(label = paste0(donor,"\n",round(median)), x = (index_start + index_end) / 2, y = 16000), size = 2) +
  scale_x_continuous(expand = c(0,0), breaks = c(1,100,200,300,400,500,545)) +
  labs(x = "Index", y = "mtDNA genomes") +
  theme(plot.margin = unit(c(0.5,0.7,0.5,0.5),"cm"))
ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/bedtools_plots/per_donor_blood_foetal_mtDNA_genomes.pdf", width = 15, height = 5)

## Blood Transplant
blood_transplant_cov <- mt_cov_cold[which(mt_cov_cold$Tissue == "blood_transplant"),]

median_blood_transplant <- vector(length = length(unique(substr(blood_transplant_cov$Sample,1,7))))
names(median_blood_transplant) <- unique(substr(blood_transplant_cov$Sample,1,7))

for(i in 1:length(median_blood_transplant)){
  median_blood_transplant[i] <- median(na.rm = T, blood_transplant_cov$bedtools_mtDNA_genomes[which(substr(blood_transplant_cov$Sample,1,7) == names(median_blood_transplant)[i])]) 
}

median_blood_transplant <- sort(median_blood_transplant, decreasing = F)

for(i in 1:nrow(blood_transplant_cov)){
  blood_transplant_cov$donor_med_rank[i] <- which(names(median_blood_transplant) == substr(blood_transplant_cov$Sample[i],1,7))
}

blood_transplant_cov_over_0_samples <- blood_transplant_cov[which(blood_transplant_cov$donor_med_rank %in% which(table(blood_transplant_cov$donor_med_rank) > 0)),]
blood_transplant_cov_over_0_samples <- blood_transplant_cov_over_0_samples[order(blood_transplant_cov_over_0_samples$donor_med_rank,blood_transplant_cov_over_0_samples$bedtools_mtDNA_genomes),]
blood_transplant_cov_over_0_samples$index <- 1:nrow(blood_transplant_cov_over_0_samples)

blood_transplant_donor_data <- as.data.frame(array(data = NA, dim = c(length(unique(substr(blood_transplant_cov_over_0_samples$Sample,1,7))), 4)))
colnames(blood_transplant_donor_data) <- c("donor","median","index_start","index_end")
blood_transplant_donor_data$donor <- unique(substr(blood_transplant_cov_over_0_samples$Sample,1,7))

for(i in 1:nrow(blood_transplant_donor_data)){
  blood_transplant_donor_data$median[i] <- median_blood_transplant[which(names(median_blood_transplant) == blood_transplant_donor_data$donor[i])]
  blood_transplant_donor_data$index_start[i] <- min(blood_transplant_cov_over_0_samples$index[which(substr(blood_transplant_cov_over_0_samples$Sample,1,7) == blood_transplant_donor_data$donor[i])])
  blood_transplant_donor_data$index_end[i] <- max(blood_transplant_cov_over_0_samples$index[which(substr(blood_transplant_cov_over_0_samples$Sample,1,7) == blood_transplant_donor_data$donor[i])])
}

ggplot() +
  geom_vline(xintercept = c(0.5,(blood_transplant_donor_data$index_end + 0.5)), color = "grey") +
  geom_point(data = blood_transplant_cov_over_0_samples, aes(x = index, y = bedtools_mtDNA_genomes, color = as.character(Study)), alpha = 0.5, size = 0.9) +
  geom_segment(aes(x = blood_transplant_donor_data$index_start,y = blood_transplant_donor_data$median,xend = blood_transplant_donor_data$index_end, yend = blood_transplant_donor_data$median), color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text(data = blood_transplant_donor_data, aes(label = paste0(donor,"\n",round(median)), x = (index_start + index_end) / 2, y = 2000), size = 2) +
  scale_x_continuous(expand = c(0,0), breaks = c(1,250,500,750,1000,1250,1405)) +
  labs(x = "Index", y = "mtDNA genomes", color = "Project") +
  theme(plot.margin = unit(c(0.5,0.7,0.5,0.5),"cm"))
ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/bedtools_plots/per_donor_blood_transplant_mtDNA_genomes.pdf", width = 15, height = 5)

## MUTYH mutant
muty_mutant_cov <- mt_cov_cold[which(mt_cov_cold$Tissue == "muty_mutant"),]

median_muty_mutant <- vector(length = length(unique(substr(muty_mutant_cov$Sample,1,7))))
names(median_muty_mutant) <- unique(substr(muty_mutant_cov$Sample,1,7))

for(i in 1:length(median_muty_mutant)){
  median_muty_mutant[i] <- median(na.rm = T, muty_mutant_cov$bedtools_mtDNA_genomes[which(substr(muty_mutant_cov$Sample,1,7) == names(median_muty_mutant)[i])]) 
}

median_muty_mutant <- sort(median_muty_mutant, decreasing = F)

for(i in 1:nrow(muty_mutant_cov)){
  muty_mutant_cov$donor_med_rank[i] <- which(names(median_muty_mutant) == substr(muty_mutant_cov$Sample[i],1,7))
}

muty_mutant_cov_over_0_samples <- muty_mutant_cov[which(muty_mutant_cov$donor_med_rank %in% which(table(muty_mutant_cov$donor_med_rank) > 0)),]
muty_mutant_cov_over_0_samples <- muty_mutant_cov_over_0_samples[order(muty_mutant_cov_over_0_samples$donor_med_rank,muty_mutant_cov_over_0_samples$bedtools_mtDNA_genomes),]
muty_mutant_cov_over_0_samples$index <- 1:nrow(muty_mutant_cov_over_0_samples)

muty_mutant_donor_data <- as.data.frame(array(data = NA, dim = c(length(unique(substr(muty_mutant_cov_over_0_samples$Sample,1,7))), 4)))
colnames(muty_mutant_donor_data) <- c("donor","median","index_start","index_end")
muty_mutant_donor_data$donor <- unique(substr(muty_mutant_cov_over_0_samples$Sample,1,7))

for(i in 1:nrow(muty_mutant_donor_data)){
  muty_mutant_donor_data$median[i] <- median_muty_mutant[which(names(median_muty_mutant) == muty_mutant_donor_data$donor[i])]
  muty_mutant_donor_data$index_start[i] <- min(muty_mutant_cov_over_0_samples$index[which(substr(muty_mutant_cov_over_0_samples$Sample,1,7) == muty_mutant_donor_data$donor[i])])
  muty_mutant_donor_data$index_end[i] <- max(muty_mutant_cov_over_0_samples$index[which(substr(muty_mutant_cov_over_0_samples$Sample,1,7) == muty_mutant_donor_data$donor[i])])
}

muty_mutant_metadata <- read.table("/Users/al28/Mounts/cancer/mtDNA/project_metadata/muty_mutant/sampleinfo_20201028.txt", sep = "\t", stringsAsFactors = F, header = T)

muty_mutant_cov_over_0_samples$tissue_type <- NA

for(i in 1:nrow(muty_mutant_cov_over_0_samples)){
  muty_mutant_cov_over_0_samples$tissue_type[i] <- muty_mutant_metadata$tissue[which(muty_mutant_metadata$sample == muty_mutant_cov_over_0_samples$Sample[i])]
}

ggplot() +
  geom_vline(xintercept = c(0.5,(muty_mutant_donor_data$index_end + 0.5)), color = "grey") +
  geom_point(data = muty_mutant_cov_over_0_samples, aes(x = index, y = bedtools_mtDNA_genomes, color = tissue_type), alpha = 0.5, size = 0.9) +
  geom_segment(aes(x = muty_mutant_donor_data$index_start,y = muty_mutant_donor_data$median,xend = muty_mutant_donor_data$index_end, yend = muty_mutant_donor_data$median), color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text(data = muty_mutant_donor_data, aes(label = paste0(donor,"\n",round(median)), x = (index_start + index_end) / 2, y = 4500), size = 2) +
  scale_x_continuous(expand = c(0,0), breaks = c(1,250,500,750,990)) +
  labs(x = "Index", y = "mtDNA genomes") +
  theme(plot.margin = unit(c(0.5,0.7,0.5,0.5),"cm"))
ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/bedtools_plots/per_donor_muty_mutant_mtDNA_genomes.pdf", width = 15, height = 5)

## Testes_tumour
testes_tumour_cov <- mt_cov_cold[which(mt_cov_cold$Tissue == "testes_tumour"),]

median_testes_tumour <- vector(length = length(unique(substr(testes_tumour_cov$Sample,1,7))))
names(median_testes_tumour) <- unique(substr(testes_tumour_cov$Sample,1,7))

for(i in 1:length(median_testes_tumour)){
  median_testes_tumour[i] <- median(na.rm = T, testes_tumour_cov$bedtools_mtDNA_genomes[which(substr(testes_tumour_cov$Sample,1,7) == names(median_testes_tumour)[i])]) 
}

median_testes_tumour <- sort(median_testes_tumour, decreasing = F)

for(i in 1:nrow(testes_tumour_cov)){
  testes_tumour_cov$donor_med_rank[i] <- which(names(median_testes_tumour) == substr(testes_tumour_cov$Sample[i],1,7))
}

testes_tumour_cov_over_1_samples <- testes_tumour_cov[which(testes_tumour_cov$donor_med_rank %in% which(table(testes_tumour_cov$donor_med_rank) > 1)),]
testes_tumour_cov_over_1_samples <- testes_tumour_cov_over_1_samples[order(testes_tumour_cov_over_1_samples$donor_med_rank,testes_tumour_cov_over_1_samples$bedtools_mtDNA_genomes),]
testes_tumour_cov_over_1_samples$index <- 1:nrow(testes_tumour_cov_over_1_samples)

testes_tumour_donor_data <- as.data.frame(array(data = NA, dim = c(length(unique(substr(testes_tumour_cov_over_1_samples$Sample,1,7))), 4)))
colnames(testes_tumour_donor_data) <- c("donor","median","index_start","index_end")
testes_tumour_donor_data$donor <- unique(substr(testes_tumour_cov_over_1_samples$Sample,1,7))

for(i in 1:nrow(testes_tumour_donor_data)){
  testes_tumour_donor_data$median[i] <- median_testes_tumour[which(names(median_testes_tumour) == testes_tumour_donor_data$donor[i])]
  testes_tumour_donor_data$index_start[i] <- min(testes_tumour_cov_over_1_samples$index[which(substr(testes_tumour_cov_over_1_samples$Sample,1,7) == testes_tumour_donor_data$donor[i])])
  testes_tumour_donor_data$index_end[i] <- max(testes_tumour_cov_over_1_samples$index[which(substr(testes_tumour_cov_over_1_samples$Sample,1,7) == testes_tumour_donor_data$donor[i])])
}

ggplot() +
  geom_vline(xintercept = c(0.5,(testes_tumour_donor_data$index_end + 0.5)), color = "grey") +
  geom_point(data = testes_tumour_cov_over_1_samples, aes(x = index, y = bedtools_mtDNA_genomes), alpha = 0.5, size = 0.9) +
  geom_segment(aes(x = testes_tumour_donor_data$index_start,y = testes_tumour_donor_data$median,xend = testes_tumour_donor_data$index_end, yend = testes_tumour_donor_data$median), color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text(data = testes_tumour_donor_data, aes(label = paste0(donor,"\n",round(median)), x = (index_start + index_end) / 2, y = 2000), size = 2) +
  labs(x = "Index", y = "mtDNA genomes") +
  scale_x_continuous(expand = c(0,0), breaks = c(1,25,50,75,100,125,150,nrow(testes_tumour_cov_over_1_samples))) +
  theme(plot.margin = unit(c(0.5,0.7,0.5,0.5),"cm"))
ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/bedtools_plots/per_donor_testes_tumour_above_1_samples_mtDNA_genomes.pdf", width = 15, height = 5)

## testes_wgs
testes_wgs_cov <- mt_cov_cold[which(mt_cov_cold$Tissue == "testes_wgs"),]

median_testes_wgs <- vector(length = length(unique(substr(testes_wgs_cov$Sample,1,7))))
names(median_testes_wgs) <- unique(substr(testes_wgs_cov$Sample,1,7))

for(i in 1:length(median_testes_wgs)){
  median_testes_wgs[i] <- median(na.rm = T, testes_wgs_cov$bedtools_mtDNA_genomes[which(substr(testes_wgs_cov$Sample,1,7) == names(median_testes_wgs)[i])]) 
}

median_testes_wgs <- sort(median_testes_wgs, decreasing = F)

for(i in 1:nrow(testes_wgs_cov)){
  testes_wgs_cov$donor_med_rank[i] <- which(names(median_testes_wgs) == substr(testes_wgs_cov$Sample[i],1,7))
}

testes_wgs_cov_over_1_samples <- testes_wgs_cov[which(testes_wgs_cov$donor_med_rank %in% which(table(testes_wgs_cov$donor_med_rank) > 1)),]
testes_wgs_cov_over_1_samples <- testes_wgs_cov_over_1_samples[order(testes_wgs_cov_over_1_samples$donor_med_rank,testes_wgs_cov_over_1_samples$bedtools_mtDNA_genomes),]
testes_wgs_cov_over_1_samples$index <- 1:nrow(testes_wgs_cov_over_1_samples)

testes_wgs_donor_data <- as.data.frame(array(data = NA, dim = c(length(unique(substr(testes_wgs_cov_over_1_samples$Sample,1,7))), 4)))
colnames(testes_wgs_donor_data) <- c("donor","median","index_start","index_end")
testes_wgs_donor_data$donor <- unique(substr(testes_wgs_cov_over_1_samples$Sample,1,7))

for(i in 1:nrow(testes_wgs_donor_data)){
  testes_wgs_donor_data$median[i] <- median_testes_wgs[which(names(median_testes_wgs) == testes_wgs_donor_data$donor[i])]
  testes_wgs_donor_data$index_start[i] <- min(testes_wgs_cov_over_1_samples$index[which(substr(testes_wgs_cov_over_1_samples$Sample,1,7) == testes_wgs_donor_data$donor[i])])
  testes_wgs_donor_data$index_end[i] <- max(testes_wgs_cov_over_1_samples$index[which(substr(testes_wgs_cov_over_1_samples$Sample,1,7) == testes_wgs_donor_data$donor[i])])
}

ggplot() +
  geom_vline(xintercept = c(0.5,(testes_wgs_donor_data$index_end + 0.5)), color = "grey") +
  geom_point(data = testes_wgs_cov_over_1_samples, aes(x = index, y = bedtools_mtDNA_genomes), alpha = 0.5, size = 0.9) +
  geom_segment(aes(x = testes_wgs_donor_data$index_start,y = testes_wgs_donor_data$median,xend = testes_wgs_donor_data$index_end, yend = testes_wgs_donor_data$median), color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text(data = testes_wgs_donor_data, aes(label = paste0(donor,"\n",round(median)), x = (index_start + index_end) / 2, y = c(rep(x = c(2600,2400), times = 10), 2600)), size = 2) +
  labs(x = "Index", y = "mtDNA genomes") +
  scale_x_continuous(expand = c(0,0), breaks = c(1,500,1000,1500,2500,nrow(testes_wgs_cov_over_1_samples))) +
  theme(plot.margin = unit(c(0.5,0.7,0.5,0.5),"cm"))
ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/bedtools_plots/per_donor_testes_wgs_above_1_samples_mtDNA_genomes.pdf", width = 15, height = 5)

