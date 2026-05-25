setwd("/Users/al28/Mounts/cancer/mtDNA/Analysis/shearwater/")

tissues <- list.files(".")

for(i in 1:length(tissues)){
  all_files <- list.files(tissues[i])
  input_files <- all_files[grepl("shearwater_calls.txt",all_files)]
  for(j in 1:length(input_files)){
    temp_data <- read.table(paste0(tissues[i],"/",input_files[j]), sep = "\t", stringsAsFactors = F, header = T)
    temp_data$donor <- gsub(pattern = "_shearwater_calls.txt", replacement = "", x = input_files[j])
    if(j == 1){
      tissue_data <- temp_data
    } else{
      tissue_data <- rbind(tissue_data,temp_data)
    }
  }
  tissue_data$tissue <- tissues[i]
  if(i == 1){
    all_data <- tissue_data
  } else{
    all_data <- rbind(all_data,tissue_data)
  }
}

all_data <- all_data[which(!(is.na(all_data$sampleID))),]
all_data$ref[which(all_data$ref == "TRUE")] <- "T"
all_data$mut[which(all_data$mut == "TRUE")] <- "T"

all_data <- all_data[-(which(all_data$sampleID == "PD45793b_lo0096" & all_data$tissue == "blood_transplant")),]
all_data <- all_data[-(which(all_data$sampleID == "P10B2_24" & all_data$tissue == "colon_ibd")),]
all_data <- all_data[-(which(all_data$sampleID == "Contamination01" & all_data$tissue == "prostate")),]

mutations <- all_data[which(all_data$vaf >= 0.01),]

wgs_coverage <- read.table("/Users/al28/Mounts/cancer/mtDNA/whole_genome_coverage_pileup_and_bedtools_annotated.csv", sep = ",", stringsAsFactors = F, header = T)

for(i in 1:nrow(wgs_coverage)){
  temp_muts <- mutations[which(mutations$sampleID == wgs_coverage$Sample[i] & mutations$tissue == wgs_coverage$Tissue[i]),]
  wgs_coverage$mut_count[i] <- nrow(temp_muts)
  wgs_coverage$SBS_count[i] <- nrow(temp_muts[which(temp_muts$ref %in% c("A","C","G","T") & temp_muts$mut %in% c("A","C","G","T")),])
  wgs_coverage$MBS_count[i] <- nrow(temp_muts[which(nchar(temp_muts$ref) > 1 & !(temp_muts$mut %in% c("INS","-"))),])
  wgs_coverage$INS_count[i] <- nrow(temp_muts[which(temp_muts$mut == "INS" | temp_muts$ref == "N"),])
  wgs_coverage$SB_Del_count[i] <- nrow(temp_muts[which(temp_muts$mut == "-" & nchar(temp_muts$ref) == 1),])
  wgs_coverage$MB_Del_count[i] <- nrow(temp_muts[which(temp_muts$mut == "-" & nchar(temp_muts$ref) > 1),])
  wgs_coverage$sum_vaf[i] <- sum(temp_muts$vaf)
}

wgs_coverage_mut_in_all_data <- wgs_coverage[which(paste(wgs_coverage$Tissue,wgs_coverage$Sample,sep = "_") %in% unique(paste(all_data$tissue,all_data$sampleID,sep = "_"))),]

library(tidyverse)

ggplot(data = wgs_coverage_mut_in_all_data) +
  geom_point(aes(x = as.numeric(wgs_coverage_mut_in_all_data$SeqX), y = as.numeric(wgs_coverage_mut_in_all_data$mtDNAgenomes), color = substr(wgs_coverage_mut_in_all_data$Tissue,1,11)), size = 0.2, alpha = 0.5) +
  theme_bw() +
  scale_color_manual(labels = c("bladder","blood_henry","blood_emily","blood_foetal","blood_transplant",
                              "breast","colon","colon_ibd","endometrium","liver","lung","pancreas",
                              "placenta","polymerase_mutant","prostate","stomach","testes","warm_autopsy"),
                     values = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")) +
  labs(x = "WGS coverage", y = "mtDNA genomes", color = "Tissue") +
  guides(colour = guide_legend(override.aes = list(size=3)))


coverage_calculated <- wgs_coverage_mut_in_all_data[which(!(is.na(wgs_coverage_mut_in_all_data$mtDNAgenomes)) & wgs_coverage_mut_in_all_data$SeqX >= 5),]

tissues <- as.character(unique(substr(coverage_calculated$Tissue,1,11)))
tissue_data <- as.data.frame(array(data = NA, dim = c(length(tissues),6)))
colnames(tissue_data) <- c("tissue","median_genomes","sample_count","cum_sample_count","cum_sample_start","ranked_median")

for(i in 1:nrow(tissue_data)){
  tissue_data$tissue[i] <- tissues[i]
  tissue_data$median_genomes[i] <- median(coverage_calculated$mtDNAgenomes[which(substr(coverage_calculated$Tissue,1,11) == tissue_data$tissue[i])])
  tissue_data$sample_count[i] <- length(which(substr(coverage_calculated$Tissue,1,11) == tissue_data$tissue[i]))
  tissue_data$median_sum_vaf[i] <- median(coverage_calculated$sum_vaf[which(substr(coverage_calculated$Tissue,1,11) == tissue_data$tissue[i])])
  tissue_data$median_mut_count[i] <- median(coverage_calculated$mut_count[which(substr(coverage_calculated$Tissue,1,11) == tissue_data$tissue[i])])
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
  coverage_calculated$tissue_ranked[i] <- tissue_data$ranked_median[which(tissue_data$tissue == substr(coverage_calculated$Tissue[i],1,11))] 
}

coverage_calculated <- coverage_calculated[order(coverage_calculated$tissue_ranked,coverage_calculated$mtDNAgenomes,decreasing = F),]
coverage_calculated$index <- 1:nrow(coverage_calculated)

pcawg_coverage <- read.table(file = "/Users/al28/Mounts/cancer/mtDNA/pcawg_mtdna_copy_number.csv", sep = ",", header = T)

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

for(i in 1:nrow(tissue_data)){
  if(i %in% 1:(nrow(tissue_data) - 1)){
    tissue_data$lab_position[i] <- (tissue_data$cum_sample_start[i] + tissue_data$cum_sample_start[i + 1]) / 2
  } else{
    tissue_data$lab_position[i] <- ((tissue_data$cum_sample_start[i] + nrow(coverage_calculated)) / 2)
  }
}

ggplot() +
  geom_point(data = pcawg_coverage, aes(x = index, y = mt_copy_number), color = "red", alpha = 0.5, size = 0.9) +
  geom_point(data = coverage_calculated, aes(x = index, y = mtDNAgenomes), color = "blue", alpha = 0.5, size = 0.9) +
  geom_vline(xintercept = c(0.5,(tissue_data$cum_sample_count + 0.5))) +
  geom_label(data = tissue_data, aes(x = lab_position, y = rep(c(8500,9000,8000), times = 6), label = tissue), size = 3) +
  geom_segment(aes(x = tissue_data$cum_sample_start,y = tissue_data$median_genomes,xend = tissue_data$cum_sample_count, yend = tissue_data$median_genomes), color = "blue") +
  geom_segment(aes(x = tissue_data$cum_sample_start,y = tissue_data$median_pcawg,xend=tissue_data$cum_sample_count,yend=tissue_data$median_pcawg), color = "red") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Sample", y = "mtDNA genomes")

tissue_data$ranked_median_sum_vaf[order(tissue_data$median_sum_vaf)] <- 1:length(tissues)
tissue_data <- tissue_data[order(tissue_data$ranked_median_sum_vaf, decreasing = F),]

for(i in 1:nrow(tissue_data)){
  tissue_data$cum_sample_count_sum_vaf[i] <- sum(tissue_data$sample_count[1:i])
  if(i == 1){
    tissue_data$cum_sample_start_sum_vaf[i] <- 1
  } else{
    tissue_data$cum_sample_start_sum_vaf[i] <- tissue_data$cum_sample_count_sum_vaf[i - 1] + 1
  }
}

coverage_calculated$tissue_ranked_sum_vaf <- NA

for(i in 1:nrow(coverage_calculated)){
  coverage_calculated$tissue_ranked_sum_vaf[i] <- tissue_data$ranked_median_sum_vaf[which(tissue_data$tissue == substr(coverage_calculated$Tissue[i],1,11))] 
}

coverage_calculated <- coverage_calculated[order(coverage_calculated$tissue_ranked_sum_vaf,coverage_calculated$sum_vaf,decreasing = F),]
coverage_calculated$index_sum_vaf <- 1:nrow(coverage_calculated)

for(i in 1:nrow(tissue_data)){
  if(i %in% 1:(nrow(tissue_data) - 1)){
    tissue_data$lab_position[i] <- (tissue_data$cum_sample_start_sum_vaf[i] + tissue_data$cum_sample_start_sum_vaf[i + 1]) / 2
  } else{
    tissue_data$lab_position[i] <- ((tissue_data$cum_sample_start_sum_vaf[i] + nrow(coverage_calculated)) / 2)
  }
}

ggplot() +
  geom_point(data = coverage_calculated, aes(x = index_sum_vaf, y = sum_vaf), color = "blue", alpha = 0.5, size = 0.9) +
  geom_vline(xintercept = c(0.5,(tissue_data$cum_sample_count_sum_vaf + 0.5))) +
  geom_label(data = tissue_data, aes(x = lab_position, y = rep(c(7,6.8,7.2), times = 6), label = tissue), size = 3) +
  geom_segment(aes(x = tissue_data$cum_sample_start_sum_vaf,y = tissue_data$median_sum_vaf,xend = tissue_data$cum_sample_count_sum_vaf, yend = tissue_data$median_sum_vaf), color = "blue") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Sample", y = "Sum of VAFs")

ggplot() +
  geom_point(data = coverage_calculated, aes(x = index_sum_vaf, y = sum_vaf), color = "blue", alpha = 0.5, size = 0.9) +
  geom_vline(xintercept = c(0.5,(tissue_data$cum_sample_count_sum_vaf + 0.5))) +
  geom_label(data = tissue_data, aes(x = lab_position, y = rep(c(7,6.8,7.2), times = 6), label = tissue), size = 3) +
  geom_segment(aes(x = tissue_data$cum_sample_start_sum_vaf,y = tissue_data$median_sum_vaf,xend = tissue_data$cum_sample_count_sum_vaf, yend = tissue_data$median_sum_vaf), color = "blue") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Sample", y = "Sum of VAFs")


emily_intercept <- c(154,
                     154 + 242,
                     154 + 242 + 366,
                     154 + 242 + 366 + 356,
                     154 + 242 + 366 + 356 + 179,
                     154 + 242 + 366 + 356 + 179 + 174,
                     154 + 242 + 366 + 356 + 179 + 174 + 94,
                     154 + 242 + 366 + 356 + 179 + 174 + 94 + 194,
                     154 + 242 + 366 + 356 + 179 + 174 + 94 + 194 + 302)

  ggplot() +
    geom_point(data = emily_blood, aes(x = index, y = sum_vaf), color = "blue", alpha = 0.5, size = 0.9) +
    geom_vline(xintercept = c(0.5,emily_intercept + 0.5)) +
    geom_segment(aes(x = 1, xend = 154, y = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD40315")]), yend = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD40315")])), color = "black") +
    geom_segment(aes(x = 154 + 1, xend = 154 + 242, y = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD45517")]), yend = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD45517")])), color = "black") +
    geom_segment(aes(x = 154 + 242 + 1, xend = 154 + 242 + 366, y = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD40521")]), yend = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD40521")])), color = "black") +
    geom_segment(aes(x = 154 + 242 + 366 + 1, xend = 154 + 242 + 366 + 356, y = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD40667")]), yend = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD40667")])), color = "black") +
    geom_segment(aes(x = 154 + 242 + 366 + 356 + 1, xend = 154 + 242 + 366 + 356 + 179, y = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD41048")]), yend = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD41048")])), color = "black") +
    geom_segment(aes(x = 154 + 242 + 366 + 356 + 179 + 1, xend = 154 + 242 + 366 + 356 + 179 + 174, y = median(emily_blood$sum_vaf[which(emily_blood$donor == "BMH1_TG")]), yend = median(emily_blood$sum_vaf[which(emily_blood$donor == "BMH1_TG")])), color = "black") +
    geom_segment(aes(x = 154 + 242 + 366 + 356 + 179 + 1, xend = 154 + 242 + 366 + 356 + 179 + 174, y = median(coverage_calculated$sum_vaf[which(coverage_calculated$Tissue == "blood")]), yend = median(coverage_calculated$sum_vaf[which(coverage_calculated$Tissue == "blood")])), color = "red") +
    geom_segment(aes(x = 154 + 242 + 366 + 356 + 179 + 174 + 1, xend = 154 + 242 + 366 + 356 + 179 + 174 + 94,  y = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD44579")]), yend = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD44579")])), color = "black") +
    geom_segment(aes(x = 154 + 242 + 366 + 356 + 179 + 174 + 94 + 1, xend = 154 + 242 + 366 + 356 + 179 + 174 + 94 + 194,  y = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD45534")]), yend = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD45534")])), color = "black") +
    geom_segment(aes(x = 154 + 242 + 366 + 356 + 179 + 174 + 94 + 194 + 1, xend = 154 + 242 + 366 + 356 + 179 + 174 + 94 + 194 + 302,  y = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD43974")]), yend = median(emily_blood$sum_vaf[which(emily_blood$donor == "PD43974")])), color = "black") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Sample", y = "Sum of VAFs")
  
liver_muts <- mutations[which(mutations$tissue == "liver"),]

liver_muts$bad <- F
liver_muts$bad[which(liver_muts$pos %in% which(prop_sample_subs$liver >= 0.25))] <- T


## MUTYH mutant
muty_mutant_cov <- wgs_coverage[which(wgs_coverage$Tissue == "muty_mutant"),]

median_muty_mutant <- vector(length = length(unique(substr(muty_mutant_cov$Sample,1,7))))
names(median_muty_mutant) <- unique(substr(muty_mutant_cov$Sample,1,7))

for(i in 1:length(median_muty_mutant)){
  median_muty_mutant[i] <- median(na.rm = T, muty_mutant_cov$sum_vaf[which(substr(muty_mutant_cov$Sample,1,7) == names(median_muty_mutant)[i])]) 
}

median_muty_mutant <- sort(median_muty_mutant, decreasing = F)

for(i in 1:nrow(muty_mutant_cov)){
  muty_mutant_cov$donor_med_rank[i] <- which(names(median_muty_mutant) == substr(muty_mutant_cov$Sample[i],1,7))
}

muty_mutant_cov_over_10_samples <- muty_mutant_cov[which(muty_mutant_cov$donor_med_rank %in% which(table(muty_mutant_cov$donor_med_rank) > 0)),]
muty_mutant_cov_over_10_samples <- muty_mutant_cov_over_10_samples[order(muty_mutant_cov_over_10_samples$donor_med_rank,muty_mutant_cov_over_10_samples$sum_vaf),]
muty_mutant_cov_over_10_samples$index <- 1:nrow(muty_mutant_cov_over_10_samples)

muty_mutant_donor_data <- as.data.frame(array(data = NA, dim = c(length(unique(substr(muty_mutant_cov_over_10_samples$Sample,1,7))), 4)))
colnames(muty_mutant_donor_data) <- c("donor","median","index_start","index_end")
muty_mutant_donor_data$donor <- unique(substr(muty_mutant_cov_over_10_samples$Sample,1,7))

for(i in 1:nrow(muty_mutant_donor_data)){
  muty_mutant_donor_data$median[i] <- median_muty_mutant[which(names(median_muty_mutant) == muty_mutant_donor_data$donor[i])]
  muty_mutant_donor_data$index_start[i] <- min(muty_mutant_cov_over_10_samples$index[which(substr(muty_mutant_cov_over_10_samples$Sample,1,7) == muty_mutant_donor_data$donor[i])])
  muty_mutant_donor_data$index_end[i] <- max(muty_mutant_cov_over_10_samples$index[which(substr(muty_mutant_cov_over_10_samples$Sample,1,7) == muty_mutant_donor_data$donor[i])])
}

muty_mutant_metadata <- read.table("/Users/al28/Mounts/cancer/mtDNA/project_metadata/muty_mutant/sampleinfo_20201028.txt", sep = "\t", stringsAsFactors = F, header = T)

muty_mutant_cov_over_10_samples$tissue_type <- NA

for(i in 1:nrow(muty_mutant_cov_over_10_samples)){
  muty_mutant_cov_over_10_samples$tissue_type[i] <- muty_mutant_metadata$tissue[which(muty_mutant_metadata$sample == muty_mutant_cov_over_10_samples$Sample[i])]
}

ggplot() +
  geom_vline(xintercept = c(0.5,(muty_mutant_donor_data$index_end + 0.5)), color = "grey") +
  geom_point(data = muty_mutant_cov_over_10_samples, aes(x = index, y = sum_vaf, color = tissue_type), alpha = 0.5, size = 0.9) +
  geom_segment(aes(x = muty_mutant_donor_data$index_start,y = muty_mutant_donor_data$median,xend = muty_mutant_donor_data$index_end, yend = muty_mutant_donor_data$median), color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text(data = muty_mutant_donor_data, aes(label = paste0(donor,"\n",round(median, digits = 3)), x = (index_start + index_end) / 2, y = 3.5), size = 2) +
  scale_x_continuous(expand = c(0,0), breaks = c(1,250,500,750,990)) +
  labs(x = "Index", y = "Sum VAF") +
  theme(plot.margin = unit(c(0.5,0.7,0.5,0.5),"cm"))


