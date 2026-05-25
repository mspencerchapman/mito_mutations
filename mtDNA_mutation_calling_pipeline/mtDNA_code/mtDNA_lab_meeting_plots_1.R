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

setwd("../shearwater_tissue_overview/")

mutations <- all_data[which(all_data$vaf >= 0.01),]
subs <- mutations[which(mutations$ref %in% c("A","C","G","T") & mutations$mut %in% c("A","C","G","T")),]
indels <- mutations[which(mutations$mut %in% c("-","INS")),]

prop_sample_subs <- as.data.frame(array(dim = c(16569,length(tissues[which(substr(tissues,1,4) != "warm")]) + 1),data = NA))
colnames(prop_sample_subs) <- c("pos",tissues[which(substr(tissues,1,4) != "warm")])
prop_sample_subs$pos <- 1:16569

for(j in 1:length(tissues[which(substr(tissues,1,4) != "warm")])){
  sample_count <- length(unique(all_data$sampleID[which(all_data$tissue == tissues[j])]))
  for(i in 1:16569){
    prop_sample_subs[i,tissues[j]] <- length(unique(subs$sampleID[which(subs$tissue == tissues[j] & subs$pos == i)]))/ sample_count
  }
}

indel_array <- as.data.frame(array(data = 0, dim = c(16569,length(unique(all_data$sampleID)))))
colnames(indel_array) <- unique(all_data$sampleID)
for(i in 1:nrow(indels)){
  indel_start <- indels$pos[i]
  indel_length <- nchar(indels$ref[i])
  indel_end <- indels$pos[i] + indel_length - 1
  indel_array[indel_start:indel_end,indels$sampleID[i]] <- 1
}

prop_sample_indels <- as.data.frame(array(dim = c(16569,length(tissues[which(substr(tissues,1,4) != "warm")]) + 1),data = NA))
colnames(prop_sample_indels) <- c("pos",tissues[which(substr(tissues,1,4) != "warm")])
prop_sample_indels$pos <- 1:16569

for(j in 1:length(tissues[which(substr(tissues,1,4) != "warm")])){
  sample_count <- length(unique(all_data$sampleID[which(all_data$tissue == tissues[j])]))
  indel_temp <- indel_array[1:16569,c(unique(all_data$sampleID[which(all_data$tissue == tissues[j])]))]
  prop_sample_indels[,tissues[j]] <- rowSums(indel_temp) / sample_count
}

mtDNA_features <- read.table("/Users/al28/Mounts/cancer/mtDNA/mtDNA_genome_feature_annotation.csv", sep = ",", stringsAsFactors = F, header = T)
mtDNA_features$colour <- NA
mtDNA_features$colour[which(mtDNA_features$feature_type == "tRNA")] <- "grey50"
mtDNA_features$colour[which(mtDNA_features$feature_type == "ND_coding_gene")] <- "firebrick"
mtDNA_features$colour[which(mtDNA_features$feature_type == "rRNA")] <- "darkorchid"
mtDNA_features$colour[which(mtDNA_features$feature_type == "ATP_coding_gene")] <- "steelblue"
mtDNA_features$colour[which(mtDNA_features$feature_type == "COX_coding_gene")] <- "springgreen3"
mtDNA_features$colour[which(mtDNA_features$feature_type == "CYTB_coding_gene")] <- "sienna1"
mtDNA_features$colour[which(mtDNA_features$feature_type == "control_element")] <- "thistle1"
mtDNA_features$colour[which(mtDNA_features$feature_type == "conserved_seq_block")] <- "slateblue"
mtDNA_features$colour[which(mtDNA_features$feature_type == "replication_origin")] <- "turquoise"
mtDNA_features$colour[which(mtDNA_features$feature_type == "promoter")] <- "skyblue3"

library("circlize")

for(j in 1:length(tissues[which(substr(tissues,1,4) != "warm")])){
  circos.par(start.degree = 90, gap.degree = 0, gap.after = 0, cell.padding = c(0,0), track.height = 0.1)
  circos.initialize(factors = c("a"), xlim = c(0,16569))
  circos.track(ylim = c(-0.2,1), bg.border = "white")
  circos.rect(xleft = 0, xright = 16569, ybottom = 0.9, ytop = 1, col = "black")
  circos.rect(xleft = mtDNA_features$start_pos, xright = mtDNA_features$end_pos, ybottom = 0, ytop = 0.9, col = mtDNA_features$colour)
  circos.axis(major.at = seq(from = 0, to = 16589, by = 1000))
  circos.text(x = (mtDNA_features$start_pos[which(grepl("coding_gene",mtDNA_features$feature_type))] + mtDNA_features$end_pos[which(grepl("coding_gene",mtDNA_features$feature_type))]) / 2, y = 0.45, labels = gsub(pattern = "MT-",replacement = "",x = mtDNA_features$map_locus[which(grepl("coding_gene",mtDNA_features$feature_type))]), cex = 0.5)
  circos.text(x = (mtDNA_features$start_pos[which(grepl("tRNA",mtDNA_features$feature_type))] + mtDNA_features$end_pos[which(grepl("tRNA",mtDNA_features$feature_type))]) / 2, y = -0.2, labels = gsub(pattern = "MT-T",replacement = "",x = mtDNA_features$map_locus[which(grepl("tRNA",mtDNA_features$feature_type))]), cex = 0.25)
  circos.track(ylim = c(0,1), bg.border = "white")
  for(i in 1:16569){
    circos.rect(xleft = prop_sample_subs$pos[i], xright = prop_sample_subs$pos[i], ybottom = 0, ytop = prop_sample_subs[i,tissues[which(substr(tissues,1,4) != "warm")][j]], col = "blue", border = "blue")
  }
  circos.track(ylim = c(0,1), bg.border = "white")
  for(i in 1:16569){
    circos.rect(xleft = prop_sample_indels$pos[i], xright = prop_sample_indels$pos[i], ybottom = 0, ytop = prop_sample_indels[i,tissues[which(substr(tissues,1,4) != "warm")][j]], col = "red", border = "red")
  }
  dev.copy2pdf(file = paste0(tissues[which(substr(tissues,1,4) != "warm")][j],"/circular_indel_sample_prop_",tissues[which(substr(tissues,1,4) != "warm")][j],"_plot.pdf"),out.type = "pdf")
  circos.clear()
}

