# library("dndscv")
# buildref(cdsfile="mtDNA_all_protein_coding_genes.txt", genomefile="/lustre/scratch125/casm/team268im/al28/bed_ref/hs37d5.fa", outfile="mtref.rda", numcode=2)
# dnds_germline <- dndscv(mutations= unique_germline, numcode=2, refdb="mtref.rda", max_muts_per_gene_per_sample=Inf, cv=NULL)
# dnds_somatic <- dndscv(mutations= unique_somatic, numcode=2, refdb="mtref.rda", max_muts_per_gene_per_sample=Inf, cv=NULL)

library("tidyverse")

setwd("/Users/al28/Mounts/cancer/mtDNA/Analysis/global_vaf_heatmaps/")

tissue_list <- list.files(".")

for(z in 1:length(tissue_list)){
  if(file.exists(paste0(tissue_list[z],"/",tissue_list[z],"_combined_somatic_indels_merged.tsv"))){
    if(!(exists("all_tissue_somatic"))){
      all_tissue_somatic <- read.table(file = paste0(tissue_list[z],"/",tissue_list[z],"_combined_somatic_indels_merged.tsv"), header = T, stringsAsFactors = F, sep = "\t")
    } else{
      all_tissue_somatic <- rbind(all_tissue_somatic, read.table(file = paste0(tissue_list[z],"/",tissue_list[z],"_combined_somatic_indels_merged.tsv"), header = T, stringsAsFactors = F, sep = "\t"))
    }
  }
  
  if(file.exists(paste0(tissue_list[z],"/",tissue_list[z],"_combined_germline_indels_merged.tsv"))){
    if(!(exists("all_tissue_germline"))){
      all_tissue_germline <- read.table(file = paste0(tissue_list[z],"/",tissue_list[z],"_combined_germline_indels_merged.tsv"), header = T, stringsAsFactors = F, sep = "\t")
    } else{
      all_tissue_germline <- rbind(all_tissue_germline, read.table(file = paste0(tissue_list[z],"/",tissue_list[z],"_combined_germline_indels_merged.tsv"), header = T, stringsAsFactors = F, sep = "\t"))
    }
  }
}

unique_germline <- all_tissue_germline[which(!(duplicated(all_tissue_germline$mut_id))),]
unique_somatic <- all_tissue_somatic[which(!(duplicated(all_tissue_somatic$mut_id))),]

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


ggplot() +
  geom_rect(mapping = aes(xmin = 1, xmax = 16569, ymin = 10, ymax = 10.1), fill = "black") +
  geom_rect(data = mtDNA_features, mapping = aes(xmin = start_pos, xmax = end_pos, ymin = 9, ymax = 10, fill = feature_type), color = "black") +
  geom_text() +
  scale_y_continuous(limits = c(0,10.5)) +
  scale_x_continuous(limits = c(0,16569), breaks = seq(from = 0, to = 16000, by = 1000), labels = seq(from = 0, to = 16000, by = 1000)) +
  coord_polar() +
  theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
  


ggplot(data = unique_somatic[which(unique_somatic$tissue == "liver"),]) +
  geom_histogram(aes(x = pos), binwidth = 50) +
  theme_bw()

liver_donor_list <- c("PD36713","PD36714","PD36715","PD36717","PD36718","PD37105","PD37107","PD37110","PD37111","PD37113","PD37114","PD37115","PD37116","PD37118","PD37230","PD37231","PD37232","PD37234","PD37237","PD37238","PD37239","PD37240","PD37243","PD37245","PD37904","PD37907","PD37909","PD37910","PD37914","PD37915","PD37917","PD37918")
blood_emily_donor_list <- c("BMH1_TG","PD40521","PD40667","PD41048","PD43974","PD44579","PD45534")
blood_transplant_donor_list <- c("PD45792","PD45793","PD45794","PD45795","PD45798","PD45804","PD45805")
colon_donor_list <- c("HLS","O174","O176","O187","O189","O193","O208","O225","O232","O244","O251","O325","O333","O337","O340","O342","O395","O401","O438","OO10","OO13","OO82","OO94","OO99","OOO8","PD28690","PD34199","PD34200","PD34201","PD37266","PD37449","PD37513","PD37514","PD37590")
colon_ibd_donor_list <- c("P1","P3","P10","P12","P15","P16","P17","P22","P25","P28","P29","P30","P31","P33","P34","P35","P36","P37","P38","P39","P40","P41","P43","P44","P45","P46","P48","P50","P51","P52","P53","P54")

prop_donor <- as.data.frame(array(dim = c(16569,4),data = NA))
colnames(prop_donor) <- c("pos","prop_liver","prop_blood","prop_colon")
prop_donor$pos <- 1:16569

for(i in 1:16569){
  prop_donor[i,"prop_liver"] <- length(unique(unique_somatic$donor[which(unique_somatic$pos == i & unique_somatic$tissue == "liver" & unique_somatic$donor %in% liver_donor_list)])) / length(liver_donor_list)
  prop_donor[i,"prop_blood"] <- length(unique(unique_somatic$donor[which(unique_somatic$pos == i & unique_somatic$tissue %in% c("blood_emily","blood_transplant") & unique_somatic$donor %in% c(blood_emily_donor_list,blood_transplant_donor_list))])) / length(c(blood_emily_donor_list,blood_transplant_donor_list))
  prop_donor[i,"prop_colon"] <- length(unique(unique_somatic$donor[which(unique_somatic$pos == i & unique_somatic$tissue %in% c("colon","colon_ibd") & unique_somatic$donor %in% c(colon_donor_list,colon_ibd_donor_list))])) / length(c(colon_donor_list,colon_ibd_donor_list))
}

library("circlize")
  circos.par(start.degree = 90, gap.degree = 0, gap.after = 0, cell.padding = c(0,0), track.height = 0.1)
  circos.initialize(factors = c("a"), xlim = c(0,16569))
  circos.track(ylim = c(0,1), bg.border = "white")
  circos.rect(xleft = 0, xright = 16569, ybottom = 0.9, ytop = 1, col = "black")
  circos.rect(xleft = mtDNA_features$start_pos, xright = mtDNA_features$end_pos, ybottom = 0, ytop = 0.9, col = mtDNA_features$colour)
  circos.axis(major.at = seq(from = 0, to = 16589, by = 1000))
  circos.track(ylim = c(0,1), bg.border = "white")
  for(i in 1:16569){
    circos.rect(xleft = prop_donor$pos[i], xright = prop_donor$pos[i], ybottom = 0, ytop = prop_donor$prop_liver[i], col = "brown", border = "brown")
  }
  circos.track(ylim = c(0,1), bg.border = "white")
  for(i in 1:16569){
    circos.rect(xleft = prop_donor$pos[i], xright = prop_donor$pos[i], ybottom = 0, ytop = prop_donor$prop_blood[i], col = "red", border = "red")
  }
  circos.track(ylim = c(0,1), bg.border = "white")
  for(i in 1:16569){
    circos.rect(xleft = prop_donor$pos[i], xright = prop_donor$pos[i], ybottom = 0, ytop = prop_donor$prop_colon[i], col = "blue", border = "blue")
  }
  dev.copy2pdf(file = "circular_prop_liver_blood_colon_plot.pdf",out.type = "pdf")
  circos.clear()

mtDNA_coverage <- read.table("../../whole_genome_and_mtDNA_coverage.csv", sep = ",", stringsAsFactors = F, header = T)
mtDNA_coverage <- mtDNA_coverage[which(!(is.na(mtDNA_coverage$mtDNAgenomes))),]
mtDNA_coverage <- mtDNA_coverage[order(mtDNA_coverage$Tissue,mtDNA_coverage$mtDNAgenomes),]
mtDNA_coverage$Index <- 1:nrow(mtDNA_coverage)

liver_coverage <- mtDNA_coverage[which(mtDNA_coverage$Tissue == "liver"),]

liver_coverage$sub_60_T_C <- 0
liver_coverage$sub_60_T_C[which(liver_coverage$Sample %in% all_tissue_somatic$sample[which(all_tissue_somatic$tissue == "liver" & all_tissue_somatic$pos == 60 & all_tissue_somatic$ref == "T" & all_tissue_somatic$mut == "C")])] <- 1

liver_coverage$sub_72_T_C <- 0
liver_coverage$sub_72_T_C[which(liver_coverage$Sample %in% all_tissue_somatic$sample[which(all_tissue_somatic$tissue == "liver" & all_tissue_somatic$pos == 72 & all_tissue_somatic$ref == "T" & all_tissue_somatic$mut == "C")])] <- 1

liver_coverage$sub_94_G_A <- 0
liver_coverage$sub_94_G_A[which(liver_coverage$Sample %in% all_tissue_somatic$sample[which(all_tissue_somatic$tissue == "liver" & all_tissue_somatic$pos == 94 & all_tissue_somatic$ref == "G" & all_tissue_somatic$mut == "A")])] <- 1



liver_coverage$indel_CSB_2 <- 0
liver_coverage$indel_CSB_2[which(liver_coverage$Sample %in% c("PD37231b_lo013", "PD37231b_lo010", "PD37231b_lo028", "PD37231b_lo038", "PD37232b_lo020",
                                                              "PD37914b_lo0021", "PD37231b_lo005", "PD37231b_lo015", "PD37234b_lo009", "PD37234b_lo020",
                                                              "PD37914b_lo0021", "PD37917a", "PD37111b_lo0056", "PD37234b_lo019", "PD37107b_lo006",
                                                              "PD37107b_lo009", "PD37231b_lo012", "PD37231b_lo014", "PD37231b_lo019","PD37231b_lo020",
                                                              "PD37231b_lo024", "PD37234b_lo006", "PD37234b_lo024", "PD37234b_lo025", "PD37234b_lo031",
                                                              "PD37917a", "PD37115b", "PD37115b_lo002", "PD37115b_lo005", "PD37115b_lo006", "PD37115b_lo010",
                                                              "PD37115b_lo011", "PD37115b_lo012", "PD37115b_lo014", "PD37115b_lo015", "PD37115b_lo018",
                                                              "PD37115b_lo020", "PD37115b_lo021", "PD37115b_lo022", "PD37115b_lo025", "PD37115b_lo026",
                                                              "PD37115b_lo027", "PD37115b_lo028", "PD37115b_lo029", "PD37115b_lo030", "PD37115b_lo031",
                                                              "PD37115b_lo032", "PD37115b_lo033", "PD37115b_lo034", "PD37115b_lo035", "PD37115b_lo036",
                                                              "PD37115b_lo037", "PD37115b_lo038", "PD37115b_lo044", "PD37230b_lo019", "PD37230b_lo024",
                                                              "PD37107b_lo034", "PD37231b_lo019", "PD37231b_lo026", "PD37231b_lo043", "PD37232b_lo002",
                                                              "PD37111b_lo0056", "PD37111b_lo006", "PD37115a", "PD37115a_lo002", "PD37115a_lo003",
                                                              "PD37115a_lo004", "PD37230b_lo015", "PD37230b_lo038", "PD37232b_lo043", "PD37232b_lo047",
                                                              "PD37917b_lo0021", "PD37231b_lo029", "PD37107b_lo023", "PD37111a", "PD37111b_lo0056",
                                                              "PD37230b_lo035", "PD37230b_lo042", "PD37231b_lo001", "PD37231b_lo003", "PD37231b_lo005",
                                                              "PD37231b_lo016", "PD37231b_lo023", "PD37231b_lo035", "PD37234b_lo027", "PD37234b_lo032",
                                                              "PD37107b_lo009", "PD37111b_lo0056"))] <- 1

liver_coverage$sub_316_to_319 <- 0
liver_coverage$sub_316_to_319[which(liver_coverage$Sample %in% c("PD37231b_lo041", "PD37234b_lo019", "PD37917a", "PD37111a","PD37234b_lo019",
                                                                 "PD37234b_lo025", "PD37231b_lo012", "PD37231b_lo013", "PD37231b_lo014", "PD37231b_lo019",
                                                                 "PD37231b_lo020", "PD37231b_lo024", "PD37231b_lo026", "PD37231b_lo041", "PD37234b_lo006",
                                                                 "PD37234b_lo024", "PD37234b_lo031", "PD37917a"))] <- 1
liver_coverage$indel_CSB_3 <- 0
liver_coverage$indel_CSB_3[which(liver_coverage$Sample %in% c("PD36717b_lo007", "PD36717b_lo020", "PD37230b_lo009", "PD37230b_lo010"))] <- 1

liver_coverage$sub_CSB_3 <- 0
liver_coverage$sub_CSB_3[which(liver_coverage$Sample %in% c("PD37230b_lo043"))] <- 1

liver_coverage$ins_432_A <- 0
liver_coverage$ins_432_A[which(liver_coverage$Sample %in% all_tissue_somatic$sample[which(all_tissue_somatic$tissue == "liver" & all_tissue_somatic$pos == 432 & all_tissue_somatic$ref == "A" & all_tissue_somatic$mut == "INS")])] <- 1

liver_coverage$combined_call <- 0
for(i in 1:nrow(liver_coverage)){
  if(liver_coverage$sub_60_T_C[i] == 1 | liver_coverage$sub_72_T_C[i] == 1 | liver_coverage$sub_94_G_A[i] == 1 | liver_coverage$indel_CSB_2[i] == 1 | liver_coverage$sub_316_to_319[i] == 1 | liver_coverage$indel_CSB_3[i] == 1 | liver_coverage$sub_CSB_3[i] == 1 | liver_coverage$ins_432_A[i] == 1){
    liver_coverage$combined_call[i] <- 1
  } 
}

for(i in 1:nrow(liver_coverage)){
  liver_coverage$number_muts[i] <- length(which(all_tissue_somatic$sample == liver_coverage$Sample[i]))
}


liver_metadata <- read.table("/Users/al28/Mounts/cancer/mtDNA/project_metadata/liver/liver.patients.info.for.andrew.csv", sep = ",", stringsAsFactors = F, header = T)
for(i in 1:nrow(liver_metadata)){
  liver_metadata$total_samples[i] <- length(which(substr(liver_coverage$Sample,1,7) == substr(liver_metadata$PD.ID[i],1,7)))
  liver_metadata$normal_samples[i] <- length(which(substr(liver_coverage$Sample,1,8) == paste0(substr(liver_metadata$PD.ID[i],1,7),"b")))
  liver_metadata$normal_mt_mut[i] <- length(which(substr(liver_coverage$Sample,1,8) == paste0(substr(liver_metadata$PD.ID[i],1,7),"b") & liver_coverage$combined_call == 1))
  liver_metadata$normal_prop_mutated[i] <- liver_metadata$normal_mt_mut[i] / liver_metadata$normal_samples[i]
  liver_metadata$cancer_samples[i] <- length(which(substr(liver_coverage$Sample,1,8) == paste0(substr(liver_metadata$PD.ID[i],1,7),"a")))
  liver_metadata$cancer_mt_mut[i] <- length(which(substr(liver_coverage$Sample,1,8) == paste0(substr(liver_metadata$PD.ID[i],1,7),"a") & liver_coverage$combined_call == 1))
  liver_metadata$cancer_prop_mutated[i] <- liver_metadata$cancer_mt_mut[i] / liver_metadata$cancer_samples[i]
  liver_metadata$median_mtDNA_genomes[i] <- median(liver_coverage$mtDNAgenomes[which(substr(liver_coverage$Sample,1,7) == substr(liver_metadata$PD.ID[i],1,7))])
  liver_metadata$mean_mut_number[i] <- mean(liver_coverage$number_muts[which(substr(liver_coverage$Sample,1,7) == substr(liver_metadata$PD.ID[i],1,7))])
}

liver_metadata <- liver_metadata[order(liver_metadata$Aetiology, liver_metadata$median_mtDNA_genomes),]
liver_metadata <- liver_metadata[!(is.na(liver_metadata$median_mtDNA_genomes)),]

for(i in 1:nrow(liver_coverage)){
  liver_coverage$donor_id[i] <- which(substr(liver_metadata$PD.ID,1,7) == substr(liver_coverage$Sample[i],1,7))
}

write.table(liver_coverage, file = "/Users/al28/Mounts/cancer/mtDNA/Analysis/liver_coverage_and_muts_per_sample.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(liver_metadata, file = "/Users/al28/Mounts/cancer/mtDNA/Analysis/liver_coverage_ane_muts_per_donor.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

ggplot(liver_coverage) +
  geom_bar(aes(x = Index, y = mtDNAgenomes, fill = as.character(liver_coverage$combined_call)), stat = "identity") +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_fill_manual(values = c("black","cyan")) +
  theme_bw() +
  labs(fill = "Control region \nmtDNA mutation") +
  theme(panel.grid = element_blank())
ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/liver_coverage_mtDNA_mutation_waterfall_plot.png", width = 15, height = 5)

ggplot() +
  geom_jitter(data = liver_coverage, aes(x = donor_id, y = mtDNAgenomes, color = paste(substr(liver_coverage$Sample,8,8),as.character(liver_coverage$combined_call),sep = "_")), width = 0.2, size = 0.5) +
  theme_bw() +
  geom_vline(xintercept = 10.5) +
  geom_vline(xintercept = 29.5) +
  geom_text(aes(x = 5, y = 4600, label = "Alcoholism")) +
  geom_text(aes(x = 20, y = 4600, label = "NAFLD")) +
  geom_text(aes(x = 32, y = 4600, label = "Normal")) +
  scale_x_continuous(breaks = c(1:34), labels = substr(liver_metadata$PD.ID,1,7), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(0,4700), expand = c(0,0)) +
  geom_text(aes(x = )) +
  scale_color_manual(values = c("black","red","deepskyblue","forestgreen"), labels = c("Cancer - no mut", "Cancer - mut", "Normal - no mut", "Normal - mut")) +
  labs(color = "Control region \nmtDNA mutation", x = "Donor") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/liver_mtDNA_coverage_and_mutations_per_donor.pdf", width = 15, height = 5)


ggplot()