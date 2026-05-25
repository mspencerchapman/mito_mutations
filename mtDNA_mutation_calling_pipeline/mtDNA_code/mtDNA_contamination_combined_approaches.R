library("tidyverse")

setwd("/Users/al28/Mounts/cancer/mtDNA/")

excluded_samples <- read.table("excluded_samples.csv", sep = ",", stringsAsFactors = F, header = T)
mtDNA_contamination <- read.table("mtDNA_contamination_estimates.tsv", sep = "\t", stringsAsFactors = F, header = T)
nuclear_contamination <- read.table("whole_genome_coverage_pileup_and_bedtools_annotated.csv", sep = ",", stringsAsFactors = F, header = T)

combined_excluded <- as.data.frame(array(data = "", dim = c(nrow(nuclear_contamination),10)))
colnames(combined_excluded) <- c("Tissue","Study","Sample","ManualSampleSwap","ManualContamination","NuclearContamination","NuclearFail","MitochondriaContamination","MitochondriaFail","CombinedExclusion")

combined_excluded$Tissue <- nuclear_contamination$Tissue
combined_excluded$Study <- nuclear_contamination$Study
combined_excluded$Sample <- nuclear_contamination$Sample

for(i in 1:nrow(excluded_samples)){
  if(excluded_samples$comment[i] == "sample_swap"){
    combined_excluded$ManualSampleSwap[which(combined_excluded$Sample == excluded_samples$sample[i] & combined_excluded$Tissue == excluded_samples$tissue[i])] <- "Y"
  }else if(excluded_samples$comment[i] == "low_level_contamination"){
    combined_excluded$ManualContamination[which(combined_excluded$Sample == excluded_samples$sample[i] & combined_excluded$Tissue == excluded_samples$tissue[i])] <- "Y"
  }
}

for(i in 1:nrow(nuclear_contamination)){
  if(is.na(nuclear_contamination$verify_bam_id_freemix[i])){
    combined_excluded$NuclearFail[which(combined_excluded$Sample == nuclear_contamination$Sample[i] & combined_excluded$Tissue == nuclear_contamination$Tissue[i])] <- "Y"
  }else if(nuclear_contamination$verify_bam_id_freemix[i] >= 0.01){
    combined_excluded$NuclearContamination[which(combined_excluded$Sample == nuclear_contamination$Sample[i] & combined_excluded$Tissue == nuclear_contamination$Tissue[i])] <- "Y"
  }
}

for(i in 1:nrow(combined_excluded)){
  if(!(combined_excluded$Sample[i] %in% mtDNA_contamination$sampleID & combined_excluded$Tissue[i] %in% mtDNA_contamination$tissue)){
    combined_excluded$MitochondriaFail[i] <- "Y"
  }else if(mtDNA_contamination$contamination[which(mtDNA_contamination$sampleID == combined_excluded$Sample[i] & mtDNA_contamination$tissue == combined_excluded$Tissue[i])] %in% c("YY","YX")){
    combined_excluded$MitochondriaContamination[i] <- "Y"
  }
}

for(i in 1:nrow(combined_excluded)){
  if(any(combined_excluded$ManualSampleSwap[i] == "Y",combined_excluded$ManualContamination[i] == "Y", combined_excluded$NuclearContamination[i] == "Y", combined_excluded$NuclearFail[i] == "Y", combined_excluded$MitochondriaContamination[i] == "Y", combined_excluded$MitochondriaFail[i] == "Y")){
    combined_excluded$CombinedExclusion[i] <- "Y"
  }
}

write.table(x = combined_excluded, file = "combined_sample_exclusion.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

mtDNA_contamination$manual_exclusion <- "N"

for(i in 1:nrow(excluded_samples)){
  mtDNA_contamination$manual_exclusion[which(mtDNA_contamination$sampleID == excluded_samples$sample[i])] <- excluded_samples$comment[i] 
}

ggplot(mtDNA_contamination) +
  geom_point(aes(x = sum_exc_top_bottom, y = second_min_diff, color = manual_exclusion)) +
  theme_bw() +
  geom_vline(xintercept = -0.05, lty = 2) +
  geom_hline(yintercept = -0.01, lty = 2) +
  labs(x = "Sum VAF difference from mean VAF in other samples across germline SNPs (exc. top and bottom SNPs)", y = "Second lowest difference from mean VAF", color = "Manual assignment")
ggsave("Analysis/contamination_analysis/manual_vs_automated_mtDNA_comparison_full_range.pdf", width = 10, height = 10)

ggplot(mtDNA_contamination) +
  geom_point(aes(x = sum_exc_top_bottom, y = second_min_diff, color = manual_exclusion), alpha = 0.8, size = 0.7) +
  theme_bw() +
  geom_vline(xintercept = -0.05, lty = 2) +
  geom_hline(yintercept = -0.01, lty = 2) +
  scale_x_continuous(limits = c(-5,0)) +
  scale_y_continuous(limits = c(-0.35,0)) +
  labs(x = "Sum VAF difference from mean VAF in other samples across germline SNPs (exc. top and bottom SNPs)", y = "Second lowest difference from mean VAF", color = "Manual assignment")
ggsave("Analysis/contamination_analysis/manual_vs_automated_mtDNA_comparison_zoom_1.pdf", width = 10, height = 10)

ggplot(mtDNA_contamination) +
  geom_point(aes(x = sum_exc_top_bottom, y = second_min_diff, color = manual_exclusion), alpha = 0.8, size = 0.7) +
  theme_bw() +
  geom_vline(xintercept = -0.05, lty = 2) +
  geom_hline(yintercept = -0.01, lty = 2) +
  scale_x_continuous(limits = c(-0.6,0)) +
  scale_y_continuous(limits = c(-0.1,0)) +
  labs(x = "Sum VAF difference from mean VAF in other samples across germline SNPs (exc. top and bottom SNPs)", y = "Second lowest difference from mean VAF", color = "Manual assignment")
ggsave("Analysis/contamination_analysis/manual_vs_automated_mtDNA_comparison_zoom_2.pdf", width = 10, height = 10)


mtDNA_contamination$nuclear_contamination <- NA

for(i in 1:nrow(mtDNA_contamination)){
  mtDNA_contamination$nuclear_contamination[i] <- nuclear_contamination$verify_bam_id_freemix[which(nuclear_contamination$Sample == mtDNA_contamination$sampleID[i] & nuclear_contamination$Tissue == mtDNA_contamination$tissue[i])]
}

mtDNA_contamination$nuclear_contamination_cutoff <- "N"
mtDNA_contamination$nuclear_contamination_cutoff[which(mtDNA_contamination$nuclear_contamination > 0.01)] <- "Y"

ggplot(mtDNA_contamination) +
  geom_point(aes(x = sum_exc_top_bottom, y = second_min_diff, color = nuclear_contamination_cutoff)) +
  theme_bw() +
  geom_vline(xintercept = -0.05, lty = 2) +
  geom_hline(yintercept = -0.01, lty = 2) +
  labs(x = "Sum VAF difference from mean VAF in other samples across germline SNPs (exc. top and bottom SNPs)", y = "Second lowest difference from mean VAF", color = "VerifyBamID > 0.01")
ggsave("Analysis/contamination_analysis/mito_vs_nuclear_mtDNA_comparison_full_range.pdf", width = 10, height = 10)


ggplot(mtDNA_contamination) +
  geom_point(aes(x = sum_exc_top_bottom, y = second_min_diff, color = nuclear_contamination_cutoff), alpha = 0.8, size = 0.7) +
  theme_bw() +
  geom_vline(xintercept = -0.05, lty = 2) +
  geom_hline(yintercept = -0.01, lty = 2) +
  scale_x_continuous(limits = c(-5,0)) +
  scale_y_continuous(limits = c(-0.35,0)) +
  labs(x = "Sum VAF difference from mean VAF in other samples across germline SNPs (exc. top and bottom SNPs)", y = "Second lowest difference from mean VAF", color = "VerifyBamID > 0.01")
ggsave("Analysis/contamination_analysis/mito_vs_nuclear_mtDNA_comparison_zoom_1.pdf", width = 10, height = 10)



nuclear_contamination$verify_neg_log <- -log(nuclear_contamination$verify_bam_id_freemix, base = 10)
nuclear_contamination <- nuclear_contamination[order(nuclear_contamination$verify_neg_log, decreasing = T),]
nuclear_contamination$index <- 1:nrow(nuclear_contamination)

nuclear_contamination_cutoff <- 0.01

ggplot(nuclear_contamination) +
  geom_point(aes(x = index, y = verify_neg_log)) +
  geom_vline(xintercept = max(nuclear_contamination$index[which(nuclear_contamination$verify_neg_log > -log(nuclear_contamination_cutoff, 10))]) + 0.5, color = "red", lty = 2) +
  theme_bw()


ggplot() +
  geom_point(data = nuclear_contamination[which(nuclear_contamination$Tissue != "gi_tps"),], aes(x = index, y = verify_neg_log), color = "black") +
  geom_point(data = nuclear_contamination[which(nuclear_contamination$Tissue == "gi_tps"),], aes(x = index, y = verify_neg_log), color = "red") +
  geom_vline(xintercept = max(nuclear_contamination$index[which(nuclear_contamination$verify_neg_log > -log(nuclear_contamination_cutoff, 10))]) + 0.5, color = "red", lty = 2) +
  theme_bw()

ggplot(nuclear_contamination) +
  geom_point(aes(x = index, y = verify_neg_log)) +
  geom_vline(xintercept = max(nuclear_contamination$index[which(nuclear_contamination$verify_neg_log > -log(nuclear_contamination_cutoff, 10))]) + 0.5, color = "red", lty = 2) +
  theme_bw() +
  scale_y_continuous(limits = c(0,5))

sample_swap_list <- excluded_samples$sample[which(excluded_samples$comment == "sample_swap")]
nuclear_con_list <- nuclear_contamination$Sample[which(nuclear_contamination$verify_bam_id_freemix >= 0.01)]
mito_con_list <- mtDNA_contamination$sampleID[which(mtDNA_contamination$contamination == "YY" | mtDNA_contamination$contamination == "YX")]

combined_exc_list <- unique(c(sample_swap_list,nuclear_con_list,mito_con_list,low_level_con_list))

new_excluded_list <- as.data.frame(array(data = NA, dim = c(length(combined_exc_list),8)))
colnames(new_excluded_list) <- c("Tissue","Study","Donor","Sample","SampleSwap","NuclearContamination","MitoContamination","Other")
