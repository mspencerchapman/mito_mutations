library("tidyverse")

setwd("/lustre/scratch125/casm/team268im/al28/mtDNA/")

whole_genome_coverage <- read.csv("/lustre/scratch125/casm/team268im/al28/mtDNA/whole_genome_coverage_pileup_and_bedtools_annotated.csv", stringsAsFactors = F, header = T)

for(i in 1:nrow(whole_genome_coverage)){
  if(is.na(whole_genome_coverage$pileup_mtDNA_coverage[i])){
    if(file.exists(paste0("Data/",whole_genome_coverage$Tissue[i],"/",whole_genome_coverage$Sample[i],"_MT_count.csv"))){
      sample_cov <- read.table(paste0("Data/",whole_genome_coverage$Tissue[i],"/",whole_genome_coverage$Sample[i],"_MT_count.csv"), sep = ",", header = T)
      whole_genome_coverage$pileup_mtDNA_coverage[i] <- sum(sample_cov) / 16569
      whole_genome_coverage$pileup_mtDNA_genomes[i] <- 2 * (whole_genome_coverage$pileup_mtDNA_coverage[i] / whole_genome_coverage$SeqX[i])
    }
  }
  if(is.na(whole_genome_coverage$bedtools_autosome_coverage[i])){
    if(file.exists(paste0("Analysis/bedtools_coverage/",whole_genome_coverage$Tissue[i],"/",whole_genome_coverage$Sample[i],".genomecov.txt.summary"))){
      bedtools_summary <- read.table(paste0("Analysis/bedtools_coverage/",whole_genome_coverage$Tissue[i],"/",whole_genome_coverage$Sample[i],".genomecov.txt.summary"), sep = "\t", header = F)
      bedtools_summary <- t(bedtools_summary)
      whole_genome_coverage[i,8:34] <- bedtools_summary
    }
  }
  if(is.na(whole_genome_coverage$verify_bam_id_freemix[i])){
    if(file.exists(paste0("Analysis/verify_bam_id/",whole_genome_coverage$Tissue[i],"/",whole_genome_coverage$Sample[i],".verifyBamID.tsv.selfSM"))){
      verify_bam_id <- read.table(paste0("Analysis/verify_bam_id/",whole_genome_coverage$Tissue[i],"/",whole_genome_coverage$Sample[i],".verifyBamID.tsv.selfSM"), header = T, stringsAsFactors = F, comment.char = "")
      whole_genome_coverage$verify_bam_id_freemix[i] <- verify_bam_id$FREEMIX
    }
  }
  if(is.na(whole_genome_coverage$haplocheck_contamination[i])){
    if(file.exists(paste0("Analysis/haplocheck/",whole_genome_coverage$Tissue[i],"/",whole_genome_coverage$Sample[i],".haplocheck"))){
      haplocheck <- read.table(paste0("Analysis/haplocheck/",whole_genome_coverage$Tissue[i],"/",whole_genome_coverage$Sample[i],".haplocheck"), header = T, stringsAsFactors = F)
      if(haplocheck$Contamination.Level == "ND"){
        whole_genome_coverage$haplocheck_contamination[i] <- 0
      }else{
        whole_genome_coverage$haplocheck_contamination[i] <- haplocheck$Contamination.Level
      }
    }
  }
  if(is.na(whole_genome_coverage$ascat_normal_contamination[i])){
    if(file.exists(paste0("/Users/al28/Mounts/nst_live_links/",whole_genome_coverage$Study[i],"/",whole_genome_coverage$Sample[i],"/",whole_genome_coverage$Sample[i],".samplestatistics.txt"))){
      ascat <- as.data.frame(t(read.table(paste0("/Users/al28/Mounts/nst_live_links/",whole_genome_coverage$Study[i],"/",whole_genome_coverage$Sample[i],"/",whole_genome_coverage$Sample[i],".samplestatistics.txt"), row.names = 1, stringsAsFactors = F)))
      whole_genome_coverage$ascat_normal_contamination[i] <- as.numeric(ascat$NormalContamination)
      whole_genome_coverage$ascat_goodness_of_fit[i] <- as.numeric(ascat$goodnessOfFit)
      whole_genome_coverage$ascat_ploidy[i] <- as.numeric(ascat$Ploidy)
    }
  }
  if(!(is.na(whole_genome_coverage$bedtools_autosome_coverage[i])) & !(is.na(whole_genome_coverage$ascat_ploidy[i])) & is.na(whole_genome_coverage$bedtools_mtDNA_genomes_ascat_corrected[i])){
    whole_genome_coverage$bedtools_mtDNA_genomes_ascat_corrected[i] <- (as.numeric(whole_genome_coverage$bedtools_mtDNA_coverage[i]) / as.numeric(whole_genome_coverage$bedtools_autosome_coverage[i])) * (((1-as.numeric(whole_genome_coverage$ascat_normal_contamination[i])) * as.numeric(whole_genome_coverage$ascat_ploidy[i])) + (2 * as.numeric(whole_genome_coverage$ascat_normal_contamination[i])))
  }
}

whole_genome_coverage <- whole_genome_coverage[order(whole_genome_coverage$Tissue, whole_genome_coverage$Sample),]

write.table(x = whole_genome_coverage, file = "/lustre/scratch125/casm/team268im/al28/mtDNA/whole_genome_coverage_pileup_and_bedtools_annotated.csv", sep = ",", col.names = T, row.names = F, quote = F)

# length(which(is.na(whole_genome_coverage$pileup_mtDNA_coverage)))
# length(which(is.na(whole_genome_coverage$bedtools_autosome_coverage)))
# length(which(is.na(whole_genome_coverage$verify_bam_id_freemix)))
# length(which(is.na(whole_genome_coverage$haplocheck_contamination)))
# length(which(is.na(whole_genome_coverage$ascat_ploidy)))
# 
# lower_bedtools_cov <- whole_genome_coverage[which(!(is.na(whole_genome_coverage$pileup_mtDNA_coverage)) & !(is.na(whole_genome_coverage$bedtools_autosome_coverage))),][which(whole_genome_coverage$pileup_mtDNA_genomes[which(!(is.na(whole_genome_coverage$pileup_mtDNA_coverage)) & !(is.na(whole_genome_coverage$bedtools_autosome_coverage)))] > whole_genome_coverage$bedtools_mtDNA_genomes[which(!(is.na(whole_genome_coverage$pileup_mtDNA_coverage)) & !(is.na(whole_genome_coverage$bedtools_autosome_coverage)))]),c("Tissue","Study","Sample","pileup_mtDNA_genomes","bedtools_mtDNA_genomes")]
# 
# ggplot(whole_genome_coverage[which(!(is.na(whole_genome_coverage$pileup_mtDNA_coverage)) & !(is.na(whole_genome_coverage$bedtools_autosome_coverage))),]) +
#   geom_point(aes(x = pileup_mtDNA_genomes, y = bedtools_mtDNA_genomes)) +
#   geom_abline(slope = 1, intercept = 0, lty = 2, color = "red") +
#   scale_x_continuous(limits = c(0,max(c(whole_genome_coverage$pileup_mtDNA_genomes[which(!(is.na(whole_genome_coverage$pileup_mtDNA_coverage)) & !(is.na(whole_genome_coverage$bedtools_autosome_coverage)))],whole_genome_coverage$bedtools_mtDNA_genomes[which(!(is.na(whole_genome_coverage$pileup_mtDNA_coverage)) & !(is.na(whole_genome_coverage$bedtools_autosome_coverage)))]), na.rm = T))) +
#   scale_y_continuous(limits = c(0,max(c(whole_genome_coverage$pileup_mtDNA_genomes[which(!(is.na(whole_genome_coverage$pileup_mtDNA_coverage)) & !(is.na(whole_genome_coverage$bedtools_autosome_coverage)))],whole_genome_coverage$bedtools_mtDNA_genomes[which(!(is.na(whole_genome_coverage$pileup_mtDNA_coverage)) & !(is.na(whole_genome_coverage$bedtools_autosome_coverage)))]), na.rm = T))) +
#   theme_bw()
# ggsave("/Users/al28/Mounts/cancer/mtDNA/Analysis/bedtools_plots/bed_tools_cov_vs_pile_up_cov.pdf", width = 10, height = 10)
# 
# ggplot() +
#   geom_point(data = whole_genome_coverage[which(whole_genome_coverage$Tissue != "liver"),], aes(x = bedtools_chr1_coverage, y = bedtools_chr22_coverage), color = "black") +
#   geom_point(data = whole_genome_coverage[which(whole_genome_coverage$Tissue == "liver"),], aes(x = bedtools_chr1_coverage, y = bedtools_chr22_coverage), color = "red") +
#   geom_abline(slope = 1, intercept = 0, lty = 2, color = "red") +
#   theme_bw()
# 
# ggplot() +
#   geom_point(data = whole_genome_coverage[which(whole_genome_coverage$Tissue != "testes_tumour"),], aes(x = (bedtools_chrX_coverage / bedtools_autosome_coverage), y = (bedtools_chrY_coverage / bedtools_autosome_coverage)), color = "black") +
#   geom_point(data = whole_genome_coverage[which(whole_genome_coverage$Tissue == "testes_tumour"),], aes(x = (bedtools_chrX_coverage / bedtools_autosome_coverage), y = (bedtools_chrY_coverage / bedtools_autosome_coverage)), color = "red") +
#   scale_x_continuous(limits = c(0,2)) +
#   scale_y_continuous(limits = c(0,2)) +
#   theme_bw()
# 
# whole_genome_coverage$Sample[which((whole_genome_coverage$bedtools_chr1_coverage + whole_genome_coverage$bedtools_chr2_coverage) / (whole_genome_coverage$bedtools_chr21_coverage + whole_genome_coverage$bedtools_chr22_coverage) > 2)]
# 
# contaminated <- whole_genome_coverage[which(whole_genome_coverage$verify_bam_id_freemix > 0.01),]
# excluded_samples <- read.table("excluded_samples.csv", sep = ",", stringsAsFactors = F, header = T)
# 
# whole_genome_coverage_ratio <- whole_genome_coverage
# whole_genome_coverage_ratio$chr1_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr1_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr2_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr2_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr3_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr3_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr4_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr4_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr5_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr5_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr6_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr6_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr7_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr7_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr8_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr8_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr9_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr9_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr10_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr10_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr11_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr11_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr12_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr12_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr13_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr13_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr14_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr14_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr15_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr15_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr16_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr16_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr17_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr17_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr18_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr18_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr19_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr19_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr20_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr20_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr21_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr21_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$chr22_residual <- (abs(whole_genome_coverage_ratio$bedtools_chr22_coverage - whole_genome_coverage_ratio$bedtools_autosome_coverage)) / whole_genome_coverage_ratio$bedtools_autosome_coverage
# whole_genome_coverage_ratio$autosome_residuals <- whole_genome_coverage_ratio$chr1_residual + whole_genome_coverage_ratio$chr2_residual + whole_genome_coverage_ratio$chr3_residual + whole_genome_coverage_ratio$chr4_residual + whole_genome_coverage_ratio$chr5_residual + whole_genome_coverage_ratio$chr6_residual + whole_genome_coverage_ratio$chr7_residual + whole_genome_coverage_ratio$chr8_residual + whole_genome_coverage_ratio$chr9_residual + whole_genome_coverage_ratio$chr10_residual + whole_genome_coverage_ratio$chr11_residual + whole_genome_coverage_ratio$chr12_residual + whole_genome_coverage_ratio$chr13_residual + whole_genome_coverage_ratio$chr14_residual + whole_genome_coverage_ratio$chr15_residual + whole_genome_coverage_ratio$chr16_residual+ whole_genome_coverage_ratio$chr17_residual + whole_genome_coverage_ratio$chr18_residual + whole_genome_coverage_ratio$chr19_residual + whole_genome_coverage_ratio$chr20_residual + whole_genome_coverage_ratio$chr21_residual + whole_genome_coverage_ratio$chr22_residual
