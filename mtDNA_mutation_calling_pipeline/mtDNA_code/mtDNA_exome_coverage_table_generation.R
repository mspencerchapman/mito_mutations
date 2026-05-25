library("tidyverse")

setwd("/Users/al28/Mounts/cancer/mtDNA/")

exome_coverage <- read.csv("/Users/al28/Mounts/cancer/mtDNA/exome_coverage_pileup_and_bedtools_annotated.csv", stringsAsFactors = F, header = T)

for(i in 1:nrow(exome_coverage)){
  if(is.na(exome_coverage$pileup_mtDNA_coverage[i])){
    if(file.exists(paste0("Data/",exome_coverage$Tissue[i],"/",exome_coverage$Sample[i],"_MT_count.csv"))){
      sample_cov <- read.table(paste0("Data/",exome_coverage$Tissue[i],"/",exome_coverage$Sample[i],"_MT_count.csv"), sep = ",", header = T)
      exome_coverage$pileup_mtDNA_coverage[i] <- sum(sample_cov) / 16569
      }
  }
  if(is.na(exome_coverage$bedtools_autosome_coverage[i])){
    if(file.exists(paste0("Analysis/bedtools_coverage/",exome_coverage$Tissue[i],"/",exome_coverage$Sample[i],".genomecov.txt.summary"))){
      bedtools_summary <- read.table(paste0("Analysis/bedtools_coverage/",exome_coverage$Tissue[i],"/",exome_coverage$Sample[i],".genomecov.txt.summary"), sep = "\t", header = F)
      bedtools_summary <- t(bedtools_summary)
      exome_coverage[i,7:32] <- bedtools_summary
    }
  }
  if(is.na(exome_coverage$verify_bam_id_freemix[i])){
    if(file.exists(paste0("Analysis/verify_bam_id/",exome_coverage$Tissue[i],"/",exome_coverage$Sample[i],".verifyBamID.tsv.selfSM"))){
      verify_bam_id <- read.table(paste0("Analysis/verify_bam_id/",exome_coverage$Tissue[i],"/",exome_coverage$Sample[i],".verifyBamID.tsv.selfSM"), header = T, stringsAsFactors = F, comment.char = "")
      exome_coverage$verify_bam_id_freemix[i] <- verify_bam_id$FREEMIX
    }
  }
  if(is.na(exome_coverage$haplocheck_contamination[i])){
    if(file.exists(paste0("Analysis/haplocheck/",exome_coverage$Tissue[i],"/",exome_coverage$Sample[i],".haplocheck"))){
      haplocheck <- read.table(paste0("Analysis/haplocheck/",exome_coverage$Tissue[i],"/",exome_coverage$Sample[i],".haplocheck"), header = T, stringsAsFactors = F)
      if(haplocheck$Contamination.Level == "ND"){
        exome_coverage$haplocheck_contamination[i] <- 0
      }else{
        exome_coverage$haplocheck_contamination[i] <- haplocheck$Contamination.Level
      }
    }
  }
}

exome_coverage <- exome_coverage[order(exome_coverage$Tissue, exome_coverage$Sample),]

write.table(x = exome_coverage, file = "/Users/al28/Mounts/cancer/mtDNA/exome_coverage_pileup_and_bedtools_annotated.csv", sep = ",", col.names = T, row.names = F, quote = F)

length(which(is.na(exome_coverage$pileup_mtDNA_coverage)))
length(which(is.na(exome_coverage$bedtools_autosome_coverage)))
length(which(is.na(exome_coverage$verify_bam_id_freemix)))
length(which(is.na(exome_coverage$haplocheck_contamination)))

synovium_excluded_samples <- c("PD43957b_lo0007","PD43957b_lo0008","PD43957b_lo0011","PD43957b_lo0014","PD43957b_lo0015","PD43957b_lo0017",
                               "PD43957b_lo0018","PD43957b_lo0026","PD43957b_lo0032","PD43957b_lo0037","PD43957b_lo0038","PD43957b_lo0039",
                               "PD43957b_lo0043","PD43957b_lo0050","PD43957b_lo0057","PD43958b_lo0004","PD43958b_lo0009","PD43958b_lo0012",
                               "PD43958b_lo0016","PD43958b_lo0017","PD43959b_lo0001","PD43959b_lo0003","PD43959b_lo0004","PD43959b_lo0006",
                               "PD43959b_lo0007","PD43959b_lo0009","PD43959b_lo0015","PD43959b_lo0021","PD43959b_lo0025","PD43959b_lo0028",
                               "PD43959b_lo0030","PD43959b_lo0031","PD43959b_lo0033","PD43959b_lo0040","PD43959b_lo0042","PD43959b_lo0044",
                               "PD43959b_lo0045","PD43959b_lo0046","PD43959b_lo0049","PD43959b_lo0052","PD43959b_lo0053","PD43959b_lo0060",
                               "PD43961b_lo0003","PD43961b_lo0024","PD43961b_lo0025","PD43961b_lo0028","PD43961b_lo0032","PD43961b_lo0040",
                               "PD43961b_lo0046")

bladder_excluded_samples <- c("PD42031e_lo0005","PD42031e_lo0006","PD42031e_lo0008","PD42031e_lo0013","PD42031f_lo0010","PD42031f_lo0011",
                              "PD42031f_lo0012","PD42031f_lo0013","PD42031f_lo0014","PD42032c_lo0005","PD42032c_lo0006","PD42032c_lo0007",
                              "PD42032c_lo0013","PD42032c_lo0023")

testes_excluded_samples <- c("PD40745c_lo0010","PD40745c_lo0012","PD40746c_lo0010","PD40746c_lo0012","PD40746c_lo0013","PD40746c_lo0014",
                             "PD40746c_lo0015","PD40746c_lo0019","PD40744c_lo0023","PD40745c_lo0014","PD40745c_lo0015","PD40746c_lo0009",
                             "PD40746c_lo0011","PD40746c_lo0020","PD40746c_lo0021","PD40746c_lo0024","PD42034b_lo0012","PD42034b_lo0013",
                             "PD42034b_lo0014","PD42034b_lo0015","PD42034b_lo0016","PD42034b_lo0017","PD42036b_lo0012","PD42036b_lo0013",
                             "PD42036b_lo0014","PD42036b_lo0017","PD42036b_lo0018","PD42036b_lo0019","PD42036b_lo0020","PD42036b_lo0021")

ggplot(data = exome_coverage[which(exome_coverage$Tissue %in% c("synovium_exome","testes_exome") & !(exome_coverage$Sample %in% synovium_excluded_samples) & !(exome_coverage$Sample %in% testes_excluded_samples)),], 
       aes(x = bedtools_chrY_coverage / bedtools_autosome_coverage, fill = Tissue)) +
  geom_histogram((aes(y = 0.005*..density..)), position = "identity", binwidth = 0.005, alpha = 0.5) +
 # scale_x_continuous(limits = c(0.5,1.5)) +
  theme_bw() +
  theme(legend.position = "none")

synovium_only <- exome_coverage[which(exome_coverage$Tissue == "synovium_exome"),]
synovium_only$label <- "excluded"

synovium_metadata <- read.table("/Users/al28/Desktop/All_Samples_20-04-2021.csv", sep = ",", stringsAsFactors = F, header = T)

synovium_only$label[which(synovium_only$Sample %in% synovium_metadata$SUPPLIER.SAMPLE.NAME[which(synovium_metadata$Submission_Order <= 384)])] <- "Sanger168"
synovium_only$label[which(synovium_only$Sample %in% synovium_metadata$SUPPLIER.SAMPLE.NAME[which(synovium_metadata$Submission_Order >= 384 & synovium_metadata$Submission_Order <= 528)])] <- "Uncertain"
synovium_only$label[which(synovium_only$Sample %in% synovium_metadata$SUPPLIER.SAMPLE.NAME[which(synovium_metadata$Submission_Order >= 529)])] <- "UDI"
synovium_only$label[which(synovium_only$Sample %in% synovium_excluded_samples)] <- "excluded"


ggplot(data = synovium_only, 
       aes(x = bedtools_chr11_coverage / bedtools_autosome_coverage, fill = label)) +
  geom_histogram((aes(y = 0.005*..density..)), position = "identity", binwidth = 0.005, alpha = 0.5) +
  scale_x_continuous(limits = c(0.5,1.5)) +
  theme_bw()
