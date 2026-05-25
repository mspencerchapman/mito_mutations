setwd("/Users/al28/Mounts/cancer/mtDNA")

coverage <- read.table("whole_genome_coverage_pileup_and_bedtools_annotated.csv", sep = ",", stringsAsFactors = F, header = T)

colon_coverage <- coverage[which(coverage$Tissue %in% c("colon","colon_ibd")),]

colon_coverage <- colon_coverage[order(colon_coverage$Tissue,colon_coverage$pileup_mtDNA_genomes),]
colon_coverage$pileup_index <- 1:nrow(colon_coverage)

ggplot(colon_coverage) +
  geom_point(aes(x = pileup_index, y = pileup_mtDNA_genomes, color = Tissue)) +
  geom_segment(aes(x = 1, xend = length(which(Tissue == "colon")), y = median(pileup_mtDNA_genomes[which(Tissue == "colon")]), yend = median(pileup_mtDNA_genomes[which(Tissue == "colon")])), size = 1) +
  geom_segment(aes(x = length(which(Tissue == "colon")) + 1, xend = nrow(colon_coverage), y = median(pileup_mtDNA_genomes[which(Tissue == "colon_ibd")]), yend = median(pileup_mtDNA_genomes[which(Tissue == "colon_ibd")])), size = 1) +
  theme_bw()
ggsave("Analysis/bedtools_plots/colon_vs_colon_ibd_pileup_mtDNA_genomes.pdf", width = 10, height = 10)

colon_coverage <- colon_coverage[order(colon_coverage$Tissue,colon_coverage$bedtools_mtDNA_genomes),]
colon_coverage$bedtools_index <- 1:nrow(colon_coverage)

ggplot(colon_coverage) +
  geom_point(aes(x = bedtools_index, y = bedtools_mtDNA_genomes, color = Tissue)) +
  geom_segment(aes(x = 1, xend = length(which(Tissue == "colon")), y = median(bedtools_mtDNA_genomes[which(Tissue == "colon")]), yend = median(bedtools_mtDNA_genomes[which(Tissue == "colon")])), size = 1) +
  geom_segment(aes(x = length(which(Tissue == "colon")) + 1, xend = nrow(colon_coverage), y = median(bedtools_mtDNA_genomes[which(Tissue == "colon_ibd")]), yend = median(bedtools_mtDNA_genomes[which(Tissue == "colon_ibd")])), size = 1) +
  theme_bw()
ggsave("Analysis/bedtools_plots/colon_vs_colon_ibd_bedtools_mtDNA_genomes.pdf", width = 10, height = 10)

pcawg_coverage <- read.table("pcawg_mtdna_copy_number.csv", sep = ",", stringsAsFactors = F, header = T)
pcawg_coverage <- pcawg_coverage[which(pcawg_coverage$combined_category == "colon"),]
pcawg_coverage <- pcawg_coverage[order(pcawg_coverage$mt_copy_number),]
pcawg_coverage$index <- 1:nrow(pcawg_coverage)

ggplot(pcawg_coverage) +
  geom_point(aes(x = index, y = mt_copy_number)) +
  geom_segment(aes(x = 1, xend = nrow(pcawg_coverage), y = median(pcawg_coverage$mt_copy_number), yend = median(pcawg_coverage$mt_copy_number)), size = 1) +
  scale_y_continuous(limits = c(0,3500)) +
  theme_bw()
