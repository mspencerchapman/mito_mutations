library("dplyr")
library("GenomicRanges")
library("Rsamtools")
library("MASS")

setwd("/Users/al28/Mounts/cancer/mtDNA/Analysis/shearwater_PCAWG_metanorm_normal_panel/")

project_list <- list.files(".")

coding_region <- 577:16023
d_loop_region <- c(1:576,16024:16569)

genomeFile <- ("/Users/al28/Mounts/cancer/bed_ref/hs37d5.fa")

sub_vec_for = c("C>A","C>G","C>T","T>A","T>C","T>G")
ctx_vec_for = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec_for = paste(rep(sub_vec_for,each=16),rep(ctx_vec_for,times=6),sep=",")

sub_vec_rev = c("G>T","G>C","G>A","A>T","A>G","A>C")
ctx_vec_rev = paste(rep(c("T","G","C","A"),times=4),rep(c("T","G","C","A"),each=4),sep="-")
full_vec_rev = paste(rep(sub_vec_rev,each=16),rep(ctx_vec_rev,times=6),sep=",")

full_vec <- c(full_vec_for,full_vec_rev)

shearwater_calls <- read.table("combined_PCAWG_calls_all_samples_no_cutoff.tsv", sep = "\t", stringsAsFactors = F, header = T)

tissue_list <- unique(shearwater_calls$tissue)

trinuc_df_tissue <- as.data.frame(array(data = NA,dim = c(length(tissue_list),192)))
rownames(trinuc_df_tissue) <- gsub(pattern = "PCAWG_", replacement = "", x = tissue_list)
colnames(trinuc_df_tissue) <- full_vec

for(j in 1:length(tissue_list)){
  tissue_calls <- shearwater_calls[which(shearwater_calls$tissue == tissue_list[j]),]
  tissue_calls <- tissue_calls[which(tissue_calls$vaf >= 0.02),]
  tissue_calls <- tissue_calls[which(tissue_calls$type == "t"),]
  
  tissue_calls <- tissue_calls[which(tissue_calls$mut %in% c("A","C","G","T")),]
  tissue_calls <- tissue_calls[which(tissue_calls$pos %in% coding_region),]
  
  trinuc_vec <- vector(length = 192)
  
  mutations <- tissue_calls[,c("chr","pos","ref","mut")]
  mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(mutations$pos-1, mutations$pos+1))))
  mutations$sub = paste(mutations$ref,mutations$mut,sep = ">")
  
  freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref,1,1),substr(mutations$trinuc_ref,3,3),sep="-"),sep=","))
  freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
  
  trinuc_df_tissue[j,] <- freqs_full 
}

write.table(trinuc_df_tissue, "/Users/al28/Mounts/cancer/mtDNA/Analysis/hdp/pcawg_shearwater_tumour_only/HDPinput_TumourOnly_VAFabove0.02_CodingOnly.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
