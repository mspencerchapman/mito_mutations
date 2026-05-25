library("dplyr")
library("GenomicRanges")
library("Rsamtools")
library("MASS")

setwd("/Users/al28/Mounts/cancer/mtDNA/Analysis/shearwater_stringent_exclusion/")

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

for(i in 1:length(project_list)){
  setwd(project_list[i])
  donor_list <- list.files(".")
  donor_list <- donor_list[grepl("_shearwater_calls.txt", donor_list)]
  
  if(length(donor_list) != 0){
    
    trinuc_df_donor <- as.data.frame(array(data = NA,dim = c(length(donor_list),192)))
    rownames(trinuc_df_donor) <- paste(project_list[i],gsub(pattern = "_shearwater_calls.txt",replacement = "",x = donor_list),sep = "_")
    colnames(trinuc_df_donor) <- full_vec
    
    for(j in 1:length(donor_list)){
      donor_calls <- read.table(donor_list[j], sep = "\t", stringsAsFactors = F, header = T)
      donor_calls <- donor_calls[which(donor_calls$vaf >= 0.01),]
      donor_calls <- donor_calls[which(!(duplicated(donor_calls$mut_site))),]
      
      donor_calls <- donor_calls[which(donor_calls$mut %in% c("A","C","G","T")),]
      donor_calls <- donor_calls[which(donor_calls$pos %in% coding_region),]
      
      trinuc_vec <- vector(length = 192)
      
      mutations <- donor_calls[,c("chr","pos","ref","mut")]
      mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(mutations$pos-1, mutations$pos+1))))
      mutations$sub = paste(mutations$ref,mutations$mut,sep = ">")
      
      freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref,1,1),substr(mutations$trinuc_ref,3,3),sep="-"),sep=","))
      freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
      
      trinuc_df_donor[j,] <- freqs_full 
    }
    if(i == 1){
      trinuc_df_combined <- trinuc_df_donor
    } else {
      trinuc_df_combined <- rbind(trinuc_df_combined,trinuc_df_donor)
    }
  }
  setwd("..")
}

write.table(trinuc_df_combined, "/Users/al28/Mounts/cancer/mtDNA/Analysis/hdp/shearwater_stringent_exclusion/HDPinput_PerDonorCollapsed_VAFabove0.01_CodingOnly.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
