library("Rsamtools")

setwd("/Users/al28/Mounts/cancer/mtDNA/Analysis/shearwater_stringent_exclusion/")

tissues <- list.files(".")

for(i in 1:length(tissues)){
  all_files <- list.files(tissues[i])
  input_files <- all_files[grepl("shearwater_calls.txt",all_files)]
  if(length(input_files) > 0){
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
}

setwd("../shearwater_tissue_overview/")

genomeFile <- ("/Users/al28/Mounts/cancer/bed_ref/hs37d5.fa")
mtdna_trinuc_freq <- readRDS("/Users/al28/Mounts/cancer/mtDNA/trinuc_freqs_coding_dloop_heavy_light.Rds")

coding_region <- 577:16023
d_loop_region <- c(1:576,16024:16569)

trinucleotide_plot = function (mutations, file_name, analysis_type, analysis_region) {
  mutations <- unique(mutations[,c("chr","pos","ref","mut","donor")])
  mutations <- mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")),]
  
  if(analysis_region == "coding"){
    mutations <- mutations[which(mutations$pos %in% coding_region),]
  }else if(analysis_region == "d_loop"){
    mutations <- mutations[which(mutations$pos %in% d_loop_region),]
  }else if(analysis_region != "all_mtDNA"){
    mutations <- NULL
  }
  
  mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(mutations$pos-1, mutations$pos+1))))
  
  # 2. Annotating the mutation from the pyrimidine base
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  
  # 3. Counting subs
  freqs_heavy = table(paste(mutations$sub[which(mutations$ref %in% c("A","G"))],paste(substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("A","G"))],1,1),substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("A","G"))],3,3),sep="-"),sep=","))
  freqs_light = table(paste(mutations$sub[which(mutations$ref %in% c("C","T"))],paste(substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("C","T"))],1,1),substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("C","T"))],3,3),sep="-"),sep=","))
  
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_heavy_full = freqs_heavy[full_vec]; freqs_heavy_full[is.na(freqs_heavy_full)] = 0; names(freqs_heavy_full) = full_vec
  freqs_light_full = freqs_light[full_vec]; freqs_light_full[is.na(freqs_light_full)] = 0; names(freqs_light_full) = full_vec
  
  
  if(analysis_type == "obs_exp"){
    if(analysis_region == "coding"){
      heavy_base_freqs <- (mtdna_trinuc_freq["coding_heavy",] / sum(mtdna_trinuc_freq["coding_heavy",]))
      exp_heavy_counts <- sum(freqs_heavy) * as.numeric(c(rep(heavy_base_freqs[1:16]/3,times = 3),rep(heavy_base_freqs[17:32]/3,times = 3)))
      freqs_heavy_full <- freqs_heavy_full / exp_heavy_counts
      
      light_base_freqs <- (mtdna_trinuc_freq["coding_light",] / sum(mtdna_trinuc_freq["coding_light",]))
      exp_light_counts <- sum(freqs_light) * as.numeric(c(rep(light_base_freqs[1:16]/3,times = 3),rep(light_base_freqs[17:32]/3,times = 3)))
      freqs_light_full <- freqs_light_full / exp_light_counts
    }else if(analysis_region == "d_loop"){
      heavy_base_freqs <- (mtdna_trinuc_freq["d_loop_heavy",] / sum(mtdna_trinuc_freq["d_loop_heavy",]))
      exp_heavy_counts <- sum(freqs_heavy) * as.numeric(c(rep(heavy_base_freqs[1:16]/3,times = 3),rep(heavy_base_freqs[17:32]/3,times = 3)))
      freqs_heavy_full <- freqs_heavy_full / exp_heavy_counts
      
      light_base_freqs <- (mtdna_trinuc_freq["d_loop_light",] / sum(mtdna_trinuc_freq["d_loop_light",]))
      exp_light_counts <- sum(freqs_light) * as.numeric(c(rep(light_base_freqs[1:16]/3,times = 3),rep(light_base_freqs[17:32]/3,times = 3)))
      freqs_light_full <- freqs_light_full / exp_light_counts
    }else if(analysis_region == "all_mtDNA"){
      heavy_base_freqs <- colSums(mtdna_trinuc_freq[c("coding_heavy","d_loop_heavy"),]) / sum(colSums(mtdna_trinuc_freq[c("coding_heavy","d_loop_heavy"),]))
      exp_heavy_counts <- sum(freqs_heavy) * as.numeric(c(rep(heavy_base_freqs[1:16]/3,times = 3),rep(heavy_base_freqs[17:32]/3,times = 3)))
      freqs_heavy_full <- freqs_heavy_full / exp_heavy_counts
      
      light_base_freqs <- colSums(mtdna_trinuc_freq[c("coding_light","d_loop_light"),]) / sum(colSums(mtdna_trinuc_freq[c("coding_light","d_loop_light"),]))
      exp_light_counts <- sum(freqs_light) * as.numeric(c(rep(light_base_freqs[1:16]/3,times = 3),rep(light_base_freqs[17:32]/3,times = 3)))
      freqs_light_full <- freqs_light_full / exp_light_counts
    }else{
      freqs_heavy_full = NULL
      freqs_light_full = NULL
    }
  }
  
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  
  dev.new(width=10,height=4)
  colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  y_heavy = freqs_heavy_full; y_light = freqs_light_full; maxy = max(c(y_heavy,y_light))
  
  if(analysis_type == "obs_exp"){
    ylab = "Mutation frequency (Obs/Exp)"
  }else{
    ylab = "Mutation count"
  }
  
  h_heavy = barplot(y_heavy, las=2, col=colvec, border=NA, ylim=c(-maxy*1.5,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab=ylab)
  h_light = barplot(-y_light, las=2, col=colvec, border=NA, ylim=c(-maxy*1.5,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab=ylab, add = T)
  
  segments(y0 = maxy*1.5, y1 = maxy*1.5, x0 = 0.5, x1 = 192.5,  col = "black")
  segments(y0 = -maxy*1.5, y1 = -maxy*1.5, x0 = 0.5, x1 = 192.5,  col = "black")
  segments(y0 = 0, y1 = 0, x0 = 0.5, x1 = 192.5,  col = "black")
  abline(v = 0.5, col = "black")
  abline(v = 32.5, col = "black")
  abline(v = 64.5, col = "black")
  abline(v = 96.5, col = "black")
  abline(v = 128.5, col = "black")
  abline(v = 160.5, col = "black")
  abline(v = 192.5, col = "black")
  
  
  for (j in 1:length(sub_vec)) {
    xpos = h_heavy[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.25, xpos[2]+0.5, maxy*1.15, border=NA, col=colvec[j*16])
    text(x=mean(xpos), y=maxy*1.15, pos=3, labels=sub_vec[j])
  }    
  dev.copy(pdf,file_name,width=12,height=5)
  dev.off()
  dev.off()
}

lower_bound <- c(0,0.0025,0.005,0.01,0.02,0.05,0.1,0.25,0.5,0.01,0)
upper_bound <- c(0.0025,0.005,0.01,0.02,0.05,0.1,0.25,0.5,1,1,0.01)

for(i in 1:length(lower_bound)){
  mutations <- all_data[which(all_data$vaf >= lower_bound[i] & all_data$vaf < upper_bound[i]),]
  trinucleotide_plot(mutations,paste("all_tissues_coding_region_counts_bin",i,"_",lower_bound[i],"to",upper_bound[i],"_pct_cutoff_trinucleotide_plot.pdf",sep=""),"counts","coding")
  trinucleotide_plot(mutations,paste("all_tissues_d_loop_counts_bin",i,"_",lower_bound[i],"to",upper_bound[i],"_pct_cutoff_trinucleotide_plot.pdf",sep=""),"counts","d_loop")
  trinucleotide_plot(mutations,paste("all_tissues_all_mtDNA_counts_bin",i,"_",lower_bound[i],"to",upper_bound[i],"_pct_cutoff_trinucleotide_plot.pdf",sep=""),"counts","all_mtDNA")
  
  trinucleotide_plot(mutations,paste("all_tissues_coding_region_obs_exp_bin",i,"_",lower_bound[i],"to",upper_bound[i],"_pct_cutoff_trinucleotide_plot.pdf",sep=""),"obs_exp","coding")
  trinucleotide_plot(mutations,paste("all_tissues_d_loop_obs_exp_bin",i,"_",lower_bound[i],"to",upper_bound[i],"_pct_cutoff_trinucleotide_plot.pdf",sep=""),"obs_exp","d_loop")
  trinucleotide_plot(mutations,paste("all_tissues_all_mtDNA_obs_exp_bin",i,"_",lower_bound[i],"to",upper_bound[i],"_pct_cutoff_trinucleotide_plot.pdf",sep=""),"obs_exp","all_mtDNA")
}

for(i in 1:length(tissues)){
  mutations <- all_data[which(all_data$tissue == tissues[i] & all_data$vaf >= 0.01),]
  
  trinucleotide_plot(mutations,paste(tissues[i],"_coding_region_counts_0.01_pct_cutoff_trinucleotide_plot.pdf",sep=""),"counts","coding")
  trinucleotide_plot(mutations,paste(tissues[i],"_d_loop_counts_0.01_pct_cutoff_trinucleotide_plot.pdf",sep=""),"counts","d_loop")
  trinucleotide_plot(mutations,paste(tissues[i],"_all_mtDNA_counts_0.01_pct_cutoff_trinucleotide_plot.pdf",sep=""),"counts","all_mtDNA")
  
  trinucleotide_plot(mutations,paste(tissues[i],"_coding_region_obs_exp_0.01_pct_cutoff_trinucleotide_plot.pdf",sep=""),"obs_exp","coding")
  trinucleotide_plot(mutations,paste(tissues[i],"_d_loop_obs_exp_0.01_pct_cutoff_trinucleotide_plot.pdf",sep=""),"obs_exp","d_loop")
  trinucleotide_plot(mutations,paste(tissues[i],"_all_mtDNA_obs_exp_0.01_pct_cutoff_trinucleotide_plot.pdf",sep=""),"obs_exp","all_mtDNA")
}

for(i in 1:length(tissues)){
  mutations <- all_data[which(all_data$tissue == tissues[i] & all_data$vaf < 0.01),]
  
  trinucleotide_plot(mutations,paste(tissues[i],"_coding_region_counts_below_0.01_pct_cutoff_trinucleotide_plot.pdf",sep=""),"counts","coding")
  trinucleotide_plot(mutations,paste(tissues[i],"_d_loop_counts_below_0.01_pct_cutoff_trinucleotide_plot.pdf",sep=""),"counts","d_loop")
  trinucleotide_plot(mutations,paste(tissues[i],"_all_mtDNA_counts_below_0.01_pct_cutoff_trinucleotide_plot.pdf",sep=""),"counts","all_mtDNA")
  
  trinucleotide_plot(mutations,paste(tissues[i],"_coding_region_obs_exp_below_0.01_pct_cutoff_trinucleotide_plot.pdf",sep=""),"obs_exp","coding")
  trinucleotide_plot(mutations,paste(tissues[i],"_d_loop_obs_exp_below_0.01_pct_cutoff_trinucleotide_plot.pdf",sep=""),"obs_exp","d_loop")
  trinucleotide_plot(mutations,paste(tissues[i],"_all_mtDNA_obs_below_exp_0.01_pct_cutoff_trinucleotide_plot.pdf",sep=""),"obs_exp","all_mtDNA")
}

mutations <- all_data[which(all_data$tissue == "bladder_wgs" & all_data$donor != "PD38778" & all_data$vaf <= 0.01),]
trinucleotide_plot(mutations,"all_bladder_except_PD38778_all_mtDNA_obs_exp_below_0.01_pct_cutoff_trinucleotide_plot.pdf","obs_exp","all_mtDNA")
mutations <- all_data[which(all_data$tissue == "bladder_wgs" & all_data$donor != "PD38778" & all_data$vaf > 0.01),]
trinucleotide_plot(mutations,"all_bladder_except_PD38778_all_mtDNA_obs_exp_above_0.01_pct_cutoff_trinucleotide_plot.pdf","obs_exp","all_mtDNA")
mutations <- all_data[which(all_data$tissue == "bladder_wgs" & all_data$donor == "PD38778" & all_data$vaf <= 0.01),]
trinucleotide_plot(mutations,"bladder_only_PD38778_all_mtDNA_obs_exp_below_0.01_pct_cutoff_trinucleotide_plot.pdf","obs_exp","all_mtDNA")
mutations <- all_data[which(all_data$tissue == "bladder_wgs" & all_data$donor == "PD38778" & all_data$vaf > 0.01),]
trinucleotide_plot(mutations,"bladder_only_PD38778_all_mtDNA_obs_exp_above_0.01_pct_cutoff_trinucleotide_plot.pdf","obs_exp","all_mtDNA")

mutations <- all_data[which(all_data$vaf > 0.01),]
trinucleotide_plot(mutations,"all_tissues_all_mtDNA_obs_exp_above_0.01_pct_cutoff_trinucleotide_plot.pdf","obs_exp","all_mtDNA")
