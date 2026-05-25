library(hdp)
library(tidyverse)

setwd("/Users/al28/Mounts/cancer/mtDNA/Analysis/hdp/pcawg_shearwater_tumour_only/")

analysis_list <- list.files(".")

for(x in 1:length(analysis_list)){
  setwd(analysis_list[x])
  
  files <- list.files(".")
  input_file <- files[grepl("HDPinput",files)]
  
  mut_count <- read.table(input_file, sep = "\t", stringsAsFactors = F)
  mut_count <- mut_count[which(as.vector(rowSums(mut_count)) != 0),]
  
  tissues <- table(sub("_[^_]+$", "", rownames(mut_count))) 
  
  mut_category <- as.data.frame(array(data = NA, dim = c(length(tissues) * 12,4)))
  colnames(mut_category) <- c("tissue","mut_type","count","fraction")
  mut_category$tissue <- rep(names(tissues), each = 12)
  
  mut_category$mut_type <- rep(c("C>A_L","C>G_L","C>T_L","T>A_L","T>C_L","T>G_L","C>A_H","C>G_H","C>T_H","T>A_H","T>C_H","T>G_H") , times = length(tissues))
  
  cum_count <- 0
  for(i in 1:length(tissues)){
    temp_count <- mut_count[((cum_count + 1):(cum_count + (tissues[i]))),]
    for(j in 1:12){
      mut_category$count[((12 *(i - 1)) + j)] <- sum(temp_count[(1:nrow(temp_count)),((16 * (j - 1)) + 1):((16 * (j - 1)) + 16)])
    }
    for(j in 1:12){
      mut_category$fraction[((12 *(i - 1)) + j)] <- mut_category$count[((12 *(i - 1)) + j)] / sum(temp_count)
    }
    cum_count <- cum_count + as.numeric(tissues[i])
  }
  
  mut_category$reverse_alphabetical <- rep(length(tissues):1, each = 12)
  mut_category$prop_heavy_T_C <- NA
  
  order_list <- order(mut_category$fraction[which(mut_category$mut_type == "T>C_H")])
  
  for(i in 1:length(tissues)){
    mut_category$prop_heavy_T_C[((i - 1) * 12) + (1:12)] <- which(order_list == i)
  }
  
  ggplot(data = mut_category) +
    geom_bar(aes(y = reorder(tissue, reverse_alphabetical), x = fraction, fill = mut_type), position = "stack", stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = c("dodgerblue","#c0eff9","black","#bfbfbf","red","#fac9cc","grey70","#f5f1f2","olivedrab3","#e7f4da","plum2","#fcf0f0")) +
    theme(axis.text.y = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = "Proportion of mutations", y = "Tissue", fill = "Mutation Type")
  ggsave("stacked_mut_type_bar_plot_alphabetical.pdf", width = 10, height = 10)
  
  ggplot(data = mut_category) +
    geom_bar(aes(y = reorder(tissue, prop_heavy_T_C), x = fraction, fill = mut_type), position = "stack", stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = c("dodgerblue","#c0eff9","black","#bfbfbf","red","#fac9cc","grey70","#f5f1f2","olivedrab3","#e7f4da","plum2","#fcf0f0")) +
    theme(axis.text.y = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = "Proportion of mutations", y = "Tissue", fill = "Mutation Type")
  ggsave("stacked_mut_type_bar_plot_prop_heavy_T_C.pdf", width = 10, height = 10)
  
  setwd("no_hierarchy")
  
  chlist <- vector("list", 10)
  
  for(i in 1:10) {
    chlist[[i]] <- readRDS(paste0("chlist_",i,".rds"))
  }
  
  mut_example_multi <- hdp_multi_chain(chlist)
  
  pdf(file = "likelihood_plot_no_hierarchy.pdf", width = 10, height = 20, useDingbats = F)
  par(mfrow=c(5,2), mar=c(4, 4, 2, 1))
  p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start = 1000)
  dev.off()
  
  pdf(file = "cluster_num_plot_no_hierarchy.pdf", useDingbats = F)
  par(mfrow=c(5,2), mar=c(4, 4, 2, 1))
  p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
  dev.off()
  
  pdf(file = "prop_data_assigned_plot_no_hierarchy.pdf", useDingbats = F)
  par(mfrow=c(5,2), mar=c(4, 4, 2, 1))
  p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")
  dev.off()
  
  mut_example_multi <- hdp_extract_components(mut_example_multi)
  
  pdf(file = "prop_data_items_assigned_plot_no_hierarchy.pdf", width = 12, height = 12, useDingbats = F)
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
  plot_comp_size(mut_example_multi, bty="L")
  dev.off()
  
  trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
  group_factor <- as.factor(c(paste("light",rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),each=16), sep = "_"),paste("heavy",rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),each=16), sep = "_")))
  mut_colours <- rep(c(RColorBrewer::brewer.pal(10, 'Paired')[seq(1,10,2)], 'grey70'),times = 2)
  
  pdf(file = "sigs_extracted_signatures_plot_no_hierarchy.pdf", height = 5, width = 15, useDingbats = F)
  plot_comp_distn(mut_example_multi, cat_names=trinuc_context,
                  grouping=group_factor, col=mut_colours,
                  col_nonsig="grey80", show_group_labels=TRUE)
  dev.off()
  
  dp_distn <- comp_dp_distn(mut_example_multi)
  
  pdf(file = "signatures_per_sample_plot_no_hierarchy_nonsig_not_shown.pdf", useDingbats = F)

    dpindices <- 1+(1:length(tissues))
   
    plot_dp_comp_exposure(mut_example_multi,
                            dpindices=dpindices, 
                            col=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
                            incl_nonsig=F,
                            dpnames = names(tissues),
                            ylab_numdata = 'SNV count', ylab_exp = 'Signature exposure',
                            leg.title = 'Signature', las = 2, mar = c(4, 4, 2, 0.5))
  dev.off()
  
  pdf(file = "signatures_per_sample_plot_no_hierarchy_nonsig_shown.pdf", useDingbats = F)
  plot_dp_comp_exposure(mut_example_multi,
                        dpindices=dpindices, 
                        col=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
                        incl_nonsig=T,
                        dpnames = names(tissues),
                        ylab_numdata = 'SNV count', ylab_exp = 'Signature exposure',
                        leg.title = 'Signature', las = 2, mar = c(4, 4, 2, 0.5))
  dev.off()
  
  write.table(comp_categ_distn(mut_example_multi)$mean, file = "trinuc_means.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
  
  mean_sig_assignment <- comp_dp_distn(mut_example_multi)$mean
  mean_sig_assignment <- mean_sig_assignment[-1,]
  rownames(mean_sig_assignment) <- rownames(mut_count)
  
  saveRDS(mean_sig_assignment, file = "mean_assignment_per_sample.rds")
  
  setwd("../..")
}
