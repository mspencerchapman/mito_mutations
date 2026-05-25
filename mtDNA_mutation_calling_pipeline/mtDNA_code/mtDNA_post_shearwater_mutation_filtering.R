library("GenomicRanges")
library("rtracklayer")
library("deepSNV")
library("ggplot2")

setwd("/Users/al28/Mounts/cancer/mtDNA/Analysis/shearwater_stringent_exclusion/")

tissues <- list.files()

for(tiss in 1:length(tissues)){
  setwd(paste0("/Users/al28/Mounts/cancer/mtDNA/Analysis/shearwater_stringent_exclusion/",tissues[tiss]))
  all_files <- list.files(".")
  input_files <- all_files[grepl(pattern = "_shearwater_mutations.tsv", all_files)]
  
  for(x in 1:length(input_files)){
    dataset_name <- gsub(pattern = "_shearwater_mutations.tsv", replacement = "", x = input_files[x])
    mutations = read.table(input_files[x], sep = "\t", stringsAsFactors = F, header = T)
    
    mutations$ref[which(mutations$ref == TRUE)] <- "T"
    mutations$mut[which(mutations$mut == TRUE)] <- "T"
    
    indels_f <- length(grep("[-I]",mutations[,"mut"]));
    subs_f   <- nrow(mutations)-indels_f;
    cat("#INITIAL_NUMBER_OF_MUTATIONS\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");
    
    mutations$disc <- substr(mutations$sampleID,8,8);
    mutations$mut_site <- paste(mutations$chr,mutations$pos,mutations$ref,mutations$mut);
    mutations$vaf = (mutations$xfw+mutations$xbw)/(mutations$nfw+mutations$nbw)
    
    L = 16569
    s = unique(mutations$sampleID)
    
    indels         <- mutations[which(mutations$mut=="-"),];
    indels         <- indels[order(indels$sampleID, indels$chr, indels$pos),];
    i              <- 1;
    group_counter  <- 1;
    while(i <= nrow(indels)) {
      indel_from  <- indels[i,"pos"];
      indel_to    <- indel_from;
      from_index  <- i;
      to_index    <- i;
      deleted_seq <- indels[i,"ref"];
      sample_f    <- indels[i,"sampleID"];
      i <- i+1;
      for(j in c(i:(nrow(indels) + 1))) {
        if(j>nrow(indels)) {
          i <- j;
          break; #finish
        }
        else if(indels[j,"sampleID"] != indels[(j-1),"sampleID"]){
          i <- j;
          break;
        } else {
          if(indels[j,"pos"] != (indels[(j-1),"pos"]+1)) {
            i <- j; #start a new indel
            break;
          } else {
            if(indels[j,"chr"] != indels[(j-1),"chr"]) {
              i <- j; #start a new indel
              break;
            } else {
              #This is a candidate, but check their VAFs are compatible with a Fishers exact test:
              mat <- matrix(nrow=2,ncol=2,0);
              #mat[1,] <- c(indels[max(i,j-1),  "xfw"]+indels[max(i,j-1),  "xbw"], indels[max(i,j-1),  "nfw"]+indels[max(i,j-1),  "nbw"])
              mat[1,] <- c(indels[j-1,  "xfw"]+indels[j-1,  "xbw"], indels[j-1,  "nfw"]+indels[j-1,  "nbw"])
              mat[2,] <- c(indels[j,    "xfw"]+indels[j,    "xbw"], indels[j,    "nfw"]+indels[j,    "nbw"])
              mat[1,2] <- mat[1,2]-mat[1,1];
              mat[2,2] <- mat[2,2]-mat[2,1];
              pvalue <- fisher.test(mat)$p.value;
              if(pvalue < 0.01) {
                cat(" Breaking up indel because VAFs do not match\n");
                cat("             ",j-1, " vs ", j, ": pval=", pvalue, " [",indels[j-1,"pos"],"-",indels[j,"pos"],"]",sep="");
                cat("    (mat=", mat[1,1],",",mat[1,2],",",mat[2,1],",",mat[2,2],")\n",sep="");
                i <- j;
                break;
              }		
              indel_to <- indels[j,"pos"];
              to_index <- j;
              deleted_seq <- paste(deleted_seq,indels[j,"ref"],sep="");
            }
          }
        }
      }
      #cat("   [Sample=",sample_f,"] Indel goes from=",indel_from,", to=", indel_to," [",deleted_seq,">-]\n",sep="");
      indels[c(from_index:to_index),"groupID"    ] <- group_counter;
      indels[c(from_index:to_index),"deleted_seq"] <- deleted_seq;
      group_counter <- group_counter + 1;
    }
    mutations$indel_group <- NA;
    mutations$deleted_seq <- NA;
    
    if(nrow(indels) > 0){
      mutations[rownames(indels),c("indel_group","deleted_seq")] = indels[rownames(indels),c("groupID","deleted_seq")]
      
      putative_indelsites = mutations[mutations$mut %in% c("-","INS"),]
      putative_indelsites$qval = p.adjust(putative_indelsites$pval, method="BH", n=L*length(s)*2)
      putative_indelsites = unique(putative_indelsites[putative_indelsites$qval<0.20, c("sampleID","chr","pos")])
      indel_flank = 10
      putative_indelsites_gr = GRanges(putative_indelsites$chr, IRanges(putative_indelsites$pos-indel_flank, putative_indelsites$pos+indel_flank))
    }
    
    mutations <- unique(mutations);
    mutations$label <- "";
    
    ggplot(data = mutations, aes(x = vaf, y = tum_globalvaf)) + 
      geom_jitter(alpha = 0.5, size = 0.5) +
      labs(title = dataset_name, x = "VAF", y = "Global sample VAF") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="none")
    ggsave(paste0(dataset_name,"_global_vaf_vs_vaf.pdf"),width = 8, height = 8)
    
    germline_muts <- mutations[which(mutations$tum_globalvaf >= 0.9),]
    mutations[which(mutations$mut_site %in% germline_muts$mut_site),"label"] <- paste(mutations[which(mutations$mut_site %in% germline_muts$mut_site),"label"],"germline;",sep="");
    
    # #For the FDR to work properly I should remove these variants before:
    mutations <- mutations[-which(mutations$label == "germline;"),]
    
    indels_f <- length(grep("[-I]",mutations[which(mutations$label == ""),"mut"]));
    subs_f   <- nrow(mutations[which(mutations$label == ""),])-indels_f;
    cat("#AFTER_GERMLINE\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");
    
    if(nrow(mutations) > 0){
      mutations$qval = p.adjust(mutations$pval, method="BH", n=L*length(s)*5)
      mutations[which(mutations$qval>=0.01),"label"] = "no-fdr;";
      
      mutations = mutations[order(mutations$chr,mutations$pos),]
      
      indels_f <- length(grep("[-I]",mutations[which(mutations$label == ""),"mut"]));
      subs_f   <- nrow(mutations[which(mutations$label == ""),])-indels_f;
      cat("#AFTER_FDR\t",nrow(mutations[which(mutations$label == ""),]),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");
      
      rmpos = (mutations$xfw==0 | mutations$xbw==0) # Asking for 1 supporting read in both strands for all mutations
      filt2 = mutations[rmpos,]; 
      if(nrow(filt2)>0) {
        filt2$filter = "Strandness"
      }
      mutations[rmpos,"label"] = paste(mutations[rmpos,"label"],"strandness;",sep="");
      #####mutations = mutations[!rmpos,]
      samples = unique(mutations$sampleID)
      
      if(nrow(indels) > 0){
        rmpos = NULL
        for (h in 1:length(samples)) {
          m = mutations[mutations$sampleID==samples[h] & !(mutations$mut %in% c("-","INS")),]
          m_gr = GRanges(m$chr, IRanges(m$pos,m$pos))
          i_gr = putative_indelsites_gr[putative_indelsites$sampleID==samples[h]]
          ol = as.matrix(suppressWarnings(findOverlaps(m_gr, i_gr, type="any", select="all")))
          rmpos = c(rmpos, rownames(m)[unique(ol[,1])])
        }
        filt3 = mutations[rmpos,]; 
        if(nrow(filt3) > 0) {
          filt3$filter = "Near_indel"
        }
        mutations[rmpos,"label"] = paste(mutations[rmpos,"label"],"near_indel;",sep=""); 
      }
      
      indels_f <- length(grep("[-I]",mutations[which(mutations$label == ""),"mut"]));
      subs_f   <- nrow(mutations[which(mutations$label == ""),])-indels_f;
      cat("#AFTER_BOTH_STRANDS_FILT\t",nrow(mutations[which(mutations$label == ""),]),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");
      
      mutations[which(mutations$label == ""),"label"] <- "OK;";
      if(length(which(mutations$label == "OK;")) > 0){
        
        ok_muts <- mutations[grep("OK",mutations$label),      ]
        sub  <- ok_muts[grep("[-I]",ok_muts[,"mut"],invert=T),]
        ins  <- ok_muts[grep("I",   ok_muts[,"mut"]),         ]
        del  <- ok_muts[grep("-",   ok_muts[,"mut"]),         ]
        sub  <- sub[order(sub$sampleID, sub$chr, sub$pos),    ]
        ins  <- ins[order(ins$sampleID, ins$chr, ins$pos),    ]
        del  <- del[order(del$sampleID, del$chr, del$pos),    ]
        
        #To store the new data
        new_mutations <- mutations[0,]
        
        # Deletions (defined in mutations$indel_group):
        # For every "OK" deletion, get its del-groupID and find all the other deletions belonging
        # to that group. Merge them and create a new entry in mutations: combine pvalues, vaf, etc
        # For every "OK" deletion first check it hasn't been already merged
        for(j in 1:nrow(del)) {
          indel_group       <- del[j,"indel_group"]
          if(nrow(new_mutations[which(new_mutations$indel_group==indel_group),]) > 0) {
            next; #we already have one from the group of indels
          }
          indels_from_group                                  <- mutations[which(mutations$indel_group==indel_group),]
          new_mutations                                      <- rbind(new_mutations,indels_from_group[1,])
          new_mutations[nrow(new_mutations),"pos"          ] <- min  (indels_from_group$pos              )
          new_mutations[nrow(new_mutations),"vaf"          ] <- mean (indels_from_group$vaf              )
          new_mutations[nrow(new_mutations),"tum_globalvaf"] <- mean (indels_from_group$tum_globalvaf    )
          new_mutations[nrow(new_mutations),"pval"         ] <- min  (indels_from_group$pval             )
          new_mutations[nrow(new_mutations),"qval"         ] <- min  (indels_from_group$qval             )
          new_mutations[nrow(new_mutations),"label"        ] <- paste(indels_from_group$label,collapse="")
          new_mutations[nrow(new_mutations),"ref"          ] <- indels_from_group[1,"deleted_seq"]
        }
        
        # Insertions. No need to look for consecutive INS. Just add them to new_mutations
        new_mutations <- rbind(new_mutations,ins);
        
        
        # Substitutions: merge consecutive... [using Iñigo's code]
        d = sub$pos-(1:nrow(sub))
        runs = rle(d)
        rmpos = rep(0,nrow(sub))
        runstarts = cumsum(runs$length)-runs$length+1
        for (h in 1:length(runs$length)) {
          if (runs$length[h]>1) { # Adjacent mutations
            mutcluster                         = runstarts[h]:(runstarts[h]+runs$lengths[h]-1)
            rmpos[mutcluster[-1]             ] = 1 # Removing all the affected rows except the first one (which we will edit to capture the complex event)
            sub[mutcluster[1],"ref"          ] = paste(sub[mutcluster,"ref"          ],collapse="")
            sub[mutcluster[1],"mut"          ] = paste(sub[mutcluster,"mut"          ],collapse="")
            sub[mutcluster[1],"mu"           ] = mean (sub[mutcluster,"mu"           ]            )
            sub[mutcluster[1],"tum_globalvaf"] = mean (sub[mutcluster,"tum_globalvaf"]            )
            sub[mutcluster[1],"vaf"          ] = mean (sub[mutcluster,"vaf"          ]            )
            sub[mutcluster[1],"pval"         ] = min  (sub[mutcluster,"pval"         ]            )
            sub[mutcluster[1],"qval"         ] = min  (sub[mutcluster,"qval"         ]            )
          }
        }
        sub = sub[!rmpos,]        
        new_mutations <- rbind(new_mutations,sub);
        
        mutations.old <- mutations
        mutations     <- new_mutations
        mutations[which(mutations$label == ""),"label"] <- "OK;";
        mutations <- rbind(mutations, mutations.old[which(mutations.old$label == "near_indel;"),])
        
        indels_f <- length(grep("[-I]",mutations[,"mut"]));
        subs_f   <- nrow(mutations)-indels_f;
        cat("#AFTER_MERGING_RUNS\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");
        
        write.table(mutations, file=sprintf("%s_shearwater_calls.txt", dataset_name), col.names=T, row.names=F, sep="\t", quote=F)
      }
    }
  }
}
