library("GenomicRanges")
library("deepSNV")
library("Rsamtools")

estimateRho_gridml = function(x, mu, prior = NULL) {
  # By Tim Coorens
  rhovec = 10 ^ seq(-6, -0.5, by = 0.05) # rho will be bounded within 1e-6 and 0.32
  mm = c(x[, 3:4])
  cov = c(x[, 1:2])
  if (is.null(prior)) {
    ll = sapply(rhovec, function(rhoj)
      sum(dbetabinom(
        x = mm,
        size = cov,
        rho = rhoj,
        prob = mu,
        log = T
      )))
  }
  else{
    ll = sapply(rhovec, function(rhoj)
      sum(dbetabinom(
        x = mm,
        size = cov,
        rho = rhoj,
        prob = mu,
        log = T
      ))) + log(prior)
  }
  rhovec[ll == max(ll)][1]
}

setwd("/Users/al28/Mounts/cancer/mtDNA/")

wgs_coverage <- read.table(file = "whole_genome_coverage.csv", sep = ",", header = T)

genome_file = "/Users/al28/Mounts/cancer/bed_ref/hs37d5.fa"
prior        <- NULL

logbb = deepSNV:::logbb
dbetabinom = VGAM::dbetabinom

min_totest = 2

tissue_list <- as.character(unique(wgs_coverage$Tissue))

ref_list <- c("A","T","C","G","-","INS")

for(i in 1:length(tissue_list)){
  tissue_samples <- as.character(wgs_coverage$Sample[which(wgs_coverage$Tissue == tissue_list[i])])
  
  if(tissue_list[i] == "blood"){
    tissue_samples <- paste0("TG01_",tissue_samples)
    donors <- "TG01"
  } else if(tissue_list[i] == "colon_ibd"){
    # metadata <- read.table("/Users/al28/Mounts/cancer/mtDNA/project_metadata/colon_ibd.txt", sep = "\t", stringsAsFactors = F, header = T)
    # donors <- length(unique(metadata$patient_ID))
    donors <- sort(unique(gsub("N.*","",unique(gsub("B.*","",tissue_samples, fixed = FALSE)),fixed = FALSE)))
  } else if(tissue_list[i] == "colon"){
    colon_metadata <- read.table("/Users/al28/Mounts/cancer/mtDNA/project_metadata/colon/colon_sample_info.txt", sep = "\t", stringsAsFactors = F, header = T)
    donors <- unique(c(unique(colon_metadata$patient), unique(substr(tissue_samples[which(!(tissue_samples %in% colon_metadata$crypt) & substr(tissue_samples,1,3) != "HLS")],1,7))))
  } else{
    donors <- unique(substr(tissue_samples,1,7))
  }
  
  for(j in 1:length(donors)){
    if(tissue_list[i] == "blood"){
      donor_samples <- as.character(wgs_coverage$Sample[which(wgs_coverage$Tissue == tissue_list[i])])
    } else if(tissue_list[i] == "colon_ibd"){
      donor_samples <- tissue_samples[which(substr(tissue_samples,1,nchar(donors[j]) + 1) %in% c(paste0(donors[j],"B"),paste0(donors[j],"N")))]
    } else if(tissue_list[i] == "colon"){
      if(donors[j] == "HLS"){
        donor_samples <- tissue_samples[which(substr(tissue_samples,1,3) == "HLS")]
      } else{
        donor_samples <- unique(c(colon_metadata$crypt[which(colon_metadata$patient == donors[j])],tissue_samples[which(substr(tissue_samples,1,7) == donors[j])]))
      }
    }else{
      donor_samples <- tissue_samples[which(substr(tissue_samples,1,7) == donors[j])]
    }
    
    tumcounts_obj_all <- list()
    tumcounts_obj_all$counts <- array(data = NA, dim = c(length(donor_samples),16569,12))
    
    for(k in 1:length(donor_samples)){
      if(file.exists(paste0("Data/",tissue_list[i],"/",donor_samples[k],"_MT_count.csv"))){
        sample_cov <- read.table(paste0("Data/",tissue_list[i],"/",donor_samples[k],"_MT_count.csv"), sep = ",", header = T)
        for(n in 1:12){
          tumcounts_obj_all$counts[k,,n] <- sample_cov[,n]
        }
      }
    }
    tum_total <- apply(tumcounts_obj_all$counts[, , 1:6] + tumcounts_obj_all$counts[, , 7:12], c(2, 3), sum)
    tum_total <- tum_total / rowSums(tum_total)
    
    normcounts_obj <- list()
    normcounts_obj$counts <- readRDS(file = "shearwater_background.Rds")
    norm_total <- apply(normcounts_obj$counts[, , 1:6] + normcounts_obj$counts[, , 7:12], c(2, 3), sum)
    norm_total <- norm_total / rowSums(norm_total)
    
    mutations_allsamples = NULL
    for (k in 1:length(donor_samples)) {
      tumcounts_obj = tumcounts_obj_all
      tumcounts_obj$counts = array(tumcounts_obj$counts[k, , ], dim = c(1, dim(tumcounts_obj$counts[k, , ])))
      
      mutsites = which((tumcounts_obj$counts[1, , 1:6] + tumcounts_obj$counts[1, , 7:12]) >= min_totest & norm_total < 0.4, arr.ind = T)
      
      if (nrow(mutsites) > 0) {
        mutations = data.frame(
          sampleID = donor_samples[k],
          chr = "MT",
          pos = c(1:16569)[mutsites[, 1]],
          ref = NA,
          mut = ref_list[mutsites[, 2]],
          xfw = NA,
          xbw = NA,
          nfw = NA,
          nbw = NA,
          mu = NA,
          pval = NA,
          rho = NA
        )
        mutations$tum_globalvaf = apply(mutsites, 1, function(x) tum_total[x[1], x[2]])
        
        l = 6
        indsm = cbind(mutsites[, 2], mutsites[, 2] + l)
        sites_gr = GRanges(mutations$chr, IRanges(mutations$pos, mutations$pos))
        seqs = as.character(scanFa(genome_file, sites_gr))
        mutations$ref = as.character(seqs)
        
        for(x in 1:nrow(mutations)){
          inds = indsm[x, ]
          norcounts = normcounts_obj$counts[, mutsites[x], ]
          tumcounts = tumcounts_obj$counts[, mutsites[x], ]
          
          pseudo = .Machine$double.eps
          x.fw = tumcounts[inds[1]]
          x.bw = tumcounts[inds[2]]
          n.fw = sum(tumcounts[1:l])
          n.bw = sum(tumcounts[(l + 1):(2 * l)])
          
          X.fw = sum(norcounts[, inds[1]])
          X.bw = sum(norcounts[, inds[2]])
          N.fw = sum(norcounts[, 1:l])
          N.bw = sum(norcounts[, (l + 1):(2 * l)])
          mu = max(X.fw + X.bw, pseudo) / max(N.fw + N.bw, pseudo)
          
          counts = cbind(rowSums(norcounts[, 1:l]), rowSums(norcounts[, (l +
                                                                           1):(2 * l)]), norcounts[, inds])
          rho = estimateRho_gridml(counts, mu, prior)
          
          rdisp = (1 - rho) / rho
          
          prob0.fw = (X.fw + x.fw) / (N.fw + n.fw)
          prob0.fw[prob0.fw == 0] = pseudo
          prob1s.fw = x.fw / (n.fw + pseudo)
          prob1s.fw[prob1s.fw == 0] = pseudo
          prob1c.fw = X.fw / (N.fw + pseudo)
          prob1c.fw[prob1c.fw == 0] = pseudo
          prob1s.fw = pmax(prob1s.fw, prob1c.fw) # Min error rate is that of the population (one-sided test)
          # Modified by fa8 to avoid p=1, which results in beta=0
          nu0.fw = prob0.fw * rdisp
          nu1s.fw = min(1 - pseudo, prob1s.fw) * rdisp
          nu1c.fw = min(1 - pseudo, prob1c.fw) * rdisp
          
          prob0.bw = (X.bw + x.bw) / (N.bw + n.bw)
          prob0.bw[prob0.bw == 0] = pseudo
          prob1s.bw = x.bw / (n.bw + pseudo)
          prob1s.bw[prob1s.bw == 0] = pseudo
          prob1c.bw = X.bw / (N.bw + pseudo)
          prob1c.bw[prob1c.bw == 0] = pseudo
          prob1s.bw = pmax(prob1s.bw, prob1c.bw) # Min error rate is that of the population (one-sided test)
          # Modified by fa8 to avoid p=1, which results in beta=0
          nu0.bw = prob0.bw * rdisp
          nu1s.bw = min(1 - pseudo, prob1s.bw) * rdisp
          nu1c.bw = min(1 - pseudo, prob1c.bw) * rdisp
          
          LL.fw = logbb(x.fw, n.fw, nu0.fw, rdisp) + logbb(X.fw, N.fw, nu0.fw, rdisp) - logbb(x.fw, n.fw, nu1s.fw, rdisp) - logbb(X.fw, N.fw, nu1c.fw, rdisp)
          pvals_fw = pchisq(-2 * LL.fw, df = 1, lower.tail = F) / 2 # We divide by 2 as we are performing a 1-sided test
          LL.bw = logbb(x.bw, n.bw, nu0.bw, rdisp) + logbb(X.bw, N.bw, nu0.bw, rdisp) - logbb(x.bw, n.bw, nu1s.bw, rdisp) - logbb(X.bw, N.bw, nu1c.bw, rdisp)
          pvals_bw = pchisq(-2 * LL.bw, df = 1, lower.tail = F) / 2 # We divide by 2 as we are performing a 1-sided test
          pvals_both = pchisq(-2 * (log(pvals_fw) + log(pvals_bw)), 4, low =
                                F) # Fisher's combined p-value
          mutations[x, 6:12] = c(x.fw, x.bw, n.fw, n.bw, mu, pvals_both, rho)
        }
        mutations_allsamples = rbind(mutations_allsamples, mutations)
      }
    }

    mutations_allsamples = mutations_allsamples[which(mutations_allsamples$pval < 1e-3), ]
    
    write.table(x = mutations_allsamples, file = paste0("Analysis/shearwater/",tissue_list[i],"/",donors[j],"_shearwater_mutations.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
  }
}