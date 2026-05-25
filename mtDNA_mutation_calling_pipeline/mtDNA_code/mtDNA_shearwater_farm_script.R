library("GenomicRanges")
library("deepSNV")
library("Rsamtools")

logbb = deepSNV:::logbb
dbetabinom = VGAM::dbetabinom

args <- commandArgs(TRUE)
tissue <- args[1]
donor <- args[2]

setwd("/lustre/scratch125/casm/team268im/al28/mtDNA/")

wgs_coverage <- read.table(file = "whole_genome_coverage_pileup_and_bedtools_annotated.csv", sep = ",", header = T, stringsAsFactors = F)
excluded_samples <- read.table(file = "combined_sample_exclusion.tsv", sep = "\t", header = T, stringsAsFactors = F)

genome_file = "/lustre/scratch125/casm/team268im/al28/bed_ref/hs37d5.fa"
prior        <- NULL

min_totest = 2

ref_list <- c("A","T","C","G","-","INS")

normcounts_obj <- list()
normcounts_obj$counts <- readRDS(file = "shearwater_background.Rds")
norm_total <- apply(normcounts_obj$counts[, , 1:6] + normcounts_obj$counts[, , 7:12], c(2, 3), sum)
norm_total <- norm_total / rowSums(norm_total)

normal_rho_panel <- read.table("/lustre/scratch125/casm/team268im/al28/mtDNA/shearwater_normal_panel_rho_mu_calculation.tsv", sep = "\t", stringsAsFactors = F, header = T)
normal_rho_panel$rho[which(normal_rho_panel$rho <= 2e-04)] <- 2e-04

tissue_samples <- as.character(wgs_coverage$Sample[which(wgs_coverage$Tissue == tissue)])
tissue_samples <- tissue_samples[which(!(tissue_samples %in% excluded_samples$Sample[which(excluded_samples$Tissue == tissue & excluded_samples$CombinedExclusion == "Y")]))]

if(tissue == "blood"){
  donor_samples <- as.character(wgs_coverage$Sample[which(wgs_coverage$Tissue == tissue)])
} else if(tissue == "colon_ibd"){
  donor_samples <- tissue_samples[which(substr(tissue_samples,1,nchar(donor) + 1) %in% c(paste0(donor,"B"),paste0(donor,"N")))]
} else if(tissue == "colon"){
  if(donor == "HLS"){
    donor_samples <- tissue_samples[which(substr(tissue_samples,1,3) == "HLS")]
  } else{
    colon_metadata <- read.table("/lustre/scratch125/casm/team268im/al28/mtDNA/project_metadata/colon/colon_sample_info.txt", sep = "\t", stringsAsFactors = F, header = T)
    donor_samples <- unique(c(colon_metadata$crypt[which(colon_metadata$patient == donor)],tissue_samples[which(substr(tissue_samples,1,7) == donor)]))
  }
} else if(tissue == "immune"){
  if(donor == "TG01"){
    donor_samples <- tissue_samples[which(substr(tissue_samples,1,2) != "PD")]
  } else{
    donor_samples <- tissue_samples[which(substr(tissue_samples,1,7) == donor)]
  }
} else{
  donor_samples <- tissue_samples[which(substr(tissue_samples,1,7) == donor)]
}

tumcounts_obj_all <- list()
tumcounts_obj_all$counts <- array(data = NA, dim = c(length(donor_samples),16569,12))

for(k in 1:length(donor_samples)){
  if(file.exists(paste0("Data/",tissue,"/",donor_samples[k],"_MT_count.csv"))){
    sample_cov <- read.table(paste0("Data/",tissue,"/",donor_samples[k],"_MT_count.csv"), sep = ",", header = T)
    for(n in 1:12){
      tumcounts_obj_all$counts[k,,n] <- sample_cov[,n]
    }
  }
}

if (length(donor_samples) > 1) {
  tum_total <- apply(tumcounts_obj_all$counts[, , 1:6] + tumcounts_obj_all$counts[, , 7:12], c(2, 3), sum)
} else{
  tum_total = tumcounts_obj_all$counts[1, , 1:6] + tumcounts_obj_all$counts[1, , 7:12]
}
tum_total <- tum_total / rowSums(tum_total)

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
      rho = normal_rho_panel$rho[which(normal_rho_panel$pos == mutsites[x,1] & normal_rho_panel$mut == ref_list[mutsites[x,2]])]
      
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

write.table(x = mutations_allsamples, file = paste0("Analysis/shearwater_stringent_exclusion/",tissue,"/",donor,"_shearwater_mutations.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
