args <- commandArgs(TRUE)
tissue <- as.character(args[1])

setwd("/lustre/scratch125/casm/team268im/al28/mtDNA/")

wgs_coverage <- read.table(file = "whole_genome_coverage_pileup_and_bedtools_annotated.csv", sep = ",", header = T)

tissue_list <- as.character(unique(wgs_coverage$Tissue))

i<-which(tissue_list==tissue)

##RUN WITH THIS INDEX TO RUN OVER THIS TISSUE ONLY
tissue_samples <- as.character(wgs_coverage$Sample[which(wgs_coverage$Tissue == tissue_list[i])])

if(tissue_list[i] == "blood"){
  tissue_samples <- paste0("TG01_",tissue_samples)
  donors <- "TG01"
} else if(tissue_list[i] == "colon_ibd"){
  # metadata <- read.table("/Users/al28/Mounts/cancer/mtDNA/project_metadata/colon_ibd.txt", sep = "\t", stringsAsFactors = F, header = T)
  # donors <- length(unique(metadata$patient_ID))
  donors <- sort(unique(gsub("N.*","",unique(gsub("B.*","",tissue_samples, fixed = FALSE)),fixed = FALSE)))
} else if(tissue_list[i] == "colon"){
  colon_metadata <- read.table("/lustre/scratch125/casm/team268im/al28/mtDNA/project_metadata/colon/colon_sample_info.txt", sep = "\t", stringsAsFactors = F, header = T)
  donors <- unique(c(unique(colon_metadata$patient), unique(substr(tissue_samples[which(!(tissue_samples %in% colon_metadata$crypt) & substr(tissue_samples,1,3) != "HLS")],1,7))))
  donors <- donors[which(!(donors %in% c("PD36813","PD36814","PD37226","PD37457","PD37265")))]
} else if(tissue_list[i] == "immune"){
  PD_samples <- tissue_samples[which(substr(tissue_samples,1,2) == "PD")]
  donors <- c("TG01",unique(substr(PD_samples,1,7)))
} else{
  donors <- unique(substr(tissue_samples,1,7))
}

for(j in 1:length(donors)){
  tissue <- tissue_list[i]
  donor <- donors[j]
  job_id <- paste(tissue,donor,"shearwater",sep = "_")
  log_path <- paste0("/lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/shearwater_stringent_exclusion/",tissue,"/logs/",job_id)
  shearwater_fun <- "/lustre/scratch125/casm/team268im/al28/mtDNA/al28_code/mtDNA_shearwater_farm_script.R"
  
  system(sprintf("bsub -J %s -q normal -R 'select[mem>=5000] rusage[mem=5000]' -M5000 -e %s.err -o %s.out Rscript %s %s %s", job_id, log_path, log_path, shearwater_fun, tissue, donor))
  
}
    
    
    
    
    
    
    
    
      
    