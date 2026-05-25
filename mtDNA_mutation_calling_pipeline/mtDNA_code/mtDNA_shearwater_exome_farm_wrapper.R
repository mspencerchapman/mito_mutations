setwd("/lustre/scratch125/casm/team268im/al28/mtDNA/")

exome_coverage <- read.table(file = "exome_coverage_pileup_and_bedtools_annotated.csv", sep = ",", header = T)

tissue_list <- as.character(unique(exome_coverage$Tissue))

for(i in 1:length(tissue_list)){
  tissue_samples <- as.character(exome_coverage$Sample[which(exome_coverage$Tissue == tissue_list[i])])
  
  donors <- unique(substr(tissue_samples,1,7))
  
  for(j in 1:length(donors)){
    tissue <- tissue_list[i]
    donor <- donors[j]
    job_id <- paste(tissue,donor,"shearwater",sep = "_")
    log_path <- paste0("/lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/shearwater_exome/",tissue,"/logs/",job_id)
    shearwater_fun <- "/lustre/scratch125/casm/team268im/al28/mtDNA/al28_code/mtDNA_shearwater_exome_farm_script.R"
    
    system(sprintf("bsub -J %s -q yesterday -R 'select[mem>=2000] rusage[mem=2000]' -M2000 -e %s.err -o %s.out /software/R-3.6.1/bin/Rscript %s %s %s", job_id, log_path, log_path, shearwater_fun, tissue, donor))
    
  }
}
    
    
    
    
    
    
    
    
      
    