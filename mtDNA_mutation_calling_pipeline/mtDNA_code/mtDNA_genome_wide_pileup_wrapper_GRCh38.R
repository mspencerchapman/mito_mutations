files <- list.files()
bam_files <- vector()

for(i in 1:length(files)){
  if(substr(files[i],nchar(files[i]) - 3,nchar(files[i])) == ".bam"){
    bam_files <- c(bam_files,files[i])
  }
}

for(i in 1:length(bam_files)){
  if(!(paste0(substr(bam_files[i],1,nchar(bam_files[i]) - 4),"_count.csv") %in% files)){
    system(sprintf("bsub -q long -M3000 -R'select[mem>3000] rusage[mem=3000] span[hosts=1]' -n1 -e logs/%s.err -o logs/%s.out Rscript /lustre/scratch126/casm/team154pc/ms56/my_programs/mtDNA_code/mtDNA_genome_wide_pileup_GRCh38.R %s",
                   bam_files[i], bam_files[i], bam_files[i]))
  }
}
