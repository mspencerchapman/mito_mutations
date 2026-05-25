all_files <- list.files(".")
bam_files <- all_files[grepl(pattern = ".bam$", x = all_files)]
bai_files <- all_files[grepl(pattern = ".bam.bai$", x = all_files)]

for(i in 1:length(bam_files)){
  if(!(paste0(bam_files[i],".bai") %in% bai_files)){
    system(sprintf("samtools index -b %s", bam_files[i]))
  }
}