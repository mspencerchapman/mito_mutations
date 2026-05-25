args <- commandArgs(TRUE)

sample_table <- as.character(args[1])

samples <- read.table(sample_table)

all_files <- list.files(".")
bam_files <- all_files[grepl(pattern = ".bam$", x = all_files)]

for (i in 1:nrow(samples)) {
  if(!(paste0(samples[i,2],"_MT.bam") %in% bam_files)){
    system(sprintf("bsub -q normal -e logs/%s.err -o logs/%s.out samtools view -h /nfs/cancer_ref01/nst_links/live/%s/%s/%s.sample.dupmarked.bam chrM -o %s_MT.bam -b",
                   samples[i,2], samples[i,2], samples[i,1], samples[i,2], samples[i,2], samples[i,2]))
  }
}
