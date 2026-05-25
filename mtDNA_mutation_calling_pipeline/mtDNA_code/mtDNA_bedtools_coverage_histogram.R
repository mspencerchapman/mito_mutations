args <- commandArgs(TRUE)

sample_table <- as.character(args[1])
sample_type <- substr(sample_table,1,nchar(sample_table) - 4)

samples <- read.table(sample_table)

files <- list.files(paste0("/lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/bedtools_coverage/",sample_type))

for (i in 1:nrow(samples)) {
  if(!(paste0(samples[i,2],".genomecov.txt") %in% files) & !(paste0(samples[i,2],".genomecov.txt.final") %in% files)){
  system(sprintf("bsub -q normal -M3000 -R'select[mem>3000] rusage[mem=3000] span[hosts=1]' -R'select[casm_cgpirods>=20] rusage[casm_cgpirods=20]' -e /lustre/scratch125/casm/team268im/al28/mtDNA/Data/%s/logs/%s.err -o /lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/bedtools_coverage/%s/%s.genomecov.txt bedtools genomecov -ibam /nfs/cancer_ref01/nst_links/live/%s/%s/%s.sample.dupmarked.bam",
                 sample_type, samples[i,2], sample_type, samples[i,2], samples[i,1], samples[i,2], samples[i,2]))
  }
}
