files <- list.files(".")
vcf_files <- files[grepl(pattern = "vcf.gz", x = files)]

for (i in 1:length(vcf_files)) {
  
  sample <- substr(vcf_files[i],1,nchar(vcf_files[i]) - 7)
  
  if(!(paste0(sample,".haplocheck") %in% files)){
    system(sprintf("bsub -q normal -o logs/%s.haplocheck.out -e logs/%s.haplocheck.err -n1 -M500 -R'select[mem>500] rusage[mem=500] span[hosts=1]' /lustre/scratch125/casm/team268im/al28/software/haplocheck --out %s.haplocheck %s.vcf.gz",
                   sample, sample, sample, sample))
  }
}