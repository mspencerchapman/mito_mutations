args <- commandArgs(TRUE)

sample_table <- as.character(args[1])

samples <- read.table(sample_table)

all_files <- list.files(".")
analysis_files <- all_files[grepl(pattern = ".verifyBamID.tsv.selfSM", x = all_files)]

for (i in 1:nrow(samples)) {
  if(!(paste0(samples[i,2],".verifyBamID.tsv.selfSM") %in% analysis_files)){
    system(sprintf('bsub -q normal -e logs/%s.err -o logs/%s.out -R"select[mem>2000]  span[hosts=1] rusage[mem=2000]" -M2000 /lustre/scratch125/casm/team268im/fa8/119/VerifyBamId-farm5/VerifyBamID --Epsilon 1e-12 --BamFile /nfs/cancer_ref01/nst_links/live/%s/%s/%s.sample.dupmarked.bam --Output %s.verifyBamID.tsv --UDPath /lustre/scratch125/casm/team268im/al28/bed_ref/striped_ref_files/ALL_500K.strictmasked.ok.vcf.UD --BedPath /lustre/scratch125/casm/team268im/al28/bed_ref/ALL_500K.strictmasked.ok.vcf.bed --MeanPath /lustre/scratch125/casm/team268im/al28/bed_ref/ALL_500K.strictmasked.ok.vcf.mu --Reference /lustre/scratch125/casm/team268im/al28/bed_ref/striped_ref_files/hs37d5.fa',
                  samples[i,2], samples[i,2], samples[i,1], samples[i,2], samples[i,2], samples[i,2]))
  }
}
