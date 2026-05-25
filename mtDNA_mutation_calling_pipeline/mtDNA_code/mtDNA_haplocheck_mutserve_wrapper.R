args <- commandArgs(TRUE)

sample_table <- as.character(args[1])
sample_type <- substr(sample_table,1,nchar(sample_table) - 4)

samples <- read.table(sample_table)

vcf_files <- list.files(paste0("/lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/haplocheck/",sample_type))
bam_files <- list.files(paste0("/lustre/scratch125/casm/team268im/al28/mtDNA/Data/",sample_type))

for (i in 1:nrow(samples)) {
  if(!(paste0(samples[i,2],".vcf.gz") %in% vcf_files) & (paste0(samples[i,2],"_MT.bam") %in% bam_files)){
    system(sprintf("bsub -q normal -o /lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/haplocheck/%s/logs/%s.mutserve.out -e /lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/haplocheck/%s/logs/%s.mutserve.err -n1 -M4000 -R'select[mem>4000] rusage[mem=4000] span[hosts=1]' java -Xmx4G -jar /lustre/scratch125/casm/team268im/al28/software/mutserve.jar call --reference /lustre/scratch125/casm/team268im/al28/software/rCRS.fasta --output /lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/haplocheck/%s/%s.vcf.gz /lustre/scratch125/casm/team268im/al28/mtDNA/Data/%s/%s_MT.bam",
                   sample_type, samples[i,2], sample_type, samples[i,2], sample_type, samples[i,2], sample_type, samples[i,2]))
  }
}
