args <- commandArgs(TRUE)

tissue <- as.character(args[1])

system(sprintf("bsub -q normal -M3000 -R'select[mem>3000] rusage[mem=3000] span[hosts=1]' -n1 -e /lustre/scratch125/casm/team268im/al28/mtDNA/Data/%s/logs/%s.err -o /lustre/scratch125/casm/team268im/al28/mtDNA/Data/%s/logs/%s.out /software/R-3.6.1/bin/Rscript /lustre/scratch125/casm/team268im/al28/mtDNA/al28_code/mtDNA_global_vaf_heatmap_exome_farm.R %s",
               tissue, tissue, tissue, tissue, tissue))
