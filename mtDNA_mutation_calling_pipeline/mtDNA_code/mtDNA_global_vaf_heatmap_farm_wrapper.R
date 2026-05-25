args <- commandArgs(TRUE)

tissue <- as.character(args[1])

system(sprintf("bsub -q long -M5000 -R'select[mem>5000] rusage[mem=5000] span[hosts=1]' -n1 -e /lustre/scratch125/casm/team268im/al28/mtDNA/Data/%s/logs/%s.err -o /lustre/scratch125/casm/team268im/al28/mtDNA/Data/%s/logs/%s.out Rscript /lustre/scratch126/casm/team154pc/ms56/my_programs/mtDNA_code/mtDNA_global_vaf_heatmap_farm.R %s",
               tissue, tissue, tissue, tissue, tissue))
