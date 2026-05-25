# for i in *; do cd $i; for j in {1..10}; do bsub -q normal -R 'select[mem>=7500] rusage[mem=7500]' -M7500 -e chlist_generation_$j.err -o chlist_generation_$j.out /software/R-3.6.1/bin/Rscript /lustre/scratch125/casm/team268im/al28/mtDNA/al28_code/mtDNA_hdp_tissue_hierarchy_farm_chlist_generation.R $j; done; cd ..; done

args <- commandArgs(TRUE)

iter <- as.numeric(args[1])

library(hdp)

files <- list.files(".")
input_file <- files[grepl("HDPinput",files)]

mut_count <- read.table(input_file, sep = "\t", stringsAsFactors = F)
rownames(mut_count)[which(rownames(mut_count) == "blood_emily_BMH1_TG")] <- "blood_emily_BMH1"
mut_count <- mut_count[which(as.vector(rowSums(mut_count)) != 0),]

tissues <- table(sub("_[^_]+$", "", rownames(mut_count)))

ppindex <- c(0, rep(1,length(tissues)))
cpindex <- c(1, rep(2, length(tissues)))

for(i in 1:length(tissues)) {
  new_pp <- max(ppindex) + 1
  new_cp <- max(cpindex) + 1
  ppindex <- c(ppindex, rep(new_pp, tissues[i]))
  cpindex <- c(cpindex, rep(new_cp, tissues[i]))
}

hdp_mut <- hdp_init(ppindex = ppindex, # index of parental node
                    cpindex = cpindex, # index of the CP to use
                    hh = rep(1, 192), # prior is uniform over 96 categories
                    alphaa = rep(1, length(unique(cpindex))), # shape hyperparameters for 2 CPs
                    alphab = rep(1, length(unique(cpindex))))  # rate hyperparameters for 2 CPs

hdp_mut <- hdp_setdata(hdp_mut, 
                       dpindex = (length(tissues) + 2):numdp(hdp_mut), # index of nodes to add data to
                       mut_count) # input data (mutation counts, sample rows match up with specified dpindex)


hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=iter*200)

chlist <- hdp_posterior(hdp_activated,
                        burnin=500000,
                        n=100,
                        space=1000,
                        cpiter=3,
                        seed=iter*1e3) 
saveRDS(chlist, paste0("chlist_",iter,".rds"))
