genome_file <- "/Users/al28/Mounts/cancer/bed_ref/hs37d5.fa"

regions <- as.data.frame(rbind(c("MT",577,16023),c("MT",1,576),c("MT",16024,16569)))
colnames(regions) <- c("chr","start","end")

regions$chr <- as.character(regions$chr)
regions$start <- as.numeric(as.character(regions$start))
regions$end <- as.numeric(as.character(regions$end))

tri_count <- array(data = NA, dim = c(3, 64))

for(i in 1:nrow(regions)) {
  sequence = as.vector(scanFa(genome_file, GRanges(regions$chr[i], IRanges(regions$start[i], regions$end[i]))))
  tri_count_temp <- t(sapply(sequence,function(x) trinucleotideFrequency( DNAString(x) )))
  tri_count[i,] <- tri_count_temp
}

tri_count <- rbind(tri_count,tri_count[2,] + tri_count[3,])
tri_count <- tri_count[-c(2,3),]

rownames(tri_count) <- c("coding","d_loop")
colnames(tri_count) <- colnames(tri_count_temp)

tri_count["coding","ACT"] <- tri_count["coding","ACT"] + 1
tri_count["coding","CTT"] <- tri_count["coding","CTT"] + 1
tri_count["coding","AGT"] <- tri_count["coding","AGT"] + 1
tri_count["coding","TGT"] <- tri_count["coding","TGT"] + 1

tri_count["d_loop","GTT"] <- tri_count["d_loop","GTT"] + 1
tri_count["d_loop","TGG"] <- tri_count["d_loop","TGG"] + 1
tri_count["d_loop","GGA"] <- tri_count["d_loop","GGA"] + 1
tri_count["d_loop","CAG"] <- tri_count["d_loop","CAG"] + 1

tri_count_py <- data.frame(array(data = NA, dim = c(4,32)))
colnames(tri_count_py) <- c("ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT",
                            "ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT")
rownames(tri_count_py) <- c("coding_heavy","coding_light","d_loop_heavy","d_loop_light")

rev_comp_trinuc <- c("TGT","GGT","CGT","AGT","TGG","GGG","CGG","AGG","TGC","GGC","CGC","AGC","TGA","GGA","CGA","AGA",
                     "TAT","GAT","CAT","AAT","TAG","GAG","CAG","AAG","TAC","GAC","CAC","AAC","TAA","GAA","CAA","AAA")

for(i in 1:ncol(tri_count_py)){
  tri_count_py["coding_light",colnames(tri_count_py)[i]] <- tri_count["coding",colnames(tri_count_py)[i]]
  tri_count_py["coding_heavy",colnames(tri_count_py)[i]] <- tri_count["coding",rev_comp_trinuc[i]]
  tri_count_py["d_loop_light",colnames(tri_count_py)[i]] <- tri_count["d_loop",colnames(tri_count_py)[i]]
  tri_count_py["d_loop_heavy",colnames(tri_count_py)[i]] <- tri_count["d_loop",rev_comp_trinuc[i]]
}

saveRDS(tri_count_py, "/Users/al28/Mounts/cancer/mtDNA/trinuc_freqs_coding_dloop_heavy_light.Rds")

