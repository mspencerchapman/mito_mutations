args = commandArgs(TRUE)
bam_file = args[1]

library("GenomicRanges")
library("deepSNV")
library("Rsamtools")

nucleotides = c("A", "T", "C", "G", "-", "INS", "a", "t", "c", "g", "_", "ins")

minPhred=30
samFlag=3844
mapq=30

test.matrix <- matrix(0, ncol = length(nucleotides), nrow = 16569)
colnames(test.matrix) <- nucleotides
for (j in 1:16569) {
  test.matrix[j, ] = bam2R(bam_file, "MT", j, j, q=minPhred, mask=samFlag, mq=mapq)[, nucleotides]
}
mode(test.matrix) = "integer"

write.table(test.matrix, file = paste0(substr(bam_file,1,nchar(bam_file) - 4),"_count.csv"),sep = ",",quote = F,row.names = F,col.names = T)
