# HDP Flow III: Combine results
# Tim Coorens, Feb 2020

options(stringsAsFactors = F)
library(hdp)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

chlist <- vector("list", 30)
for (i in 1:30){
  if(file.exists(paste0("hdp_chain_",i,".Rdata"))){
    chlist[[i]] <- readRDS(paste0("hdp_chain_",i,".Rdata"))
  }
}
if(any(unlist(lapply(chlist,is.null)))) chlist=chlist[-which(unlist(lapply(chlist,is.null)))]

mut_example_multi <- hdp_multi_chain(chlist)
pdf("QC_plots_chain.pdf") 
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")
dev.off()

mut_example_multi <- hdp_extract_components(mut_example_multi) #This step can take a while. If too long, submit R script as job
saveRDS(mut_example_multi,"HDP_multi_chain.Rdata")

pdf("muts_attributed.pdf")
plot_comp_size(mut_example_multi, bty="L")
dev.off()

mutations=read.table("trinuc_mut_mat.txt")
trinuc_context <- sapply(strsplit(colnames(mutations),split="\\."),function(string) paste0(string[3],string[1],substr(string[4],1,1)))
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16,times=2))
mut_colours=c("dodgerblue","black","red","grey70","olivedrab3","plum2")

#dev.new(width=12,height=4)
#par(mfrow=c(3,4))


for (i in 0:mut_example_multi@numcomp){
  pdf(paste0("hdp_component_",i,".pdf"),width=12,height=4)
  
  plot_comp_distn(mut_example_multi, cat_names=trinuc_context,
                  grouping=group_factor, col=mut_colours,comp=i,
                  col_nonsig="grey80", show_group_labels=TRUE)
  dev.off()
}

plot_dp_comp_exposure(mut_example_multi,
                      dpindices=2:4, incl_numdata_plot=FALSE,
                      col=c(RColorBrewer::brewer.pal(12, "Paired"),"magenta","firebrick"),
                      incl_nonsig=TRUE, cex.names=0.8,
                      ylab_exp = 'Signature exposure', leg.title = 'Signature')

pdf("signature_attribution.pdf",width=10,height=8)

key_table=read.table("key_table.txt")

plot_dp_comp_exposure(mut_example_multi, dpindices=(length(unique(key_table$Patient))+2):length(mut_example_multi@comp_dp_counts), incl_nonsig = T,ylab_exp = 'Signature exposure', leg.title = 'Signature',
                      col=c(RColorBrewer::brewer.pal(12, "Set3"),"magenta","firebrick"))
dev.off()

library(dplyr)
library(ggplot2)
library(tidyr)

mutations=read.table("trinuc_mut_mat.txt")
key_table=read.table("key_table.txt")
sample_remove=rownames(mutations)[rowSums(mutations)<40]
mutations=mutations[!rownames(mutations)%in%sample_remove,]
key_table=key_table[!key_table$Sample%in%sample_remove,]
freq=nrow(mutations)
#freq=table(mutations)

dp_distn <- comp_dp_distn(mut_example_multi)
ndp <- nrow(dp_distn$mean)
ncomp <- ncol(dp_distn$mean)
exposures <- t(dp_distn$mean[length(freq)+1+1:nrow(mutations),,drop=FALSE])
colnames(exposures)=rownames(mutations)

mut_numbers=data.frame(SampleID=rownames(mutations),nmuts=rowSums(mutations))

VAF_groups=unique(sapply(stringr::str_split(colnames(exposures),pattern="_"),function(vec) paste(vec[2:3],collapse="_")))
if(length(VAF_groups)==11) {
  new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")
} else if(length(VAF_groups)==5) {
  new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%",">0.8%")
} else if(length(VAF_groups)==4) {
  new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%",">0.4%")
}

names(new_VAF_groups)<-VAF_groups

mito_sigs_plot<-t(exposures)%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="SampleID")%>%
  left_join(mut_numbers)%>%
  mutate(exp_ID=stringr::str_split(SampleID,pattern="_",simplify=T)[,1])%>%
  mutate(VAF_range=sapply(stringr::str_split(SampleID,pattern="_"),function(vec) paste(vec[2:3],collapse="_")))%>%
  mutate(VAF_range=factor(new_VAF_groups[VAF_range],levels=new_VAF_groups))%>%
  dplyr::select(-SampleID)%>%
  gather(-VAF_range,-exp_ID,-nmuts,key="Signature",value="Exposure")%>%
  mutate(abs_muts=Exposure*nmuts,Signature=paste0("N",Signature))%>%
  mutate(Signature=factor(Signature,levels=c(paste0("N",as.numeric(rownames(exposures)[nrow(exposures)]):3),"N0","N1","N2")))%>%
  ggplot(aes(x=VAF_range,y=abs_muts,fill=Signature))+
  geom_bar(stat="identity",position="stack",col="black")+
  theme_bw()+
  scale_fill_brewer(palette="Paired")+
  theme(axis.text.x = element_text(size=8,angle=90))+
  facet_grid(cols=vars(exp_ID),scales="free_y")

ggsave(filename = "Mut_sigs_by_vaf_range.pdf",plot=mito_sigs_plot,width = 20,height=8)

true_sig="N5"
true_mito_sigs_plot<-t(exposures)%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="SampleID")%>%
  left_join(mut_numbers)%>%
  mutate(exp_ID=stringr::str_split(SampleID,pattern="_",simplify=T)[,1])%>%
  mutate(VAF_range=sapply(stringr::str_split(SampleID,pattern="_"),function(vec) paste(vec[2:3],collapse="_")))%>%
  mutate(VAF_range=factor(new_VAF_groups[VAF_range],levels=new_VAF_groups))%>%
  dplyr::select(-SampleID)%>%
  gather(-VAF_range,-exp_ID,-nmuts,key="Signature",value="Exposure")%>%
  mutate(abs_muts=Exposure*nmuts,Signature=paste0("N",Signature))%>%
  mutate(Signature=factor(Signature,levels=c(paste0("N",as.numeric(rownames(exposures)[nrow(exposures)]):3),"N0","N1","N2")))%>%
  filter(Signature==true_sig)%>%
  ggplot(aes(x=VAF_range,y=abs_muts,fill=Signature))+
  geom_bar(stat="identity",position="stack",col="black")+
  theme_bw()+
  scale_fill_brewer(palette="Paired")+
  theme(axis.text.x = element_text(size=8,angle=90))+
  facet_grid(cols=vars(exp_ID),scales="free_y")

ggsave(filename = "True_sig_by_vaf_range.pdf",plot=true_mito_sigs_plot,width = 20,height=8)

#Get posterior likelihoods of a signature being from any category
components <- comp_categ_distn(mut_example_multi)$mean

sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")

colnames(components)<-c(paste(full_vec,"H",sep="_"),paste(full_vec,"L",sep="_"))
rownames(components)<-paste0("N",0:(nrow(components)-1))
rownames(exposures)<-paste0("N",0:(nrow(exposures)-1))

library(readr)
write_csv(as.data.frame(exposures),"exposures.csv")
write_csv(as.data.frame(components),"components.csv")


