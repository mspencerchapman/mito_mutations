library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(phylosignal)
library(GenomicRanges)
library(Rsamtools)
options(stringsAsFactors = F)

#Set these file paths before running the script
genomeFile="~/R_work/reference_files/genome.fa"
root_dir="~/R_work/mito_mutations_blood"
source(paste0(root_dir,"/data/mito_mutations_blood_functions.R"))

#Set the plotting theme for ggplot2
my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))+
  theme(legend.key.height=unit(3,"mm"),
        legend.title = element_text(size=8))

#Set the key file paths using the root dir
tree_file_paths = list.files(paste0(root_dir,"/data/tree_files"),pattern=".tree",full.names = T)
ref_file=paste0(root_dir,"/data/Samples_metadata_ref.csv")
raw_data_adult_folder=paste0(root_dir,"/data/blood_adult/")
raw_data_foetal_folder=paste0(root_dir,"/data/blood_foetal/")

#Get individual level metadata
ref_file=paste0(root_dir,"/data/Samples_metadata_ref.csv")
mito_data_file=paste0(root_dir,"/data/mito_data.Rds")
ref_df<-read.csv(ref_file)

#Read in the main "mito_data" object & list of "CN correlating muts" that get excluded from analysis
mito_data<-readRDS(mito_data_file)
CN_correlating_muts<-readRDS(paste0(root_dir,"/data/CN_correlation.RDS"))
length(CN_correlating_muts)
exclude_muts=c("MT_16183_A_C") #This is a specific artefact primarily in the 8pcw foetus due to a germline SNP and consequent frequent PCR slippage

##----------------PART 1: PRE-SIGNATURE EXTRACTION----------------
##DIVIDE MUTATION SETS INTO VAF BINS PER SAMPLE

#Define the VAF bins to use for signature extraction
bins=2^c(-10:0)

#Divide the mutations set by VAF/ SAMPLE
mutCN_cutoff=25 #CN correlating muts with a mutation copy number less than this cut-off are censored as NUMTs/ contamination
df_tidy_full<-dplyr::bind_rows(Map(list=mito_data,exp_ID=names(mito_data),function(list,exp_ID) {
  if(is.null(list)){stop(return(NULL))}
  CN_correlating_mut_removal_mat=list$matrices$implied_mutCN>mutCN_cutoff|(matrix((!rownames(list$matrices$vaf)%in%CN_correlating_muts),ncol=1)%*%matrix(rep(1,ncol(list$matrices$vaf)),nrow=1))
  implied_mut_CN_tidy<-list$matrices$implied_mutCN%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var="mut_ref")%>%
    dplyr::select(-global)%>%
    tidyr::gather(key="Sample",value="implied_mut_CN",-mut_ref)
  
  df_tidy<-(list$matrices$vaf*list$matrices$SW*CN_correlating_mut_removal_mat)%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var="mut_ref")%>%
    mutate(rho_val=list$rho_vals)%>%
    dplyr::select(-global)%>%
    tidyr::gather(key="Sample",value="vaf",-mut_ref,-rho_val)
  
  comb_tidy<-left_join(df_tidy,implied_mut_CN_tidy,by=c("Sample","mut_ref"))%>%
    dplyr::filter(!grepl("DEL|INS",mut_ref) & !mut_ref%in%exclude_muts)%>%
    dplyr::filter(vaf>0)%>%
    mutate(exp_ID=exp_ID)
  return(comb_tidy)
}))

all_muts<-lapply(names(mito_data),function(this_exp_ID) {
  print(this_exp_ID)
  exp_bins<-lapply(1:length(bins),function(i) {
    lower_limit=ifelse(i==1,0,bins[i-1])
    upper_limit=bins[i]
    print(paste(lower_limit,"to",upper_limit))
    bin_exp_muts<-df_tidy_full%>%
      dplyr::filter(exp_ID==this_exp_ID & vaf>lower_limit & vaf<= upper_limit)%>%
      dplyr::select(mut_ref)%>%
      tidyr::separate(mut_ref,into=c("chr","pos","ref","mut"))%>%
      mutate(donor=paste(this_exp_ID,lower_limit,upper_limit,sep = "_"))%>%
      dplyr::filter(!duplicated(.))
    
    print(nrow(bin_exp_muts))
    return(bin_exp_muts)
  })
  return(exp_bins)
})

#Bind together mutation sets & annotate the trinucleotide context
mutations<-dplyr::bind_rows(all_muts)%>%
  mutate(pos=as.numeric(pos))
mutations$trinuc_ref = as.vector(Rsamtools::scanFa(genomeFile, GRanges(mutations$chr, IRanges(mutations$pos-1, mutations$pos+1))))

#Annotate the mutation from the pyrimidine base
ntcomp = c(T="A",G="C",C="G",A="T")
mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
mutations$trinuc_ref_py = mutations$trinuc_ref
for (j in 1:nrow(mutations)) {
  if (mutations$ref[j] %in% c("A","G")) { # Purine base
    mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
    mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
  }
}

mutation_profiles<-lapply(unique(mutations$donor),function(SampleID) {
  sample_mutations<-mutations%>%dplyr::filter(donor==SampleID)
  
  #Counting subs
  freqs_heavy = table(paste(sample_mutations$sub[which(sample_mutations$ref %in% c("A","G"))],paste(substr(sample_mutations$trinuc_ref_py[which(sample_mutations$ref %in% c("A","G"))],1,1),substr(sample_mutations$trinuc_ref_py[which(sample_mutations$ref %in% c("A","G"))],3,3),sep="-"),sep=","))
  freqs_light = table(paste(sample_mutations$sub[which(sample_mutations$ref %in% c("C","T"))],paste(substr(sample_mutations$trinuc_ref_py[which(sample_mutations$ref %in% c("C","T"))],1,1),substr(sample_mutations$trinuc_ref_py[which(sample_mutations$ref %in% c("C","T"))],3,3),sep="-"),sep=","))
  
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_heavy_full = freqs_heavy[full_vec]; freqs_heavy_full[is.na(freqs_heavy_full)] = 0; names(freqs_heavy_full) = paste(full_vec,"H",sep="_")
  freqs_light_full = freqs_light[full_vec]; freqs_light_full[is.na(freqs_light_full)] = 0; names(freqs_light_full) = paste(full_vec,"L",sep="_")
  
  output=matrix(c(freqs_heavy_full,freqs_light_full),nrow=1,dimnames = list(SampleID,c(names(freqs_heavy_full),names(freqs_light_full))))
  return(output)
})

mutation_profiles_mat<-Reduce(rbind,mutation_profiles)
dim(mutation_profiles_mat)

key_table=data.frame(Sample=rownames(mutation_profiles_mat))%>%
  dplyr::mutate(Patient=stringr::str_split(Sample,pattern = "_",simplify=T)[,1])%>%
  dplyr::select(Patient,Sample)

write.table(mutation_profiles_mat,paste0(root_dir,"/data/mutational_signatures/trinuc_mut_mat.txt"))
write.table(key_table,paste0(root_dir,"/data/mutational_signatures/key_table.txt"))

##NOW PERFORM SIGNATURE EXTRACTION USING HDP - THIS IS DONE IN THE COMMAND LINE SO THAT CAN BE DONE IN PARALLEL

##----------------PART 2a: POST-SIGNATURE EXTRACTION (HDP)----------------
###Go back to local drive
exposures<-read.csv(paste0(root_dir,"/data/mutational_signatures/exposures.csv"))
components<-read.csv(paste0(root_dir,"/data/mutational_signatures/components.csv"))

sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
colnames(components)<-c(paste(full_vec,"H",sep="_"),paste(full_vec,"L",sep="_"))
rownames(components)<-paste0("N",0:(nrow(components)-1))
rownames(exposures)<-paste0("N",0:(nrow(exposures)-1))
colnames(exposures)<-gsub("^X","",colnames(exposures))
colnames(exposures)<-gsub("\\.pcw","pcw",colnames(exposures))

mutations$mut_type<-sapply(1:nrow(mutations),function(i) {
  chain<-ifelse(mutations$ref[i]%in%c("A","G"),"H","L")
  mut_type=paste0(mutations$sub[i],",",substr(mutations$trinuc_ref_py[i],1,1),"-",substr(mutations$trinuc_ref_py[i],3,3),"_",chain)
  return(mut_type)
})

sig.post.probs<-lapply(1:nrow(mutations),function(i) {
  #cat(i,sep="\n")
  sample_cat<-mutations$donor[i]
  mut_type<-mutations$mut_type[i]
  
  if(!sample_cat%in%colnames(exposures)){
    exp_ID=stringr::str_split(sample_cat,"_",simplify=T)[,1]
    vaf_bin_min=stringr::str_split(sample_cat,"_",simplify=T)[,2]
    sample_cat<-colnames(exposures)[max(grep(exp_ID,colnames(exposures)))]
  }
  
  if(!mut_type%in%colnames(components)){
    res=rep(NA,nrow(components)+1)
    names(res)<-c(rownames(components),"ML_Sig")
    return(res)
  } else {
    likelihoods=exposures[,sample_cat]*components[,mut_type]
    
    post.probs=likelihoods/sum(likelihoods)
    ML_sig<-rownames(components)[which.max(post.probs)]
    res<-c(post.probs,ML_sig)
    names(res)<-c(rownames(components),"ML_Sig")
    return(res)
  }
  
})%>%bind_rows()

mutations<-bind_cols(mutations,sig.post.probs)
sig_ref<-mutations%>%
  tidyr::separate(donor,into=c("exp_ID","lower_VAF","upper_VAF"),sep="_",remove=F)%>%
  mutate_at(c("lower_VAF","upper_VAF"),as.numeric)%>%
  tidyr::unite(col="mut_ref",chr,pos,ref,mut,remove=F)%>%
  dplyr::select(mut_ref,exp_ID,lower_VAF,upper_VAF,ML_Sig)

saveRDS(sig_ref,file=paste0(root_dir,"/data/mutational_signatures/sig_ref_file.Rds"))

sigs_1<-mitochondrial_extracted_signature_plot(components[1:4,])
sigs_2<-mitochondrial_extracted_signature_plot(components[5:7,])
real_sig<-mitochondrial_extracted_signature_plot(components[2,])

##----------------SUMMARY PLOTS OF EXTRACTION RESULTS (HDP)----------------
sig_ref<-readRDS(file=paste0(root_dir,"/data/mutational_signatures/sig_ref_file.Rds"))
VAF_groups=unique(sapply(stringr::str_split(colnames(exposures),pattern="_"),function(vec) paste(vec[2:3],collapse="_")))
new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")
names(new_VAF_groups)<-VAF_groups

mut_numbers=data.frame(SampleID=rownames(mutation_profiles_mat),nmuts=rowSums(mutation_profiles_mat))

#Fill in the bins with too low mutation numbers - assume all are N1, as these are high VAF mutations
excluded_cats<-rownames(mutation_profiles_mat)[!rownames(mutation_profiles_mat)%in%colnames(exposures)]
all_N1=c(0,1,rep(0,nrow(exposures)-2))
excluded_mat<-matrix(rep(all_N1,times=length(excluded_cats)),nrow=length(all_N1),dimnames = list(rownames(exposures),excluded_cats))

mito_sigs_plot<-t(cbind(exposures,excluded_mat))%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="SampleID")%>%
  left_join(mut_numbers)%>%
  mutate(exp_ID=stringr::str_split(SampleID,pattern="_",simplify=T)[,1])%>%
  mutate(exp_ID=factor(exp_ID,levels=c("8pcw","18pcw","CB001" ,"CB002", "KX001", "KX002", "SX001","AX001","KX007", "KX008", "KX004", "KX003")))%>%
  mutate(VAF_range=sapply(stringr::str_split(SampleID,pattern="_"),function(vec) paste(vec[2:3],collapse="_")))%>%
  mutate(VAF_range=factor(new_VAF_groups[VAF_range],levels=new_VAF_groups))%>%
  dplyr::select(-SampleID)%>%
  gather(-VAF_range,-exp_ID,-nmuts,key="Signature",value="Exposure")%>%
  mutate(abs_muts=Exposure*nmuts)%>%
  mutate(Signature=factor(Signature,levels=c(rownames(exposures)[nrow(exposures):3],"N0","N1")))%>%
  ggplot(aes(x=VAF_range,y=abs_muts,fill=Signature))+
  geom_bar(stat="identity",position="stack",col="black",size=0.25)+
  theme_bw()+
  scale_fill_brewer(palette="Paired")+
  my_theme+
  theme(axis.text.x = element_text(size=5,angle=90),strip.text.x = element_text(size=7,margin = unit(c(1,0,1,0),"mm")))+
  facet_wrap(~exp_ID,ncol=6)+
  labs(x="VAF range",y="Number of mutations assigned")

muts_per_bin_per_sample<-t(cbind(exposures,excluded_mat))%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="SampleID")%>%
  left_join(mut_numbers)%>%
  mutate(exp_ID=stringr::str_split(SampleID,pattern="_",simplify=T)[,1])%>%
  mutate(exp_ID=factor(exp_ID,levels=c("8pcw","18pcw","CB001" ,"CB002", "KX001", "KX002", "SX001","AX001","KX007", "KX008", "KX004", "KX003")))%>%
  mutate(VAF_range=sapply(stringr::str_split(SampleID,pattern="_"),function(vec) paste(vec[2:3],collapse="_")))%>%
  mutate(VAF_range=factor(new_VAF_groups[VAF_range],levels=new_VAF_groups))%>%
  dplyr::select(-SampleID)%>%
  gather(-VAF_range,-exp_ID,-nmuts,key="Signature",value="Exposure")%>%
  mutate(abs_muts=Exposure*nmuts)

##----------------PART 2b: POST-SIGNATURE EXTRACTION (SigProfiler)----------------
###Go back to local drive
exposures_SP<-read.delim(paste0(root_dir,"/data/mutational_signatures/SigProfiler/CH192_S7_NMF_Activities.txt"))%>%tibble::column_to_rownames(var="Samples")%>%t()
components_SP<-read.delim(paste0(root_dir,"/data/mutational_signatures/SigProfiler/CH192_S7_Signatures.txt"))%>%tibble::column_to_rownames(var="MutationType")%>%t()%>%as.data.frame()

sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
colnames(components_SP)<-c(paste(full_vec,"H",sep="_"),paste(full_vec,"L",sep="_"))
rownames(components_SP)<-paste0("N",0:(nrow(components_SP)-1))

exposures_SP<-apply(exposures_SP,2,function(x) x/sum(x))
rownames(exposures_SP)<-paste0("N",0:(nrow(exposures_SP)-1))
colnames(exposures_SP)<-gsub("^X","",colnames(exposures_SP))
colnames(exposures_SP)<-gsub("\\ pcw","pcw",colnames(exposures_SP))

mutations$mut_type<-sapply(1:nrow(mutations),function(i) {
  chain<-ifelse(mutations$ref[i]%in%c("A","G"),"H","L")
  mut_type=paste0(mutations$sub[i],",",substr(mutations$trinuc_ref_py[i],1,1),"-",substr(mutations$trinuc_ref_py[i],3,3),"_",chain)
  return(mut_type)
})

sig.post.probs_SP<-lapply(1:nrow(mutations),function(i) {
  #cat(i,sep="\n")
  sample_cat<-mutations$donor[i]
  mut_type<-mutations$mut_type[i]
  
  if(!sample_cat%in%colnames(exposures_SP)){
    exp_ID=stringr::str_split(sample_cat,"_",simplify=T)[,1]
    vaf_bin_min=stringr::str_split(sample_cat,"_",simplify=T)[,2]
    sample_cat<-colnames(exposures_SP)[max(grep(exp_ID,colnames(exposures_SP)))]
  }
  
  if(!mut_type%in%colnames(components_SP)){
    res=rep(NA,nrow(components_SP)+1)
    names(res)<-c(rownames(components_SP),"ML_Sig")
    return(res)
  } else {
    likelihoods=exposures_SP[,sample_cat]*components_SP[,mut_type]
    
    post.probs=likelihoods/sum(likelihoods)
    ML_sig<-rownames(components_SP)[which.max(post.probs)]
    res<-c(post.probs,ML_sig)
    names(res)<-c(rownames(components_SP),"ML_Sig")
    return(res)
  }
  
})%>%bind_rows()

mutations_SP<-bind_cols(mutations,sig.post.probs_SP)
sig_ref_SP<-mutations_SP%>%
  tidyr::separate(donor,into=c("exp_ID","lower_VAF","upper_VAF"),sep="_",remove=F)%>%
  mutate_at(c("lower_VAF","upper_VAF"),as.numeric)%>%
  tidyr::unite(col="mut_ref",chr,pos,ref,mut,remove=F)%>%
  dplyr::select(mut_ref,exp_ID,lower_VAF,upper_VAF,ML_Sig)

saveRDS(sig_ref_SP,file=paste0(root_dir,"/data/mutational_signatures/SigProfiler/sig_ref_file_SP.Rds"))

sigs_1<-mitochondrial_extracted_signature_plot(components_SP[1:4,])
sigs_2<-mitochondrial_extracted_signature_plot(components_SP[5:7,])
real_sig<-mitochondrial_extracted_signature_plot(components_SP[1,])

comb_sig<-gridExtra::arrangeGrob(grobs=list(sigs_1,sigs_2),ncol=2)
ggsave(filename = paste0(root_dir,"/rebuttal_plots/SigProfiler_sigs.pdf"),plot = comb_sig,width=7,height=3)

pdf(file=paste0(root_dir,"/rebuttal_plots/HDP_vs_SigProfiler_cosine_sim.pdf"),height=3,width=3)
pheatmap(MutationalPatterns::cos_sim_matrix(mut_matrix1 = t(components),mut_matrix2 = t(components_SP)),cluster_rows = F,cluster_cols = F,scale = "none")
dev.off()

##----------------SUMMARY PLOTS OF EXTRACTION RESULTS (SigProfiler) ----------------
sig_ref_SP<-readRDS(file=paste0(root_dir,"/data/mutational_signatures/SigProfiler/sig_ref_file_SP.Rds"))
VAF_groups=unique(sapply(stringr::str_split(colnames(exposures_SP),pattern="_"),function(vec) paste(vec[2:3],collapse="_")))
new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")
names(new_VAF_groups)<-VAF_groups

mut_numbers=data.frame(SampleID=rownames(mutation_profiles_mat),nmuts=rowSums(mutation_profiles_mat))

#Fill in the bins with too low mutation numbers - assume all are N1, as these are high VAF mutations
excluded_cats<-rownames(mutation_profiles_mat)[!rownames(mutation_profiles_mat)%in%colnames(exposures_SP)]
all_N0=c(1,0,rep(0,nrow(exposures_SP)-2))
excluded_mat<-matrix(rep(all_N0,times=length(excluded_cats)),nrow=length(all_N0),dimnames = list(rownames(exposures_SP),excluded_cats))

mito_sigs_plot_SP<-t(cbind(exposures_SP,excluded_mat))%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="SampleID")%>%
  left_join(mut_numbers)%>%
  mutate(exp_ID=stringr::str_split(SampleID,pattern="_",simplify=T)[,1])%>%
  mutate(exp_ID=factor(exp_ID,levels=c("8pcw","18pcw","CB001" ,"CB002", "KX001", "KX002", "SX001","AX001","KX007", "KX008", "KX004", "KX003")))%>%
  mutate(VAF_range=sapply(stringr::str_split(SampleID,pattern="_"),function(vec) paste(vec[2:3],collapse="_")))%>%
  mutate(VAF_range=factor(new_VAF_groups[VAF_range],levels=new_VAF_groups))%>%
  dplyr::select(-SampleID)%>%
  gather(-VAF_range,-exp_ID,-nmuts,key="Signature",value="Exposure")%>%
  mutate(abs_muts=Exposure*nmuts)%>%
  mutate(Signature=factor(Signature,levels=c(rownames(exposures_SP)[nrow(exposures_SP):3],"N1","N0")))%>%
  ggplot(aes(x=VAF_range,y=abs_muts,fill=Signature))+
  geom_bar(stat="identity",position="stack",col="black",size=0.25)+
  theme_bw()+
  scale_fill_brewer(palette="Paired")+
  my_theme+
  theme(axis.text.x = element_text(size=5,angle=90),strip.text.x = element_text(size=7,margin = unit(c(1,0,1,0),"mm")))+
  facet_wrap(~exp_ID,ncol=6)+
  labs(x="VAF range",y="Number of mutations assigned")

ggsave(filename = paste0(root_dir,"/rebuttal_plots/SigProfiler_contributions_by_VAF_bin.pdf"),plot = mito_sigs_plot_SP,width=7,height=3.5)

#Compare the numbers of 'real' mutations in each bin between SigProfiler & HDP
muts_per_bin_per_sample_SP<-t(cbind(exposures_SP,excluded_mat))%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="SampleID")%>%
  left_join(mut_numbers)%>%
  mutate(exp_ID=stringr::str_split(SampleID,pattern="_",simplify=T)[,1])%>%
  mutate(exp_ID=factor(exp_ID,levels=c("8pcw","18pcw","CB001" ,"CB002", "KX001", "KX002", "SX001","AX001","KX007", "KX008", "KX004", "KX003")))%>%
  mutate(VAF_range=sapply(stringr::str_split(SampleID,pattern="_"),function(vec) paste(vec[2:3],collapse="_")))%>%
  mutate(VAF_range=factor(new_VAF_groups[VAF_range],levels=new_VAF_groups))%>%
  dplyr::select(-SampleID)%>%
  gather(-VAF_range,-exp_ID,-nmuts,key="Signature",value="Exposure")%>%
  mutate(abs_muts=Exposure*nmuts)

HDP_Sigprofiler_correlation<-bind_rows(muts_per_bin_per_sample_SP%>%filter(Signature=="N0" &!is.na(nmuts))%>%mutate(Extraction="SigProfiler"),
          muts_per_bin_per_sample%>%filter(Signature=="N1")%>%mutate(Extraction="HDP"))%>%
  pivot_wider(id_cols=c("exp_ID","VAF_range"),names_from = "Extraction",values_from = "abs_muts")%>%
  ggplot(aes(x=SigProfiler,y=HDP))+
  geom_point(size=0.3)+
  geom_smooth(method="lm",linetype=2,col="blue",linewidth=0.5)+
  geom_abline(linetype=1,col="black")+
  theme_classic()+
  my_theme+
  ggpmisc::stat_poly_eq(ggpmisc::use_label("eq"))
  

ggsave(filename = paste0(root_dir,"/rebuttal_plots/HDP_Sigprofiler_correlation.pdf"),plot = HDP_Sigprofiler_correlation,width=2.5,height=2)


##----------------PRODUCE THE 'ML_Sig' MATRIX----------------
#Enhance data with a matched "ML_Sig" matrix matching the VAF/ SW matrixes in the mito_data object
sig_ref<-readRDS(paste0(root_dir,"/data/mutational_signatures/sig_ref_file.Rds"))

#The "MT_16519_T_T" mutation is real, but is at a site where there is a germline SNP (MT_16519_T_C)
#Therefore actual somatic mutation is a MT_16519_C_T - need to ensure this is marked as the real "N1" signature so that is retained in mutaiton set
sig_ref<-rbind(sig_ref,data.frame(mut_ref=rep("MT_16519_T_T",6),
                                  exp_ID=rep("KX003",6),
                                  lower_VAF=bins[5:10],
                                  upper_VAF=bins[6:11],
                                  ML_Sig=rep("N1",6)))

if(is.null(mito_data[[1]]$matrices$ML_Sig)) {
  cat("Adding the maximum likelihood mutational signature matrix for all mutations/ samples.",sep="\n")
  #Annotate specific mutations with their most likely signature using the 'sig_ref' dataframe
  mito_data<-Map(list=mito_data,this_exp_ID=names(mito_data),function(list,this_exp_ID) {
    cat(this_exp_ID,sep="\n")
    if(!is.null(list$matrices$ML_Sig)) {
      cat("ML_sig matrix already exists",sep="\n")
      return(list)
    } else {
      #this_exp_ID<-gsub("8pcw","8 pcw",this_exp_ID)
      mat<-list$matrices$vaf*list$matrices$SW
      sig_list<-lapply(1:nrow(mat),function(i) {
        ref<-sig_ref%>%filter(mut_ref==rownames(mat)[i] & exp_ID==this_exp_ID)
        mut_res<-mat[i,]
        for(j in 1:nrow(ref)){
          mut_res[mut_res>ref$lower_VAF[j] & mut_res<=ref$upper_VAF[j]]<-ref$ML_Sig[j]
        }
        return(mut_res%>%mutate_all(as.character))
      })
      sig_mat<-sig_list%>%bind_rows()%>%as.matrix()
      list$matrices$ML_Sig<-sig_mat
      return(list)
    }
    
  })
  cat("Analysis completed. Saving updated mito_data file.",sep="\n")
  saveRDS(mito_data,file=mito_data_file)
}

##---------------PRODUCE A "VAF.FILT" MATRIX-------------------
#Finally, combine information from Shearwater, mutational signature extraction and whether mutation may be contaminant to
#create a matrix of mutation VAFs but with likely artefacts and contaminants excluded (i.e. VAF set to 0)

if(is.null(mito_data[[1]]$matrices$vaf.filt)) {
  cat("Adding the filtered VAF matrix.",sep="\n")
  #Annotate specific mutations with their most likely signature using the 'sig_ref' dataframe
  mito_data<-Map(list=mito_data,this_exp_ID=names(mito_data),function(list,this_exp_ID) {
    cat(this_exp_ID,sep="\n")
    
    mutCN_cutoff=25 #If the mitochondrial copy number is over 25, retain mutation even if is in the "CN correlating muts" list
    CN_correlating_mut_removal_mat=list$matrices$implied_mutCN>mutCN_cutoff|(matrix((!rownames(list$matrices$vaf)%in%CN_correlating_muts),ncol=1)%*%matrix(rep(1,ncol(list$matrices$vaf)),nrow=1))
    vaf.filt<-(list$matrices$vaf[,list$tree$tip.label]*list$matrices$SW[,list$tree$tip.label]*(list$matrices$ML_Sig=="N1")[,list$tree$tip.label]*CN_correlating_mut_removal_mat[,list$tree$tip.label])
    
    list$matrices$CN_correlating_mut_removal_mat<-CN_correlating_mut_removal_mat
    list$matrices$vaf.filt<-vaf.filt
    return(list)
  })
  cat("Analysis completed. Saving updated mito_data file.",sep="\n")
  saveRDS(mito_data,file=mito_data_file)
}
