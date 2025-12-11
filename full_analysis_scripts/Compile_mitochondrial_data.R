library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(phylosignal)
options(stringsAsFactors = F)

#Set these file paths before running the script
genomeFile="~/Documents/Reference_files/hs37d5.fa"
root_dir="~/R_work/mito_mutations_blood/"
source(paste0(root_dir,"/data/mito_mutations_blood_functions.R"))

#Set the key file paths using the root dir
tree_file_paths = list.files(paste0(root_dir,"/data/tree_files"),pattern=".tree",full.names = T)
ref_file=paste0(root_dir,"/data/Samples_metadata_ref.csv")
raw_data_adult_folder=paste0(root_dir,"/data/blood_adult/")
raw_data_foetal_folder=paste0(root_dir,"/data/blood_foetal/")

#Get individual level metadata
ref_file=paste0(root_dir,"/data/Samples_metadata_ref.csv")
mito_data_file=paste0(root_dir,"/data/mito_data.Rds")
ref_df<-read.csv(ref_file)
figures_dir=paste0(root_dir,"/figures/")

##---------------------COMPILE THE ADULT DATA---------------------
#Samples from initial sequencing analysis that showed high levels of specific artefacts & may skew attribution algorithms
high_artefact_samples=c("PD43974af2","PD43974ag2","PD43974ck2","PD43974bg2","PD43974am2","PD43974an2","PD41048b_lo0061")

mito_data_adult<-lapply(ref_df%>%filter(Dataset=="Adult_Cord")%>%pull(Sample),function(exp_ID) {
  print(exp_ID)
  PDID=ref_df$PDID[ref_df$Sample==exp_ID]
  if(PDID=="PD43976"){PDID<-"BMH1_TG"} #AX001 files are named with old ID
  res<-generate_mito_matrices(PD_number = PDID,
                              tree_file_path = grep(exp_ID,tree_file_paths,value = T),
                              shearwater_calls_file=paste0(raw_data_adult_folder,"shearwater_calls/",PDID,"_shearwater_mutations.tsv"),
                              pileup_folder=paste0(raw_data_adult_folder,"pileup_count/"),
                              run_bb = T,
                              exclude_samples = high_artefact_samples)
  return(res)
})
names(mito_data_adult)<-ref_df%>%filter(Dataset=="Adult_Cord")%>%pull(Sample)

#Add the ultrametric trees to the data lists
mito_data_adult<-lapply(mito_data_adult,function(list) {
  tree.ultra<-make.ultrametric.tree(list$tree)
  tree.ultra<-add_ancestral_outgroup(tree.ultra)
  tree.ultra$coords<-NULL
  tree.ultra$edge.length<-tree.ultra$edge.length*mean(get_mut_burden(list$tree))
  list$tree.ultra<-tree.ultra
  list$tree.ultra<-plot_tree(list$tree.ultra)
  return(list)
})

##---------------------COMPILE THE FOETAL DATA---------------------
mito_data_foetal=lapply(ref_df%>%filter(Dataset=="Foetal")%>%pull(Sample),function(exp_ID) {
  print(exp_ID)
  PDID=ref_df$PDID[ref_df$Sample==exp_ID]
  print(PDID)
  mito_data=generate_mito_matrices(PD_number = PDID,
                                   tree_file_path = grep(paste0("_",exp_ID),tree_file_paths,value = T),
                                   pileup_folder = paste0(raw_data_foetal_folder,"pileup_counts/"),
                                   shearwater_calls_file = paste0(raw_data_foetal_folder,"shearwater_calls/",PDID,"_shearwater_mutations.tsv"),
                                   rev_germline = T,
                                   run_bb = T)
  cat(paste("Generate the ultrametric tree for",exp_ID),sep="\n")
  mito_data$tree.ultra<-make.ultrametric.tree(mito_data$tree)
  mito_data$tree.ultra$edge.length<-mito_data$tree.ultra$edge.length*mean(get_mut_burden(mito_data$tree))
  mito_data$tree.ultra<-plot_tree(mito_data$tree.ultra,cex.label = 0,plot_axis = F)
  return(mito_data)
})
names(mito_data_foetal)<-ref_df%>%filter(Dataset=="Foetal")%>%pull(Sample)

##---------------------COMBINE FOETAL AND ADULT INTO SINGLE LIST---------------------
mito_data<-c(mito_data_adult,mito_data_foetal)

##-----------ADD PHYLOGENETIC SIGNAL ANALYSIS FOR EACH MITOCHONDRIAL MUTATION-----------
cat("Starting calculation of phylogenetic signal for all mutations (this can take >12 hours)",sep="\n")
mito_data<-Map(list=mito_data,Exp_ID=names(mito_data),function(list,Exp_ID) {
  cat(Exp_ID,sep="\n")
  if(is.null(list$phylosignal)) {
    tree.noancestral<-drop.tip(list$tree.ultra,"Ancestral")
    tree.4d<-phylobase::phylo4d(tree.noancestral,tip.data=t(list$matrices$vaf[,tree.noancestral$tip.label]))
    list$phylosignal<-phyloSignal(tree.4d)
  }
  return(list)
})

#SAVE THE COMBINED FILE
saveRDS(mito_data,file=mito_data_file)

##------------ENHANCE WITH A "implied_mutCN" MATRIX PER SAMPLE/ MUTATION-----------

#Read in the copy number data
mito_cn=read.csv(paste0(root_dir,"/data/whole_genome_coverage_pileup_and_bedtools.csv"),header=T)

#This is the calculated absolute number of mutated mtDNA copies per cell, using the VAF and mtDNA copy number
if(is.null(mito_data[[1]]$matrices$implied_mutCN)) {
  cat("Starting calculation of the implied mutation mtDNA copy number matrix",sep="\n")
  mito_data<-Map(list=mito_data,this_exp_ID=names(mito_data),function(list,this_exp_ID) {
    cat(this_exp_ID,sep="\n")
    this_exp_ID<-gsub("8pcw","8 pcw",this_exp_ID)
    mat<-list$matrices$vaf*list$matrices$SW
    sample_mtDNA_CNs<-sapply(colnames(mat),function(SampleID) {
      if(SampleID=="global") {
        return(0)
      } else {
        return(mito_cn%>%filter(Sample==SampleID)%>%slice_head(n=1)%>%pull(bedtools_mtDNA_genomes))
      }
    })
    mutCN_list<-lapply(1:nrow(mat),function(i) {
      mut_res<-mat[i,] * sample_mtDNA_CNs
      return(mut_res)
    })
    mutCN_mat<-mutCN_list%>%bind_rows()%>%as.matrix()
    list$matrices$implied_mutCN<-mutCN_mat
    return(list)
  })
  saveRDS(mito_data,file=mito_data_file)
}

##------------TEST MUTATIONS FOR INVERSE CORRELATION BETWEEN VAF AND COPY NUMBER -----------
#Mutations with this relationship are likely to be artefacts caused either by contamination or NUMTs

df_tidy_full<-dplyr::bind_rows(Map(list=mito_data,exp_ID=names(mito_data),function(list,exp_ID) {
  if(is.null(list)){stop(return(NULL))}
  df_tidy<-(list$matrices$vaf*list$matrices$SW)%>%
    tibble::rownames_to_column(var="mut_ref")%>%
    mutate(rho_vals=list$rho_vals)%>%
    dplyr::select(-global)%>%
    tidyr::gather(key="Sample",value="vaf",-mut_ref,-rho_vals)%>%
    #left_join(list$matrices$ML_Sig%>%as.data.frame()%>%tibble::rownames_to_column(var="mut_ref")%>%dplyr::select(-global)%>%tidyr::gather(key="Sample",value="ML_Sig",-mut_ref))%>%
    dplyr::filter(!grepl("DEL|INS",mut_ref))%>%
    dplyr::filter(vaf>0)%>%
    mutate(exp_ID=exp_ID)
  return(df_tidy)
}))

test_muts<-df_tidy_full%>%
  mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
  left_join(mito_cn,by="Sample")%>%
  filter(!is.na(bedtools_mtDNA_genomes))%>%
  pull(mut_ref)%>%
  unique()
length(test_muts)

pval_df<-lapply(test_muts,function(test_mut) {
  print(test_mut)
  if(which(test_muts==test_mut)%%1000==0){print(which(test_muts==test_mut))}
  temp=df_tidy_full%>%
    filter(mut_ref==test_mut)%>%
    mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
    left_join(mito_cn%>%dplyr::filter(!duplicated(Sample)))%>%
    mutate(implied_mut_CN=vaf*bedtools_mtDNA_genomes)%>%
    dplyr::select(Sample,
                  exp_ID,
                  mut_ref,
                  vaf,
                  bedtools_mtDNA_genomes,
                  #ML_Sig,
                  implied_mut_CN)
  if(nrow(temp%>%filter(vaf!=0&!is.na(bedtools_mtDNA_genomes)))<=2){
    return(data.frame(mut_ref=test_mut,coef=NA,pval=NA,max_vaf=NA))
  } else {
    max_prop=max(table(temp$exp_ID))/length(temp$exp_ID)
    if(max_prop>0.95) {
      which_max_prop=names(table(temp$exp_ID))[which.max(table(temp$exp_ID))]
      temp<-temp%>%filter(exp_ID==which_max_prop)
      assessment_type="single_sample"
      sample_assessed=which_max_prop
    } else {
      temp<-temp%>%filter(implied_mut_CN<=25)
      assessment_type="multi_sample"
      sample_assessed=which_max_prop=names(table(temp$exp_ID))[which.max(table(temp$exp_ID))]
      if(nrow(temp%>%filter(vaf!=0&!is.na(bedtools_mtDNA_genomes)))<=2){
        return(data.frame(mut_ref=test_mut,coef=NA,pval=NA,max_vaf=NA))
      }
    }
    
    #peak_sig=names(table(temp$ML_Sig))[which.max(table(temp$ML_Sig))]
    
    #Do linear regression for relationship of 1/sample mitochondrial copy number against sample VAF
    lm.summary<-summary(lm(1/bedtools_mtDNA_genomes~vaf,data=temp%>%filter(vaf!=0&!is.na(bedtools_mtDNA_genomes))))
    vaf_coefs<-lm.summary$coefficients['vaf',]
    
    #Return a summary data frame of the linear model results and some other statistics
    return(data.frame(mut_ref=test_mut,
                      #peak_sig=peak_sig,
                      sum_of_vaf=sum(temp$vaf),
                      assessment_type=assessment_type,
                      sample_assessed=sample_assessed,coef=vaf_coefs[1],
                      pval=vaf_coefs[4],
                      r2=lm.summary$r.squared,
                      max_vaf=temp%>%dplyr::slice_max(vaf,n=1)%>%pull(vaf)))
  }
})%>%dplyr::bind_rows()

#Remove values with negative coefficients or NA/NaN p-values - then perform multiple hypothesis testing correction (FDR method)
pval_df_filt<-pval_df%>%dplyr::filter(!is.na(pval)&!is.nan(pval)&coef>0)
pval_df_filt$qval<-p.adjust(pval_df_filt$pval,method="BH")
qval_cutoff=1e-2
sum(pval_df_filt$qval<qval_cutoff) #Number of correlating mutations (FDR set to 1%)
CN_correlating_muts<-pval_df_filt%>%dplyr::filter(qval<qval_cutoff)%>%pull(mut_ref)%>%sort()
saveRDS(CN_correlating_muts,file=paste0(root_dir,"data/CN_correlation.RDS"))

#Visualize this correlation for the top 20 strongest correlation mutations
Individual_cols=RColorBrewer::brewer.pal(12,"Paired")
names(Individual_cols)<-ref_df$Sample[order(ref_df$Age)]
CN_correlation_plot_20<-df_tidy_full%>%
  filter(mut_ref%in%(pval_df_filt%>%dplyr::filter(qval<1e-2)%>%slice_min(order_by = qval,n=20)%>%pull(mut_ref)%>%sort()))%>%
  mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
  left_join(mito_cn)%>%
  mutate(implied_mut_CN=vaf*bedtools_mtDNA_genomes)%>%
  dplyr::select(exp_ID,Sample,mut_ref,vaf,bedtools_mtDNA_genomes,implied_mut_CN)%>%
  filter(implied_mut_CN<8 & vaf!=0)%>%
  mutate(mut_ref=gsub("MT_","",mut_ref))%>%
  ggplot(aes(x=vaf, y=1/bedtools_mtDNA_genomes,col=exp_ID))+
  geom_point(alpha=0.25,size=0.2)+
  geom_smooth(aes(x=vaf,y=1/bedtools_mtDNA_genomes),col="black",size=0.5,method="lm",inherit.aes = F)+
  facet_wrap(~mut_ref,ncol = 10)+
  scale_color_manual(values=Individual_cols)+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  my_theme+
  theme(strip.text.x = element_text(size=5,margin = unit(c(1,0,1,0),"mm")),
        axis.text.x=element_text(angle=90),
        legend.key.height = unit(3,"mm"))+
  labs(x="Mutation VAF",y="1/mtDNA genomes",col="Individual")+
  guides(colour = guide_legend(override.aes = list(size=2,alpha=1)))
ggsave(filename = paste0(figures_dir,"Supp_Figure_01/CN_correlation_plot_20.pdf"),plot=CN_correlation_plot_20,height=3,width=7)
