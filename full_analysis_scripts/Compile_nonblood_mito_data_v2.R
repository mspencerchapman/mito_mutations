library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(phylosignal)
options(stringsAsFactors = F)

#Source functions needed for the script
my_working_directory<-getwd()
R_function_files = list.files("/lustre/scratch126/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)
genomeFile="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"
root_dir="/lustre/scratch125/casm/team268im/al28/mtDNA"

#Import the mitochondrial copy number data
mito_cn=read.csv("/lustre/scratch125/casm/team268im/al28/mtDNA/whole_genome_coverage_pileup_and_bedtools_annotated.csv",header=T)

# #Set these file paths before running the script
# genomeFile="~/Documents/Reference_files/hs37d5.fa"
# root_dir="~/R_work/mito_mutations_blood/"
# source(paste0(root_dir,"/data/mito_mutations_blood_functions.R"))

translate=data.frame(dataset=c("KY","HL","SO","PR","LM","NW","lymph","CML"),
                     al_ref=c("lung_organoid","colon","colon_ibd","muty_mutant","endometrium","blood_MPN","immune","blood_CML"))

muty_samples=c("PD44887","PD44888","PD44889","PD44890","PD44891")

#Import metadata relating to all indviduals studied
ref_df<-readxl::read_excel("/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/non_blood_metadata.xlsx")

for(dataset in translate$dataset) {
  
  al_ref=translate$al_ref[translate$dataset==dataset]
  
  if(dataset=="lymph") {
    lymph_metadata=read.csv("/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Lymph_metadata.csv",stringsAsFactors = F)
    dataset_Donors=unique(lymph_metadata$Donor)
    lymph_samples_list=lapply(dataset_Donors,function(this_Donor) {lymph_metadata%>%filter(Donor==this_Donor)%>%pull(colony)})
    names(lymph_samples_list)<-dataset_Donors
    
    ref_file=paste0("/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Samples_metadata_ref.csv")
    ref_df<-read.csv(ref_file)
    ref_df$PDID[ref_df$Sample=="AX001"]<-"TG01"
  } else {
    samples_file=paste0("/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"_samples.txt")
    trees_folder=paste0("/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data/",dataset)
    dataset_samples=readLines(samples_file)
    
    #Set the key file paths using the root dir
    dataset_tree_file_paths = list.files(trees_folder,pattern=".tree",full.names = T,recursive=T)
    if(dataset=="PR") {
      dataset_samples=muty_samples
      dataset_tree_file_paths=grep("polytom",dataset_tree_file_paths,value=T)
    }
  }

  MT_count_file_paths = list.files(paste0("/lustre/scratch125/casm/team268im/al28/mtDNA/Data/",al_ref),pattern = "MT_count.csv")
  
  ###---------IMPORT THE DATA - PILEUP DATA AND SHEARWATER DATA INTO A MATRIX---------########
  if(dataset=="lymph") {
    dataset_Donors <-c("AX001","KX001","KX002","KX003","TX001","TX002")
    dataset_mito_data<-lapply(dataset_Donors,function(exp_ID) {
      print(exp_ID)
      PD_ID=ref_df$PDID[ref_df$Sample==exp_ID]
      
      #This is a version of the function that doesn't need a tree
      res<-generate_mito_matrices_mod(PD_number = PD_ID,
                                      sample_vector = lymph_samples_list[[exp_ID]],
                                      shearwater_calls_file=paste0("/lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/shearwater_stringent_exclusion/",al_ref,"/",PD_ID,"_shearwater_mutations.tsv"),
                                      pileup_folder=paste0("/lustre/scratch125/casm/team268im/al28/mtDNA/Data/",al_ref),
                                      run_bb = T,
                                      exclude_samples = NULL)
      return(res)
    })
    names(dataset_mito_data)<-dataset_Donors
  } else {
    dataset_mito_data<-lapply(dataset_samples,function(exp_ID) {
      print(exp_ID)
      if(dataset=="SO") {
        tree_file_path=grep(paste0(exp_ID,"_"),dataset_tree_file_paths,value = T)
      } else if(dataset=="CML") {
        tree_file_path=grep(paste0("moltime_",exp_ID,".tree"),dataset_tree_file_paths,value = T)
      } else {
        tree_file_path=grep(exp_ID,dataset_tree_file_paths,value = T)
      }
      res<-generate_mito_matrices(PD_number = exp_ID,
                                  tree_file_path = tree_file_path,
                                  shearwater_calls_file=paste0("/lustre/scratch125/casm/team268im/al28/mtDNA/Analysis/shearwater_stringent_exclusion/",al_ref,"/",exp_ID,"_shearwater_mutations.tsv"),
                                  pileup_folder=paste0("/lustre/scratch125/casm/team268im/al28/mtDNA/Data/",al_ref),
                                  run_bb = T,
                                  rev_germline=F,
                                  exclude_samples = NULL)
      return(res)
    })
    names(dataset_mito_data)<-dataset_samples
    
    if(any(sapply(dataset_mito_data, is.null))) {
      dataset_mito_data<-dataset_mito_data[-which(sapply(dataset_mito_data, is.null))]
    }
    
    #Add the ultrametric trees to the data lists
    if(dataset=="CML") {
      dataset_mito_data<-Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID) {
        ultra_tree_file_path=grep(paste0("tree_",exp_ID,".tree"),dataset_tree_file_paths,value = T)
        list$tree.ultra<-ape::read.tree(ultra_tree_file_path)
        list$tree.ultra<-plot_tree(list$tree.ultra,b_do_not_plot=T)
        return(list)
      })
    } else {
      dataset_mito_data<-lapply(dataset_mito_data,function(list) {
        list$tree<-ape::di2multi(list$tree)
        list$tree$edge.length[which(list$tree$edge.length==0)]<-1
        tree.ultra<-make.ultrametric.tree(list$tree)
        tree.ultra<-add_ancestral_outgroup(tree.ultra)
        tree.ultra$coords<-NULL
        tree.ultra$edge.length<-tree.ultra$edge.length*mean(get_mut_burden(list$tree))
        list$tree.ultra<-tree.ultra
        list$tree.ultra<-plot_tree(list$tree.ultra,b_do_not_plot=T)
        return(list)
      })
    }
    
  }
  
  #Remove samples from the tree that were excluded from shearwater due to contamination OR are not in the mitochondrial CN dataframe (this should include all sapmles run through shearwater)
  dataset_mito_data<-lapply(dataset_mito_data,function(list) {
    
    #Get list of samples that are both included in the shearwater calls & mitochondrial CN data
    nonexcluded_samples<-unique(list$sample_shearwater_calls$sampleID)
    nonexcluded_samples<-nonexcluded_samples[nonexcluded_samples%in%mito_cn$Sample]
    
    #Update the matrices to include only these samples
    list$matrices<-lapply(list$matrices,function(mat) {
      return(mat[,which(colnames(mat)%in%nonexcluded_samples)])
    })
    
    #Update the trees to include only these samples (if present)
    if(!is.null(list$tree)) {list$tree<-keep.tip(phy=list$tree,tip = list$tree$tip.label[which(list$tree$tip.label%in%nonexcluded_samples)])}
    if(!is.null(list$tree.ultra)) {list$tree.ultra<-keep.tip(phy=list$tree.ultra,tip = list$tree.ultra$tip.label[which(list$tree.ultra$tip.label%in%nonexcluded_samples)])}
    
    #Return the updated list
    return(list)
  })
  
  #This is the calculated absolute number of mutated mtDNA copies per cell, using the VAF and mtDNA copy number
  if(is.null(dataset_mito_data[[1]]$matrices$implied_mutCN)) {
    cat("Starting calculation of the implied mutation mtDNA copy number matrix",sep="\n")
    dataset_mito_data<-Map(list=dataset_mito_data,this_exp_ID=names(dataset_mito_data),function(list,this_exp_ID) {
      cat(this_exp_ID,sep="\n")
      this_exp_ID<-gsub("8pcw","8 pcw",this_exp_ID)
      mat<-list$matrices$vaf
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
    saveRDS(dataset_mito_data,file=paste0("mito_mutation_data_",dataset,".RDS"))
  }
  
  df_tidy_full<-dplyr::bind_rows(Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID) {
    if(is.null(list)){stop(return(NULL))}
    df_tidy<-(list$matrices$vaf)%>%
      tibble::rownames_to_column(var="mut_ref")%>%
      mutate(rho_vals=list$rho_vals)%>%
      dplyr::select_if(!names(.) %in% "global")%>%
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
    #print(test_mut)
    if(which(test_muts==test_mut)%%1000==0){print(which(test_muts==test_mut))}
    temp=df_tidy_full%>%
      filter(mut_ref==test_mut)%>%
      mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
      left_join(mito_cn%>%dplyr::filter(!duplicated(Sample)),by="Sample")%>%
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
  saveRDS(CN_correlating_muts,file=paste0("CN_correlation_",dataset,".RDS"))
}

##PULL OUT THE SHARED MUTATIONS - embed as an additional object within the mito_data list
#Find those mitochondrial mutations that are present in more than one sample at the specified cut_off
#Calculate the "phylogenetic signal" of each of these mutations i.e. the degree to which they follow the phylogeny

##CHOOSE DATASET AND START ANALYSIS
all_cohorts=c("CML","KY","HL","SO","PR","LM","NW","lymph","blood")
rho_cut_off=0
exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_456_C_T","MT_567_A_C","MT_574_A_C","MT_8270_C_T","MT_16170_A_C","MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C") #These are the ones to exclude for the endometrial analysis
all_mito_datasets<-lapply(all_cohorts,function(dataset) {
  if(dataset=="blood") {
    dataset_mito_data<-readRDS("../mito_data.Rds")
    CN_correlating_muts<-readRDS("../CN_correlation.RDS")
  } else {
    CN_correlating_muts=readRDS(paste0("CN_correlation_",dataset,".RDS"))
    dataset_mito_data=readRDS(paste0("mito_mutation_data_",dataset,".RDS"))
  }
  
  #Remove samples that aren't in the shearwater calls dataframe - these have failed the 'haplocheck' or 'verifyBamID' tests
  if(dataset=="lymph") {
    dataset_mito_data<-lapply(dataset_mito_data, function(list) {
      list$tree<-list$tree.ultra<-ape::stree(n = length(unique(list$sample_shearwater_calls$sampleID)),tip.label = unique(list$sample_shearwater_calls$sampleID))
      list$tree$edge.length<-rep(10,nrow(list$tree$edge))
      list$tree.ultra$edge.length<-rep(1,nrow(list$tree$edge))
      return(list)
    })
  }

  #Deal with the germline mutations (if not done already)
  dataset_mito_data<-Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID) {
    if(is.null(list$germline_muts)) {
      list$germline_muts<-identify_germline(matrices=list$matrices,threshold=0.9)
      cat(list$germline_muts,sep="\n")
      list$matrices<-reverse_germline(matrices=list$matrices,threshold=0.9)
    }
    return(list)
  })
  
  #Sort out the names of the mutCN column
  dataset_mito_data<-lapply(dataset_mito_data,function(list) {
    colnames(list$matrices$implied_mutCN)<-stringr::str_split(colnames(list$matrices$implied_mutCN),pattern = '\\.\\.\\.',simplify=T)[,1]
    return(list)
  })
  
  #Add the CN correlating muts vector to each patient (slightly inefficient but programatically more straight-forward)
  dataset_mito_data<-lapply(dataset_mito_data,function(list) {
    list$CN_correlating_muts<-CN_correlating_muts
    return(list)
  })
  
  cat("Adding the filtered VAF matrix.",sep="\n")
  #Annotate specific mutations with their most likely signature using the 'sig_ref' dataframe
  dataset_mito_data<-Map(list=dataset_mito_data,this_exp_ID=names(dataset_mito_data),function(list,this_exp_ID) {
    cat(this_exp_ID,sep="\n")
    
    if(dataset=="lymph") {list$tree$tip.label<-unique(list$sample_shearwater_calls$sampleID)}
    
    mutCN_cutoff=25 #If the mitochondrial copy number is over 25, retain mutation even if is in the "CN correlating muts" list
    CN_correlating_mut_removal_mat=list$matrices$implied_mutCN>mutCN_cutoff|(matrix((!rownames(list$matrices$vaf)%in%CN_correlating_muts),ncol=1)%*%matrix(rep(1,ncol(list$matrices$vaf)),nrow=1))
    vaf.filt<-(list$matrices$vaf[,list$tree$tip.label]*list$matrices$SW[,list$tree$tip.label][,list$tree$tip.label]*CN_correlating_mut_removal_mat[,list$tree$tip.label])[!grepl("DEL|INS",rownames(list$matrices$vaf))&
                                                                                                                                                                            list$rho_vals>rho_cut_off&
                                                                                                                                                                            !rownames(list$matrices$vaf)%in%exclude_muts,]
    
    list$matrices$CN_correlating_mut_removal_mat<-CN_correlating_mut_removal_mat
    list$matrices$vaf.filt<-vaf.filt
    return(list)
  })
  
  if(dataset=="lymph") {
    dataset_mito_data<-Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID) {
      if(!exp_ID=="TX001") {list$het_oocyte_muts<-c()}
      return(list)
    })
  }
  return(dataset_mito_data)
})
names(all_mito_datasets)<-all_cohorts

temp=Map(dataset_mito_data=all_mito_datasets,dataset=all_cohorts,function(dataset_mito_data,dataset) {
  cat(dataset,sep="\n")

  CN_correlating_muts=dataset_mito_data[[1]]$CN_correlating_muts
  figures_dir=paste0(dataset,"/")
  
  dataset_mito_data<-Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data), function(list,exp_ID) {
    cat(exp_ID,sep="\n")
    mut_vaf_cutoff=0.03
    
    vaf.filt<-list$matrices$vaf.filt
    tree.ultra<-list$tree.ultra
    
    ##This doesn't work if the trees are too small
    if(length(tree.ultra$tip.label)<=5) {
      list$shared_muts_df<-NULL
      stop(return(list))
    }
    
    shared_muts<-rownames(vaf.filt)[rowSums(vaf.filt>mut_vaf_cutoff,na.rm = T)>1]
    
    if(exp_ID=="PD37590") {
      excl_muts<-c("MT_955_A_C","MT_966_A_C","MT_953_T_C","MT_961_T_C")
      shared_muts<-shared_muts[!shared_muts%in%excl_muts]
    }
    n_pos<-rowSums(vaf.filt[shared_muts,]>mut_vaf_cutoff,na.rm = T)
    print(paste("There are",length(shared_muts),"shared mutations"))
    
    if(length(shared_muts)==0) {
      list$shared_muts_df<-NULL
      stop(return(list))
    }
    
    #Test shared muts with phylosignal
    tree4d<-phylobase::phylo4d(drop.tip(tree.ultra,"Ancestral"),tip.data=t(vaf.filt[shared_muts,list$tree$tip.label]))
    
    if(!dataset=="lymph") {res_cor=phyloSignal(tree4d)} else {res_cor=list(stat=matrix(NA,nrow=length(shared_muts),ncol=2,dimnames = list(shared_muts,c("Lambda","Cmean"))),pval=matrix(NA,nrow=length(shared_muts),ncol=2,dimnames = list(shared_muts,c("Lambda","Cmean"))))}
    
    #Add phylosignal info onto the dataframe
    shared_muts_df<-data.frame(mut=shared_muts,n_pos=n_pos)%>%
      tidyr::separate(mut,into=c("Chrom","Pos","Ref","Alt"),sep="_",remove=F)%>%
      mutate(Pos=as.numeric(Pos))%>%
      mutate(global_VAF=list$matrices$vaf[mut,"global"])%>%
      mutate(lambda_pval=res_cor$pval[mut,"Lambda"],Cmean_pval=res_cor$pval[mut,"Cmean"])%>%
      mutate(lambda=res_cor$stat[mut,"Lambda"],Cmean=res_cor$stat[mut,"Cmean"])%>%
      arrange(Pos)
    list$shared_muts_df<-shared_muts_df
    return(list)
  })
  
  #Visualize the shared mutations
  col_scheme<-c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9])(100))
  lowest_VAF_to_show=0.03
  tree_to_show="tree.ultra"
  temp=Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),f=function(list,exp_ID){
    cat(exp_ID,sep="\n"); Age<-ref_df%>%filter(ID==exp_ID)%>%pull(Age)
    if(is.null(list)) {stop(return(NULL))}
    
    #Don't bother printing if the tree is <10 samples
    if(length(list[[tree_to_show]]$tip.label)<8) {stop(return(NULL))}
    
    if(is.null(list$shared_muts_df)) {
      pdf(file = paste0(figures_dir,exp_ID,"_shared_muts.pdf"),width = 7,height=4)
      plot_tree(tree = list[[tree_to_show]],cex.label = 0,plot_axis=F,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    if(nrow(list$shared_muts_df)==0) {
      pdf(file = paste0(figures_dir,exp_ID,"_shared_muts.pdf"),width = 7,height=4)
      plot_tree(tree = list[[tree_to_show]],cex.label = 0,plot_axis=F,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    bb_df<-data.frame(mut=rownames(list$matrices$vaf),rho=list$rho_vals)
    
    vaf.mtx<-list$matrices$vaf.filt*list$matrices$SW[rownames(list$matrices$vaf.filt),]
    
    plot_muts<-list$shared_muts_df%>%
      filter(mut%in%rownames(vaf.mtx))%>%
      pull(mut)
    
    #Remove samples if not included at all in the shearwater table - these samples are excluded due to contamination
    list[[tree_to_show]]<-keep.tip(list[[tree_to_show]],list[[tree_to_show]]$tip.label[list[[tree_to_show]]$tip.label%in%list$sample_shearwater_calls$sampleID])
    list[[tree_to_show]]$coords<-NULL
    vaf.mtx<-vaf.mtx[,list[[tree_to_show]]$tip.label]
    
    names(col_scheme)<-seq(0,1,0.01)
    hm<-matrix(0,nrow=length(plot_muts),ncol=length(list[[tree_to_show]]$tip.label),dimnames = list(plot_muts,list[[tree_to_show]]$tip.label))
    for(i in 1:length(plot_muts)) {
      mut<-plot_muts[i]
      mut_vafs<-vaf.mtx[mut,list[[tree_to_show]]$tip.label]
      mut_vafs[mut_vafs<lowest_VAF_to_show]<-0
      mut_vafs<-round(mut_vafs,digits=2)
      hm[i,]<-col_scheme[as.character(mut_vafs)]
    }
    if(length(plot_muts)>=2) {
      plot_muts.clustered<-hclust(dist(vaf.mtx[plot_muts,]))
    } else {
      plot_muts.clustered=list(order=1)
    }
    
    par(mfrow=c(1,1))
    pdf(file = paste0(figures_dir,exp_ID,"_shared_muts.pdf"),width = 7,height=4)
    list[[tree_to_show]]=plot_tree(tree = list[[tree_to_show]],cex.label = 0,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
    add_mito_mut_heatmap(tree=list[[tree_to_show]],heatmap=hm[plot_muts.clustered$order,,drop=F],border="gray",heatmap_bar_height=0.1,cex.label = 0.25)
    dev.off()
  })
  
  
  tree_to_show="tree"
  col_scheme<-c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9])(100))
  lowest_VAF_to_show=0.03
  temp=Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),f=function(list,exp_ID){
    cat(exp_ID,sep="\n");Age<-ref_df%>%filter(ID==exp_ID)%>%pull(Age)
    if(is.null(list)) {stop(return(NULL))}
    
    #Don't bother printing if the tree is <8 samples
    if(length(list$tree$tip.label)<8) {stop(return(NULL))}
    
    #Plot the trees even if there are no shared muts
    if(is.null(list$shared_muts_df)) {
      pdf(file = paste0(figures_dir,exp_ID,"_raw_shared_muts.pdf"),width = 7,height=4)
      list[[tree_to_show]]=plot_tree(tree = list[[tree_to_show]],cex.label = 0,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    if(nrow(list$shared_muts_df)==0) {
      pdf(file = paste0(figures_dir,exp_ID,"_raw_shared_muts.pdf"),width = 7,height=4)
      list[[tree_to_show]]=plot_tree(tree = list[[tree_to_show]],cex.label = 0,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    bb_df<-data.frame(mut=rownames(list$matrices$vaf),rho=list$rho_vals)
    plot_muts<-list$shared_muts_df%>%
      #filter(Cmean_pval<0.05)%>%
      pull(mut)
    
    if(length(plot_muts)==0) {
      pdf(file = paste0(figures_dir,exp_ID,"_raw_shared_muts.pdf"),width = 7,height=4)
      list[[tree_to_show]]=plot_tree(tree = list[[tree_to_show]],cex.label = 0,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    
    vaf.mtx<-list$matrices$vaf[rownames(list$matrices$vaf.filt),]
    
    names(col_scheme)<-seq(0,1,0.01)
    hm<-matrix(0,nrow=length(plot_muts),ncol=length(list[[tree_to_show]]$tip.label),dimnames = list(plot_muts,list[[tree_to_show]]$tip.label))
    for(i in 1:length(plot_muts)) {
      mut<-plot_muts[i]
      mut_vafs<-vaf.mtx[mut,list[[tree_to_show]]$tip.label]
      mut_vafs[mut_vafs<lowest_VAF_to_show]<-0
      mut_vafs<-round(mut_vafs,digits=2)
      hm[i,]<-col_scheme[as.character(mut_vafs)]
    }
    if(length(plot_muts)>=2) {
      plot_muts.clustered<-hclust(dist((vaf.mtx^(1/4))[plot_muts,]))
    } else {
      plot_muts.clustered=list(order=1)
    }
    
    par(mfrow=c(1,1))
    #fig<-ifelse(exp_ID=="KX004","Figure_05/","Supp_Figure_05/")
    pdf(file = paste0(figures_dir,exp_ID,"_raw_shared_muts.pdf"),width = 7,height=4)
    list[[tree_to_show]]=plot_tree(tree = list[[tree_to_show]],cex.label = 0,cex.axis=0.5,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
    add_mito_mut_heatmap(tree=list[[tree_to_show]],heatmap=hm[plot_muts.clustered$order,,drop=F],border="gray",heatmap_bar_height=0.1,cex.label = 0.25)
    dev.off()
  })
  
  #PLOT ALL MUTS THAT PASS SHEARWATER IN AT LEAST ONE SAMPLE - BIG HEATMAPS
  #Higher granularity visualization - show down to 0.001 VAF
  lowest_VAF_to_show=0.001
  col_scheme_high<-c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9])(1001))
  names(col_scheme_high)<-seq(0,1,0.001)
  temp=Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),f=function(list,exp_ID){
    cat(exp_ID,sep="\n"); Age<-ref_df%>%filter(ID==exp_ID)%>%pull(Age)
    if(is.null(list)) {stop(return(NULL))}
    bb_df<-data.frame(mut=rownames(list$matrices$vaf),rho=list$rho_vals)
    plot_muts<-rownames(list$matrices$vaf.filt)
    
    if(length(plot_muts)==0) {
      pdf(file = paste0(figures_dir,exp_ID,"_raw_shared_muts.pdf"),width = 7,height=4)
      list$tree=plot_tree(tree = list$tree,cex.label = 0,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : All mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    
    #Showing the raw VAF here - not filtered by SW
    vaf.mtx<-list$matrices$vaf.filt[rownames(list$matrices$vaf.filt),]
    
    hm<-matrix(0,nrow=length(plot_muts),ncol=length(list$tree$tip.label),dimnames = list(plot_muts,list$tree$tip.label))
    for(i in 1:length(plot_muts)) {
      mut<-plot_muts[i]
      mut_vafs<-vaf.mtx[mut,list$tree$tip.label]
      mut_vafs[mut_vafs<lowest_VAF_to_show]<-0
      mut_vafs<-round(mut_vafs,digits=3)
      hm[i,]<-col_scheme_high[as.character(mut_vafs)]
    }
    if(length(plot_muts)>=2) {
      plot_muts.clustered<-hclust(dist((vaf.mtx^(1/4))[plot_muts,]))
    } else {
      plot_muts.clustered=list(order=1)
    }
    
    par(mfrow=c(1,1))
    pdf(file = paste0(figures_dir,exp_ID,"_all_muts_filt.pdf"),width = 7,height=8)
    list$tree=plot_tree(tree = list$tree,cex.label = 0,plot_axis=T,vspace.reserve = 4,title=paste0(exp_ID," (",Age," years old) : All mutations (VAF>",100*lowest_VAF_to_show,"%)"))
    add_mito_mut_heatmap(tree=list$tree,heatmap=hm[plot_muts.clustered$order,,drop=F],border="gray",heatmap_bar_height=0.01,cex.label = 0.25)
    dev.off()
  })
  
  #PLOT SUSPECTED GERMLINE MUTATIONS - CHECK THEM
  #Higher granularity visualization - show down to 0.001 VAF
  col_scheme<-c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9])(100))
  lowest_VAF_to_show=0.03
  tree_to_show="tree.ultra"
  temp=Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),f=function(list,exp_ID){
    cat(exp_ID,sep="\n"); Age<-ref_df%>%filter(ID==exp_ID)%>%pull(Age)
    if(is.null(list)) {stop(return(NULL))}
    bb_df<-data.frame(mut=rownames(list$matrices$vaf),rho=list$rho_vals)
    plot_muts<-grep(pattern = "\\*",rownames(list$matrices$vaf),value=T)
    
    if(length(plot_muts)==0) {
      pdf(file = paste0(figures_dir,exp_ID,"_germline_muts.pdf"),width = 7,height=4)
      list$tree=plot_tree(tree = list$tree,cex.label = 0,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : germline mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    
    #Showing the raw VAF here - not filtered by SW
    vaf.mtx<-list$matrices$vaf
    
    hm<-matrix(0,nrow=length(plot_muts),ncol=length(list$tree$tip.label),dimnames = list(plot_muts,list$tree$tip.label))
    for(i in 1:length(plot_muts)) {
      mut<-plot_muts[i]
      mut_vafs<-vaf.mtx[mut,list$tree$tip.label]
      mut_vafs[mut_vafs<lowest_VAF_to_show]<-0
      mut_vafs<-round(mut_vafs,digits=3)
      hm[i,]<-col_scheme_high[as.character(mut_vafs)]
    }
    if(length(plot_muts)>=2) {
      plot_muts.clustered<-hclust(dist((vaf.mtx^(1/4))[plot_muts,]))
    } else {
      plot_muts.clustered=list(order=1)
    }
    
    par(mfrow=c(1,1))
    pdf(file = paste0(figures_dir,exp_ID,"_germline_muts.pdf"),width = 7,height=8)
    list$tree=plot_tree(tree = list$tree,cex.label = 0,plot_axis=T,vspace.reserve = 4,title=paste0(exp_ID," (",Age," years old) : germline mutations (VAF>",100*lowest_VAF_to_show,"%)"))
    add_mito_mut_heatmap(tree=list$tree,heatmap=hm[plot_muts.clustered$order,,drop=F],border="gray",heatmap_bar_height=0.01,cex.label = 0.25)
    dev.off()
  })
})

###--------MEASURE 'CLONAL MARKING' OF EXPANDED CLADES BY MITOCHONDRIAL MUTATIONS--------------------
#Define function used to recognise the clonal expansions
get_expanded_clade_nodes=function(tree,height_cut_off=100,min_clonal_fraction=0.02,min_samples=1){
  nodeheights=nodeHeights(tree)
  
  #This pulls out nodes that fulfill on the criteria: branches cross the cut-off & contain the minimum proportion of samples
  nodes=tree$edge[,2][nodeheights[,1] < height_cut_off &
                        !nodeheights[,2] < height_cut_off &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))/length(tree$tip.label)})>min_clonal_fraction &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))})>=min_samples]
  df=data.frame(nodes=nodes,n_samples=sapply(nodes,function(node) {length(getTips(tree,node))}),MRCA_time=sapply(nodes,function(node) {nodeheight(tree,node)}),clonal_fraction=sapply(nodes,function(node) {length(getTips(tree,node))/length(tree$tip.label)}))
  return(df)
}

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 7),
                axis.title = element_text(size=8),
                axis.line = element_line(linewidth = 0.4),
                axis.ticks = element_line(linewidth = 0.3),
                legend.text = element_text(size=6),
                legend.title = element_text(size=8),
                strip.text = element_text(size=7),
                strip.background = element_rect(fill="lightgray",linewidth = 0.4),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))


expanded.clades.plot<-dplyr::bind_rows(Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID){
  exp_nodes<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off = ifelse(dataset=="NW",50,ifelse(dataset=="CML",2,200)),min_clonal_fraction=0.01,min_samples = 2)
  if(nrow(exp_nodes)==0) {stop(return(NULL))}
  exp_nodes$exp_ID<-exp_ID
  return(exp_nodes)
}))%>%
  arrange(desc(n_samples))%>%
  ggplot(aes(x=exp_ID,y=clonal_fraction,fill=n_samples))+
  geom_bar(stat="identity",position="stack",col="black")+
  scale_fill_gradient(low="lightgrey",high="darkred")+
  labs(x="Individual",y="Clonal fraction",fill="Number of\n samples\n in clone")+
  theme_bw()+
  my_theme+
  theme(axis.text.x=element_text(angle=90))
ggsave(filename=paste0(figures_dir,dataset,"_expanded_clades_plot.pdf"),expanded.clades.plot,width=4,height=2.5)

expanded_clades_df<-Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID){
  cat(paste0(exp_ID,"\n"))
  
  vaf.filt<-list$matrices$vaf.filt
  
  #Review how many expanded clades have reliable mitochondrial marker mutations
  marker_mut_cutoff<-0.03
  pos_mut_cutoff<-0.03
  exp_nodes<-get_expanded_clade_nodes(list$tree,height_cut_off = ifelse(dataset=="NW",50,200),min_clonal_fraction=0.01,min_samples = 2)
  if(nrow(exp_nodes)==0|is.null(list$shared_muts_df)) {stop(return(NULL))}
  full_df<-dplyr::bind_cols(data.frame(exp_ID=rep(exp_ID,nrow(exp_nodes))),
                            exp_nodes,
                            data.frame(marker_mut_cutoff=rep(marker_mut_cutoff,nrow(exp_nodes))),
                            exp_nodes_muts<-dplyr::bind_rows(lapply(exp_nodes$nodes,function(node) {
                              cat(node,sep="\n")
                              node_samples=getTips(list$tree,node)
                              if(any(vaf.filt[,node_samples]>marker_mut_cutoff)){
                                
                                #Get mutations that are present in at least one sample in the clade above the 'marker_mut' VAF cutoff
                                node_homo_muts<-names(rowSums(vaf.filt[,node_samples,drop=F]>marker_mut_cutoff)[rowSums(vaf.filt[,node_samples,drop=F]>marker_mut_cutoff)>0])
                                
                                #Test samples outside the node - exclude if mutation is also present in a significant proportion of samples outside the node in question
                                #This will primarily exclude mutations that are high level oocyte heteroplasmies that are not informative
                                outside_node_samples<-list$tree$tip.label[!list$tree$tip.label%in%node_samples]
                                pos_outside_node_samples<-rowSums(vaf.filt[node_homo_muts,outside_node_samples,drop=F]>marker_mut_cutoff)
                                prop_pos_outside_node_samples<-pos_outside_node_samples/length(outside_node_samples)
                                node_homo_muts<-node_homo_muts[prop_pos_outside_node_samples<0.2]
                                
                                if(length(node_homo_muts)==0) {
                                  stop(return(data.frame(nmuts=0,homo_muts=NA,pos_samples_per_mut=NA,max_pos_samples=0,max_pos_prop=0,mean_heteroplasmy=NA)))
                                }
                                
                                pos_samples_per_mut<-sapply(node_homo_muts,function(mut){
                                  n_samples<-sum(vaf.filt[mut,node_samples,drop=F]>pos_mut_cutoff)
                                  return(n_samples)
                                })
                                
                                mean_het_of_pos<-sapply(node_homo_muts,function(mut){
                                  pos_samples<-node_samples[vaf.filt[mut,node_samples]>pos_mut_cutoff]
                                  return(mean(as.numeric(vaf.filt[mut,pos_samples])))
                                })
                                
                                return(data.frame(nmuts=length(node_homo_muts),
                                                  homo_muts=paste0(node_homo_muts,collapse=","),
                                                  BMM=node_homo_muts[which.max(pos_samples_per_mut)],
                                                  pos_samples_per_mut=paste0(pos_samples_per_mut,collapse=","),
                                                  max_pos_samples=max(pos_samples_per_mut),
                                                  max_pos_prop=max(pos_samples_per_mut)/length(node_samples),
                                                  mean_heteroplasmy=mean(as.numeric(vaf.filt[node_homo_muts[which.max(pos_samples_per_mut)],node_samples])),
                                                  mean_het_of_pos=mean_het_of_pos[which.max(pos_samples_per_mut)]))
                              } else {
                                return(data.frame(nmuts=0,homo_muts=NA,pos_samples_per_mut=NA,max_pos_samples=0,max_pos_prop=0,mean_heteroplasmy=NA))
                              }
                            })))
  #return(dplyr::select(full_df,-homo_muts,-pos_samples_per_mut))
  return(full_df)
})%>%dplyr::bind_rows()

node_factor_levels=expanded_clades_df%>%arrange(n_samples)%>%mutate(levels=stringr::str_c(exp_ID,nodes,sep = "_"))%>%pull(levels)
expanded.clades.marking.plot<-expanded_clades_df%>%
  mutate(n_neg_samples=n_samples-max_pos_samples)%>%
  dplyr::select(exp_ID,nodes,n_samples,max_pos_samples,n_neg_samples)%>%
  mutate(max_pos_samples=ifelse(max_pos_samples<=1,0,max_pos_samples))%>%
  mutate(n_neg_samples=n_samples-max_pos_samples)%>%
  tidyr::gather(-exp_ID,-nodes,-n_samples,key="Pos_or_neg",value="n_samples")%>%
  mutate(Pos_or_neg=ifelse(Pos_or_neg=="n_neg_samples","Absent","Present"))%>%
  mutate(Pos_or_neg=factor(Pos_or_neg,levels=c("Present","Absent")))%>%
  mutate(levels=str_c(exp_ID,nodes,sep = "_"))%>%
  ggplot(aes(x=factor(levels,levels=node_factor_levels),y=n_samples,fill=Pos_or_neg))+
  geom_bar(position="stack",stat="identity",col="black",linewidth=0.15)+
  scale_fill_brewer(palette="Set2")+
  facet_grid(~exp_ID,scales="free",space = "free")+
  #scale_y_continuous(breaks=seq(0,10,2))+
  theme_bw()+
  my_theme+
  theme(panel.spacing=unit(0.5,"mm"),axis.text.x = element_blank(),strip.text.x = element_text(angle=90))+
  labs(x="Clonal expansion",y=str_wrap("Number of samples within expansion",width=15),fill=str_wrap("Best mitochondrial marker mutation",width=10))
ggsave(filename=paste0("expanded_clades_marking_plot_",dataset,".pdf"),expanded.clades.marking.plot,width=7,height=2.3)

expanded_clades_df$n_marker<-sapply(1:nrow(expanded_clades_df),function(i) {
  pos_samples_per_mut<-as.numeric(stringr::str_split(expanded_clades_df$pos_samples_per_mut[i],pattern = ",")[[1]])
  n_perfect_markers<-sum(pos_samples_per_mut==expanded_clades_df$n_samples[i])
  return(n_perfect_markers)
})

expanded_clades_df%>%
  filter(!is.na(n_marker))%>%
  ggplot(aes(x=factor(n_marker),y=MRCA_time,col=exp_ID))+
  geom_jitter(height=0,width=0.2)+
  theme_classic()+
  my_theme