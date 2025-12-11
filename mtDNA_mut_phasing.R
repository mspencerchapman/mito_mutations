### IMPORT ALL THE DATA
library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(phylosignal)
options(stringsAsFactors = F)

#Set the basic plotting theme for ggplot2
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

#Source functions needed for the script
my_working_directory<-"/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood"
R_function_files = list.files("/lustre/scratch126/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)
genomeFile="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"
root_dir="/lustre/scratch125/casm/team268im/al28/mtDNA"

#Import the mitochondrial copy number data
mito_cn=read.csv("/lustre/scratch125/casm/team268im/al28/mtDNA/whole_genome_coverage_pileup_and_bedtools_annotated.csv",header=T)

translate=data.frame(dataset=c("KY","HL","SO","PR","LM","NW","lymph"),
                     al_ref=c("lung_organoid","colon","colon_ibd","muty_mutant","endometrium","blood_MPN","immune"))

muty_samples=c("PD44887","PD44888","PD44889","PD44890","PD44891")

#Import metadata relating to all indviduals studied
ref_df<-readxl::read_excel("/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/non_blood_metadata.xlsx")

#Import colony-level data for the lymphoid dataset as there are multiple different cell types
colony_info<-read.delim("colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt")
exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_567_A_C","MT_574_A_C","MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C") #These are the ones to exclude for the endometrial analysis

##CHOOSE DATASET AND START ANALYSIS
all_cohorts=c("NW","blood","HL","SO","PR","LM","lymph","KY")
all_mito_datasets<-lapply(all_cohorts,function(dataset) {
  if(dataset=="blood") {
    dataset_mito_data<-readRDS("../mito_data.Rds")
    CN_correlating_muts<-readRDS("../CN_correlation.RDS")
  } else {
    CN_correlating_muts=readRDS(paste0("CN_correlation_",dataset,".RDS"))
    dataset_mito_data=readRDS(paste0("mito_mutation_data_",dataset,".RDS"))
  }
  
  rho_cut_off=0
  
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


## DEFINE BESPOKE FUNCTIONS FOR THIS SCRIPT

#### Set up function to check the phasing between two mutations
check_mtDNA_mut_phasing=function(mut1,mut2,vaf_mat,tol1=0.03,tol2=0.03) {
  
  if(!mut1%in%rownames(vaf_mat)|!mut2%in%rownames(vaf_mat)) {stop(cat("Either/ both mutations not found in the vaf matrix. Please check mutation references."))}
  if(mut1==mut2) {stop(return(NA))}
    
  #Perform test 1: Do the vafs of the two mutations sum to >1 in any individual cell
  sum_of_vaf<-sum(colSums(vaf_mat[c(mut1,mut2),])>(1+tol1))
  test1_evidence=max(colSums(vaf_mat[c(mut1,mut2),])-1) #This is the maximum value by which the sum of VAFs exceeds 1
  same_haplo_res<-sum_of_vaf>0
  
  #Perform test 2: is there evidence of one mutation being 'sub-clonal' to the other i.e. one VAF is always at or equal to the other VAF
  mut1_higher=sum((vaf_mat[mut1,]-vaf_mat[mut2,])>tol2)
  mut2_higher=sum((vaf_mat[mut2,]-vaf_mat[mut1,])>tol2)
  test2_evidence<-min(max((vaf_mat[mut1,]-vaf_mat[mut2,])),max((vaf_mat[mut2,]-vaf_mat[mut1,]))) #This is the minimum value by which, for each mutation, the maximum by which the VAF of that mutation is exceeded by that of the other in a given cell.
  diff_haplo_res<-mut1_higher>0 & mut2_higher>0
  
  if(same_haplo_res&diff_haplo_res) {
    if(test1_evidence>test2_evidence) {
      diff_haplo_res<-F
    } else {
      same_haplo_res<-F
    }
  }
  
  res<-ifelse(same_haplo_res&!diff_haplo_res,"Same haplotype",ifelse(diff_haplo_res&!same_haplo_res,"Different haplotype",ifelse(same_haplo_res & diff_haplo_res,"Inconsistent","Inconclusive")))
  
  if(res=="Same haplotype") {
    mut1_mean_vaf<-mean(unlist(vaf_mat[mut1,]))
    mut2_mean_vaf<-mean(unlist(vaf_mat[mut2,]))
    
    res<-ifelse(mut1_mean_vaf>mut2_mean_vaf,"Ancestral","Descendant")
  }
  return(res)
}

get_descendant_mut_names=function(tree,node) {
  descendant_node_numbers<-phangorn::Descendants(tree,node=node,type = "all")
  descendant_tips<-descendant_node_numbers[descendant_node_numbers<=length(tree$tip.label)]
  descendant_nodes<-descendant_node_numbers[descendant_node_numbers>length(tree$tip.label)]
  descendant_muts<-c(tree$node.label[descendant_nodes-length(tree$tip.label)],tree$tip.label[descendant_tips])
  return(descendant_muts)
}

get_ancestral_mut_names=function(tree,node) {
  ancestral_node_numbers<-phangorn::Ancestors(tree,node=node,type = "all")
  ancestral_muts<-tree$node.label[ancestral_node_numbers-length(tree$tip.label)]
  return(ancestral_muts)
}

get_direct_descendants<-function(phasing_mat,mut) {
  descendant_muts<-colnames(phasing_mat)[phasing_mat[mut,]=="Ancestral"]
  descendant_muts<-descendant_muts[!is.na(descendant_muts)]
  sub_phasing_mat<-phasing_mat[descendant_muts,descendant_muts,drop=F]
  rownames(sub_phasing_mat)[rowSums(sub_phasing_mat=="Descendant",na.rm = T)==0]
}

get_all_descendants<-function(phasing_mat,mut) {
  descendant_muts<-colnames(phasing_mat)[phasing_mat[mut,]=="Ancestral"]
  return(descendant_muts[!is.na(descendant_muts)])
}

rename_node=function(tree,tip_set,new_name) {
  matches_sub_tree=sapply(unique(tree$edge[,1]),function(node) {
    descendant_nodes<-phangorn::Descendants(tree,node=node,type="children")
    if(length(descendant_nodes)!=length(tip_set)) {stop(return(F))}
    all(tree$tip.label[descendant_nodes] == tip_set)
  })
  mut_node=unique(tree$edge[,1])[matches_sub_tree]
  tree$node.label[mut_node-length(tree$tip.label)]<-new_name
  return(tree)
}

#Build a tree based on these phasing results
#Iterative function to build a genotype tree from the phasing results
build_mito_haplotype_tree=function(tree=NULL,muts,vaf_mat,phasing_mat) {
  
  if(is.null(tree)) {
    cat("Initiating tree",sep="\n")
    tree<-ape::stree(length(muts),type="star",tip.label=muts)
    tree$edge.length<-rep(1,length(muts))
    tree$node.label<-"ROOT"
  }
  
  for(mut in muts) {
    cat(mut,sep="\n")
    direct_descendant_muts<-get_direct_descendants(phasing_mat,mut)
    if(length(direct_descendant_muts)==0) {next}
    
    #Create subtree of the direct descendants of this mut
    if(length(direct_descendant_muts)==1) {
      #If there is only 1 descendant mutation need to create a dummy blank tip, otherwise it will get collapsed
      sub_tree<-stree(2,type="star",tip.label=c(direct_descendant_muts,""))
      sub_tree$edge.length<-c(1,0) #Add edge lengths
    } else {
      sub_tree<-stree(length(direct_descendant_muts),type="star",tip.label=direct_descendant_muts)
      sub_tree$edge.length<-rep(1,length(direct_descendant_muts)) #Add edge lengths
    }
    
    #Bind this to the old mutation
    new_tree<-bind.tree(x = tree,y=sub_tree,where=which(tree$tip.label==mut))
    
    #Label the newly created node with the mut_ref
    new_tree<-rename_node(new_tree,tip_set = sub_tree$tip.label,new_name = mut)
    
    #Remaining muts to include
    remaining_muts<-c(mut,get_all_descendants(phasing_mat = phasing_mat,mut = mut))
    sub_phasing_mat<-phasing_mat[remaining_muts,remaining_muts,drop=F]
    tree<-build_mito_haplotype_tree(tree=new_tree,muts=direct_descendant_muts,vaf_mat=vaf_mat,phasing_mat = sub_phasing_mat)
  }
  return(tree)
}


temp=lapply(all_cohorts,function(cohort) {
  cat(cohort,sep="\n")
  dataset_mito_data<-all_mito_datasets[[cohort]]
  dataset_samples<-names(dataset_mito_data)
  temp2=lapply(dataset_samples,function(exp_ID) {
    cat(exp_ID,sep="\n")
    sample_mito_data<-dataset_mito_data[[exp_ID]]
    vaf.filt<-sample_mito_data$matrices$vaf.filt
    
    #Get mut_refs of mutations that are >3% in at least one sample
    high_vaf_mut_refs<-rownames(vaf.filt)[rowSums(vaf.filt>0.03)>0]
    
    #Set file paths for saving
    phasing_file_dir="/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/phasing_files/"
    phasing_mat_path<-paste0(phasing_file_dir,"phasing_mat_",exp_ID,".Rds")
    tree_file_path<-paste0(phasing_file_dir,"tree_mtDNA_genotype_",exp_ID,".tree")
    tree_pdf_path<-paste0(phasing_file_dir,"genotype_tree_",exp_ID,".pdf")
    
    force_rerun=F
    if(file.exists(phasing_mat_path)&!force_rerun) {
      mat<-readRDS(phasing_mat_path)
    } else {
      mat=matrix(NA,nrow=length(high_vaf_mut_refs),ncol = length(high_vaf_mut_refs),dimnames = list(high_vaf_mut_refs,high_vaf_mut_refs))
      
      switcheroo=c("Ancestral","Descendant")
      names(switcheroo)<-rev(switcheroo)
      
      #Build the phasing matrix
      for(i in 1:nrow(mat)) {
        cat(i,sep="\t")
        for(j in 1:ncol(mat)) {
          if(i>=j) {
            next
          } else {
            mat[i,j]<-check_mtDNA_mut_phasing(mut1=rownames(mat)[i],mut2=rownames(mat)[j],vaf_mat=vaf.filt)
            mat[j,i]<-ifelse(mat[i,j]%in%c("Ancestral","Descendant"),switcheroo[mat[i,j]],mat[i,j])
          }
        }
      }
      
      #Save the phasing mat (this takes a long time to make)
      saveRDS(mat,file=phasing_mat_path)
    }
    
    #Start with all the mutations that aren't descendants of any other branches (i.e. that must therefore come from the 'root' genotype)
    muts_ancestral<-rownames(mat)[rowSums(mat=="Descendant",na.rm=T)==0]
    
    #Then build the tree from here & save
    tree<-build_mito_haplotype_tree(muts=muts_ancestral,vaf_mat=vaf.filt,phasing_mat=mat)
    write.tree(tree,file=tree_file_path)
    
    #plot the genotype phylogeny
    pdf(tree_pdf_path,height=15,width=4)
    plot(tree,show.node.label=T,cex=0.3)
    dev.off()
    
  })
})

# mat_col<-mat
# mat_col[mat=="Same haplotype"]<-1
# mat_col[mat=="Different haplotype"]<-0
# mat_col[mat=="Inconclusive"]<-0.5
# mat_col[mat=="Inconsistent"]<-NA
# for(i in 1:nrow(mat_col)) {mat_col[i,i] <-NA}
# 
# mat_col<-apply(mat_col, 2, as.numeric)
# dimnames(mat_col)<-dimnames(mat)
# 
# library(pheatmap)
# pdf("temp.pdf",width=20,height=20)
# pheatmap(mat_col,scale = "none",na_col = "lightgrey")
# dev.off()
# 



pdf("test.pdf")
plot(tree,show.node.label=T,cex=0.4)
dev.off()





# ###### PREVIOUS ATTEMPT - ABANDONED #################
# 
# for(i in 1:length(muts_ancestral_first)) {
#   mut<-muts_ancestral_first[i]
#   cat(mut,sep="\n")
#   if(is.null(tree)) {
#     cat("Initiating tree")
#     tree<-ape::stree(2,type="star",tip.label=c("Ancestral",mut))
#     tree$edge.length<-c(1,1)
#     tree$node.label<-"ROOT"
#   }
#   
#   #If the mutation is not yet in the tree, start by adding it to the root
#   if(!mut%in%c(tree$tip.label,tree$node.label)) {
#     cat("Adding mutation to root as starting point",sep="\n")
#     tree<-TreeTools::AddTip(tree=tree,where = 1+length(tree$tip.label),label = mut,edgeLength = 1)
#   }
#   
#   same_haplo_muts<-rownames(mat)[mat[mut,]%in%c("Ancestral","Descendant")]
#   same_haplo_muts<-same_haplo_muts[!is.na(same_haplo_muts)]
#   
#   if(length(same_haplo_muts)>0) {
#     test_mut_meanvaf<-mean(unlist(vaf_mat[mut,]))
#     for(j in 1:length(same_haplo_muts)) {
#       other_mut=same_haplo_muts[j]
#       cat(paste("Testing against",other_mut),sep="\n")
#       comparator_mut_meanvaf<-mean(unlist(vaf_mat[other_mut,]))
#       relationship=ifelse(test_mut_meanvaf>comparator_mut_meanvaf,"A","D")
#       
#       
#       #Get the current node/ tip number of the test mut
#       if(mut%in%tree$node.label) {
#         current_mut_node<-length(tree$tip.label)+which(tree$node.label==mut)
#       } else if(mut%in%tree$tip.label) {
#         current_mut_node<-which(tree$tip.label==mut)
#       }
#       
#       #If the 'other mut' is already in the tree, check the relationship to the test mut is correct
#       #If it isn't correct, need to rearrange the tree accordingly
#       
#       if(other_mut%in%c(tree$tip.label,tree$node.label)) {
#         descendant_muts<-get_descendant_mut_names(tree,node=current_mut_node)
#         ancestral_muts<-get_ancestral_mut_names(tree,node=current_mut_node)
#         
#         if(other_mut%in%descendant_muts) {
#           current_relationship<-"A"
#         } else if(other_mut%in%ancestral_muts) {
#           current_relationship<-"D"
#         } else {
#           current_relationship<-"DIFFERENT HAPLOTYPES"
#         }
#         
#         if(current_relationship==relationship) {
#           print("CORRECT RELATIONSHIP IN TREE ALREADY")
#           next
#         } else {
#           print("INCORRECT CURRENT RELATIONSHIP IN TREE - need to investigate")
#           next
#         }
#           
#       } else {
#         
#         
#         
#         #If the 'other mut' is not already in the tree, add it in the appropriate place
#         if(mut%in%tree$tip.label) {
#           if(relationship=="A") {
#             
#             #If new mut is a descendant of the test mut, bind it as a new tip
#             updated_tree<-TreeTools::AddTip(tree=tree,where = current_mut_node,label = other_mut,nodeLabel=mut,lengthBelow = 1,edgeLength = 1)
#             updated_tree$tip.label[updated_tree$tip.label==mut]<-""
#             #updated_tree$edge.length[which(updated_tree$tip.label=="")]<-0
#             updated_tree$node.label[which(updated_tree$node.label=="")]<-mut
#             
#             tree<-updated_tree
#           } else if(relationship=="D") {
#             #Bit more complicated if need to 'insert' the new mutant ancestral to the test mut
#             cat("This bit of code not yet sorted")
#           }
#         } else if (mut%in%tree$node.label) {
#           
#           if(relationship=="A") {
#             #If other mut is a descendant of the test mut, just bind it as a new tip
#             updated_tree<-phytools::bind.tip(tree,where=current_mut_node,tip.label = other_mut,edge.length = 1)
#             
#             tree<-updated_tree
#           } else if(relationship=="D") {
#             #Bit more complicated if need to 'insert' the new mutant ancestral to the test mut
#             cat("This bit of code not yet sorted")
#           }
#         }
#       }
#     }
#   } else {
#     #No other mutations on the same haplo type
#     next
#   }
# }
# 
