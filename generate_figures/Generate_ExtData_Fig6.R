#-----------------------------------------------------------------------------------#
# Generate_ExtData_Fig6.R
#
# Extended Data Fig. 6 | Heteroplasmic oocyte mutations.
#
#   a  Proportion of individuals inferred to carry at least one heteroplasmic
#      oocyte mutation above the heteroplasmy threshold on the x-axis. Because
#      sensitivity to low heteroplasmy is limited, values below 10% are lower
#      bounds.
#   b  Number of individuals against the number of heteroplasmic oocyte mutations
#      inferred at >1% heteroplasmy.
#   c  Distribution of heteroplasmic oocyte mutations across functional
#      categories - NOT built here, see the note at the foot of this script.
#
# Panel code follows mtDNA_mutations_comparator_tissues.Rmd, chunk
# heteroplasmic_oocyte_mutations.
#
# The principle: a mutation shared by clones whose common ancestor lies back in
# early development is more likely to have been present in the oocyte than to
# have arisen independently in each lineage.
#-----------------------------------------------------------------------------------#

cran_packages=c("dplyr","tidyr","ggplot2","ape","phytools")
for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

options(stringsAsFactors = F)

source(here::here("config.R")) #sets root_dir, genomeFile, project paths and my_theme
#plots_dir, rebuttal_figs_dir and my_theme all come from config.R
dir.create(paste0(plots_dir,"Extended_Data_Figure_06"),showWarnings=FALSE,recursive=TRUE)
ed6_dir=paste0(plots_dir,"Extended_Data_Figure_06/")

source(paste0(root_dir,"/data/mito_mutations_blood_functions.R")) #detect_het_oocyte_mutation, find_latest_acquisition_node, get_expanded_clade_nodes, getTips

#One independent sample is drawn at random per post-developmental clade, so fix
#the seed to make the inferred oocyte VAFs reproducible.
set.seed(42)

all_cohorts<-c("KY","HL","SO","PR","LM","NW","lymph","blood")

#Blacklisted recurrent artefacts
exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_456_C_T","MT_567_A_C","MT_574_A_C",
               "MT_8270_C_T","MT_16170_A_C","MT_16181_A_C","MT_16182_A_C",
               "MT_16183_A_C","MT_16189_T_C")

#The foetal and cord donors are analysed separately and excluded from this figure
foetal_ids<-c("8 pcw","18 pcw","8pcw","18pcw","CB001","CB002")

#-----------------------------------------------------------------------------------#
# Identify heteroplasmic oocyte mutations across the cohorts
#-----------------------------------------------------------------------------------#

all_het_oocyte_mut_df<-dplyr::bind_rows(lapply(all_cohorts,function(dataset) {
  f<-paste0(root_dir,"/data/nonblood/mito_mutation_data_",dataset,".RDS")
  if(!file.exists(f)) {cat("  no data for",dataset,"- skipped\n"); return(NULL)}
  dataset_mito_data<-readRDS(f)

  dplyr::bind_rows(Map(list=dataset_mito_data,Exp_ID=names(dataset_mito_data),function(list,Exp_ID){
    #Only assess donors with at least 8 colonies - fewer gives too little power to
    #tell a shared oocyte mutation from independent acquisition
    if(length(list$tree$tip.label)<8) return(NULL)

    #Thresholds differ for the foetal donors, which are analysed separately
    if(Exp_ID%in%foetal_ids) {
      threshold_vaf<-0.025; threshold_molecular_time<-5; node_height_for_independence<-20
    } else {
      threshold_vaf<-0.25; threshold_molecular_time<-30; node_height_for_independence<-100
    }

    het_oocyte_muts<-detect_het_oocyte_mutation(matrices=list$matrices,tree=list$tree,vaf_cutoff=threshold_vaf)
    het_oocyte_muts<-het_oocyte_muts[!het_oocyte_muts%in%exclude_muts & !grepl("DEL|INS",het_oocyte_muts)]
    if(length(het_oocyte_muts)==0) return(NULL)

    ho_vaf_mat<-list$matrices$vaf[het_oocyte_muts,list$tree$tip.label,drop=F]

    #How far back the most recent common ancestor of the positive samples sits.
    #An early (low molecular time) MRCA is what marks a mutation as oocyte-derived.
    MRCA_mol_times<-sapply(1:nrow(ho_vaf_mat),function(i) {
      pos_samples<-colnames(ho_vaf_mat)[ho_vaf_mat[i,]>threshold_vaf]
      MRCA_node<-find_latest_acquisition_node(tree=list$tree,pos_samples=pos_samples)
      phytools::nodeheight(tree=list$tree,node=MRCA_node)
    })

    keep<-MRCA_mol_times<threshold_molecular_time
    high_vaf_het_oocyte_muts<-het_oocyte_muts[keep]
    if(length(high_vaf_het_oocyte_muts)==0) return(NULL)
    ho_vaf_mat<-ho_vaf_mat[keep,,drop=FALSE]

    #Estimate the oocyte VAF from one sample per independent post-developmental
    #clade, so that a single expanded clone cannot dominate the average
    if(dataset=="lymph") {
      samples_to_include<-list$tree$tip.label
    } else {
      post_dev_nodes<-get_expanded_clade_nodes(tree=list$tree,height_cut_off=node_height_for_independence,
                                               min_samples=1,min_clonal_fraction=0)
      samples_to_include<-unlist(lapply(post_dev_nodes$nodes,function(node) sample(size=1,x=getTips(list$tree,node))))
    }
    samples_to_include<-intersect(samples_to_include,colnames(ho_vaf_mat))
    if(!length(samples_to_include)) return(NULL)

    ml_vaf<-sapply(high_vaf_het_oocyte_muts,function(mut_ref) {
      mean(unlist(ho_vaf_mat[mut_ref,samples_to_include,drop=T]))
    })

    data.frame(dataset=dataset,exp_ID=Exp_ID,
               mut_ref=high_vaf_het_oocyte_muts,
               n_pos=sapply(1:nrow(ho_vaf_mat),function(i) sum(ho_vaf_mat[i,]>threshold_vaf)),
               MRCA_molecular_time=MRCA_mol_times[keep],
               ml_vaf=ml_vaf)
  }))
}))

#Every donor that met the >=8 colony bar - the denominator for both panels
individuals_assessed<-dplyr::bind_rows(lapply(all_cohorts,function(dataset) {
  f<-paste0(root_dir,"/data/nonblood/mito_mutation_data_",dataset,".RDS")
  if(!file.exists(f)) return(NULL)
  d<-readRDS(f)
  dplyr::bind_rows(Map(list2=d,exp_ID=names(d),function(list2,exp_ID){
    if(length(list2$tree$tip.label)<8) return(NULL)
    data.frame(exp_ID=exp_ID,n_samp=length(list2$tree$tip.label))
  }))
}))

all_het_oocyte_mut_df<-all_het_oocyte_mut_df%>%filter(!exp_ID%in%foetal_ids)
individuals_assessed<-individuals_assessed%>%filter(!exp_ID%in%foetal_ids)

cat("Heteroplasmic oocyte mutations:",nrow(all_het_oocyte_mut_df),
    "across",length(unique(all_het_oocyte_mut_df$exp_ID)),"individuals\n")
cat("Individuals assessed:",nrow(individuals_assessed),"\n")

#-----------------------------------------------------------------------------------#
# Fig ED6a | Proportion of individuals with an oocyte mutation above a threshold
#-----------------------------------------------------------------------------------#

thresholds<-seq(0.005,0.9,0.005)
n_samples_with_ho_mut_over_threshold<-sapply(thresholds,function(threshold){
  all_het_oocyte_mut_df%>%filter(ml_vaf>threshold)%>%pull(exp_ID)%>%unique()%>%length()
})

prop_of_samples_with_ho_mut<-data.frame(threshold=thresholds,n=n_samples_with_ho_mut_over_threshold)%>%
  ggplot(aes(x=threshold,y=n/nrow(individuals_assessed)))+
  geom_line()+
  theme_classic()+
  scale_y_continuous(limits=c(0,0.35))+
  scale_x_log10(limits=c(0.005,1),breaks=c(0.005,0.01,0.02,0.05,0.1,0.2,0.5))+
  my_theme+
  labs(x="Heteroplasmy level",
       y="Proportion of individuals with\n at least one heteroplasmic oocyte\nmutation above threshold")

ggsave(filename=paste0(ed6_dir,"ExtDataFig6a.prop_of_samples_with_ho_mut.pdf"),prop_of_samples_with_ho_mut,width=2.5,height=2)
cat("ED6a: written\n")

#-----------------------------------------------------------------------------------#
# Fig ED6b | Number of oocyte mutations per individual at >1% heteroplasmy
#
# Individuals with none still count, hence the right_join onto the full assessed
# list followed by filling in zero.
#-----------------------------------------------------------------------------------#

n_of_het_oocyte_muts_plot<-all_het_oocyte_mut_df%>%
  tidyr::separate(mut_ref,into=c("chr","pos","ref","mut"),sep="_",remove=FALSE)%>%
  dplyr::filter(!grepl("\\*",mut) & ml_vaf>0.01)%>%
  group_by(exp_ID)%>%
  dplyr::summarise(nmut=n(),.groups="drop")%>%
  right_join(individuals_assessed,by="exp_ID")%>%
  tidyr::replace_na(list(nmut=0))%>%
  dplyr::group_by(nmut)%>%
  dplyr::summarise(n_samples=n(),.groups="drop")%>%
  ggplot(aes(x=nmut,y=n_samples))+
  geom_bar(stat="identity",col="black",linewidth=0.2,fill="lightblue")+
  theme_classic()+
  my_theme+
  labs(x="Number of heteroplasmic oocyte mutations\nwith VAF inferred >1%",y="Number of individuals")

ggsave(filename=paste0(ed6_dir,"ExtDataFig6b.n_of_het_oocyte_muts_plot.pdf"),n_of_het_oocyte_muts_plot,width=2,height=1.8)
cat("ED6b: written\n")

#-----------------------------------------------------------------------------------#
# PANEL c | not generated here
#
# The functional-category comparison (Mut_cat_comparison in the cross-tissue
# notebook) sets the oocyte mutations against ALL somatic mutations, so it needs
# complete.annotated.mutation.table - the full cross-tissue somatic mutation
# table, annotated with dndscv against the mtDNA reference (data/mtref.rda) and
# with mitovizR region coordinates. That table is the main product of the
# notebook rather than a cached object, so it is not lifted here.
#-----------------------------------------------------------------------------------#

cat("\nExtended Data Fig. 6 panels a and b written to",ed6_dir,"\n")
cat("Panel c not generated - see note in this script.\n")
