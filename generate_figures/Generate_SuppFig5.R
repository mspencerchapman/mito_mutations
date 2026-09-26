#-----------------------------------------------------------------------------------#
# Generate_SuppFig5.R
#
# Supplementary Fig. 5a, b, e | Heteroplasmic oocyte mutations (Supplementary Note 6).
# Panels c, d and f are not built here - see the notes further down.
#
#   a  Proportion of individuals inferred to carry at least one heteroplasmic
#      oocyte mutation above the heteroplasmy threshold on the x-axis. Because
#      sensitivity to low heteroplasmy is limited, values below 10% are lower
#      bounds.
#   b  Number of individuals against the number of heteroplasmic oocyte mutations
#      inferred at >1% heteroplasmy.
#   c  Distribution of heteroplasmic oocyte mutations across functional
#      categories - NOT built here, see the note at the foot of this script.
#   e  Phylogenies of four example individuals with a heatmap of their inferred
#      heteroplasmic oocyte mutations beneath.
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
dir.create(paste0(plots_dir,"Supplementary_Figures"),showWarnings=FALSE,recursive=TRUE)
supp_dir=paste0(plots_dir,"Supplementary_Figures/")

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
# Supp Fig 5a | Proportion of individuals with an oocyte mutation above a threshold
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

ggsave(filename=paste0(supp_dir,"SuppFig5a.prop_of_samples_with_ho_mut.pdf"),prop_of_samples_with_ho_mut,width=2.5,height=2)
cat("Supp Fig 5a: written\n")

#-----------------------------------------------------------------------------------#
# Supp Fig 5b | Number of oocyte mutations per individual at >1% heteroplasmy
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

ggsave(filename=paste0(supp_dir,"SuppFig5b.n_of_het_oocyte_muts_plot.pdf"),n_of_het_oocyte_muts_plot,width=2,height=1.8)
cat("Supp Fig 5b: written\n")

#-----------------------------------------------------------------------------------#
# Annotate the oocyte mutations | coding consequence and genome region
#
# Needed by panels c and d. dndscv supplies the protein consequence; mitovizR
# supplies the region boundaries (tRNA, rRNA, coding, D-loop).
#-----------------------------------------------------------------------------------#

suppressMessages({library(dndscv); library(mitovizR)})

all_mtDNA_genes<-c("MT-CYB","MT-ND5","MT-ND2","MT-ND4","MT-ND1","MT-CO3","MT-ATP6",
                   "MT-ND3","MT-ATP8","MT-ND4L","MT-CO2","MT-CO1","MT-ND6")
mtref_rda_path<-paste0(root_dir,"/data/mtref.rda")

oocyte_variant_data<-all_het_oocyte_mut_df%>%
  tidyr::separate(mut_ref,c("chr","pos","ref","mut"),"_")%>%
  mutate(pos=as.numeric(pos))%>%
  dplyr::select("sampleID"=exp_ID,chr,pos,ref,mut)%>%
  dplyr::filter(!duplicated(.))%>%
  arrange(sampleID,pos)

oocyte_dndsout<-suppressMessages(dndscv(oocyte_variant_data,gene_list=all_mtDNA_genes,
  refdb=mtref_rda_path,numcode=2,max_coding_muts_per_sample=Inf,max_muts_per_gene_per_sample=Inf))

#Region boundaries from mitovizR, used for both the functional categories in c
#and the circular layout in d
mito_domain_convert<-c("tRNA","rRNA","Coding","D-loop")
names(mito_domain_convert)<-c("trna","rrna","cds","reg")
mito_coords_reference_df<-mitovizR:::mito_df()%>%
  dplyr::select(type,ymin,ymax)%>%
  dplyr::mutate(Mutation_type=mito_domain_convert[type])
region_of<-function(pos) sapply(pos,function(x)
  mito_coords_reference_df$Mutation_type[mito_coords_reference_df$ymin<x & mito_coords_reference_df$ymax>=x])

all_het_oocyte_mut_df_annotated<-all_het_oocyte_mut_df%>%
  tidyr::separate(mut_ref,into=c("chr","pos","ref","mut"),sep="_")%>%
  mutate(pos=as.numeric(pos))%>%
  left_join(oocyte_dndsout$annotmuts,by=c("chr","pos","ref","mut"),relationship="many-to-many")%>%
  tidyr::replace_na(replace=list(impact="Non-coding"))
all_het_oocyte_mut_df_annotated$Mutation_type<-region_of(all_het_oocyte_mut_df_annotated$pos)

#-----------------------------------------------------------------------------------#
# PANEL c | Functional categories, oocyte mutations against all somatic mutations
#
# The somatic comparator is the complete cross-tissue annotated mutation table.
# Rebuilding it here would mean repeating the whole per-tissue dndscv chain, so
# Generate_ExtData_Fig5.R caches it and this script reads it back.
#-----------------------------------------------------------------------------------#

somatic_tbl_path<-paste0(root_dir,"/data/complete_annotated_mutation_table.Rds")
if(!file.exists(somatic_tbl_path)) {
  cat("Supp Fig 5c: no cached somatic mutation table - run Generate_ExtData_Fig5.R first - skipped\n")
} else {
  complete.annotated.mutation.table<-readRDS(somatic_tbl_path)
  complete.annotated.mutation.table$Mutation_type<-region_of(complete.annotated.mutation.table$pos)

  het_oocyte_mut_type_summary<-all_het_oocyte_mut_df_annotated%>%
    dplyr::filter(!grepl("\\*",mut))%>%
    dplyr::mutate(final_type=ifelse(is.na(impact)|impact=="Non-coding",Mutation_type,impact),
                  cat=ifelse(ml_vaf>=0.01,"Heteroplasmic\noocyte\n(VAF≥1%)","Heteroplasmic\noocyte\n(VAF<1%)"))%>%
    group_by(cat,final_type)%>%
    dplyr::summarise(n=n(),.groups="drop_last")%>%
    dplyr::mutate(prop=n/sum(n))%>%
    ungroup()

  Mut_type_levels<-c("D-loop","Synonymous","tRNA","rRNA","Missense","Stop_loss","Nonsense","Inconsistent\nannotation")
  somatic_summary<-complete.annotated.mutation.table%>%
    dplyr::mutate(final_type=ifelse(impact=="Non-Coding" & Mutation_type=="Coding","Inconsistent\nannotation",
                             ifelse(is.na(impact)|impact=="Non-Coding",Mutation_type,impact)))%>%
    group_by(final_type)%>%
    dplyr::summarise(n=n(),.groups="drop")%>%
    dplyr::mutate(prop=n/sum(n),cat="Somatic")

  Mut_cat_comparison<-bind_rows(het_oocyte_mut_type_summary,somatic_summary)%>%
    mutate(final_type=factor(final_type,levels=Mut_type_levels))%>%
    ggplot(aes(x=cat,y=prop,fill=forcats::fct_rev(final_type)))+
    geom_bar(stat="identity",position="fill",col="black",linewidth=0.2)+
    theme_classic()+
    scale_y_continuous(breaks=seq(0,1,0.1))+
    scale_fill_brewer(palette="Set2",direction=-1)+
    my_theme+
    labs(fill="Mutation\ncategory",x="",y="Proportion")
  ggsave(filename=paste0(supp_dir,"SuppFig5c.mut_category_comparison.pdf"),Mut_cat_comparison,width=3.6,height=2.5)

  cat("Supp Fig 5c: written -",
      paste(sprintf("%s n=%d",
        c("oocyte <1%","oocyte >=1%","somatic"),
        c(sum(het_oocyte_mut_type_summary$n[grepl("<1",het_oocyte_mut_type_summary$cat)]),
          sum(het_oocyte_mut_type_summary$n[grepl("≥1",het_oocyte_mut_type_summary$cat)]),
          sum(somatic_summary$n))),collapse=", "),"\n")
}

#-----------------------------------------------------------------------------------#
# PANEL d | Positions of the oocyte mutations around the mitochondrial genome
#-----------------------------------------------------------------------------------#

mtDNA_mut_pos<-mitovizR::plot_df(all_het_oocyte_mut_df_annotated%>%
    dplyr::filter(!grepl("\\*",mut))%>%
    dplyr::mutate(SAMPLE="all")%>%
    dplyr::select(SAMPLE,pos,ref,mut,"HF"=ml_vaf),
  pos_col="pos",ref_col="ref",alt_col="mut")
ggsave(filename=paste0(supp_dir,"SuppFig5d.mtDNA_mut_positions.pdf"),mtDNA_mut_pos,width=5,height=5)
cat("Supp Fig 5d: written\n")

#-----------------------------------------------------------------------------------#
# PANEL e | Phylogenies with a heatmap of the heteroplasmic oocyte mutations
#
# Four example individuals. The cross-tissue notebook writes an equivalent plot
# for every donor as <exp_ID>_shared_muts.pdf, but showing ALL shared mutations;
# the published panel shows only the inferred heteroplasmic oocyte mutations, so
# it was previously subset by hand. Restricting the heatmap here makes the panel
# reproducible from the repository.
#-----------------------------------------------------------------------------------#

panel_e_donors<-data.frame(
  exp_ID =c("PD44890","PD44887","PD41857","PD5182"),
  dataset=c("PR","PR","LM","NW"),
  tissue =c("MUTYH-mutant colorectal crypts","MUTYH-mutant colorectal crypts",
            "endometrial glands","MPN blood colonies"))

#Matches the shared-mutation heatmaps in the cross-tissue notebook
lowest_VAF_to_show<-0.03
col_scheme<-c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9])(100))
names(col_scheme)<-seq(0,1,0.01)

for(i in 1:nrow(panel_e_donors)) {
  this_id<-panel_e_donors$exp_ID[i]
  list<-readRDS(paste0(root_dir,"/data/nonblood/mito_mutation_data_",
                       panel_e_donors$dataset[i],".RDS"))[[this_id]]

  #The oocyte mutations inferred for this donor, highest heteroplasmy first
  plot_muts<-all_het_oocyte_mut_df%>%filter(exp_ID==this_id)%>%
    arrange(desc(ml_vaf))%>%pull(mut_ref)
  plot_muts<-plot_muts[plot_muts%in%rownames(list$matrices$vaf)]
  if(!length(plot_muts)) {cat("Supp Fig 5e:",this_id,"- no oocyte mutations - skipped\n"); next}

  #Drop samples absent from the shearwater table (excluded for contamination)
  tree<-ape::keep.tip(list$tree,
          list$tree$tip.label[list$tree$tip.label%in%list$sample_shearwater_calls$sampleID])
  tree$coords<-NULL
  vaf.mtx<-(list$matrices$vaf*list$matrices$SW)[plot_muts,tree$tip.label,drop=FALSE]

  hm<-matrix(0,nrow=length(plot_muts),ncol=length(tree$tip.label),
             dimnames=list(plot_muts,tree$tip.label))
  for(j in seq_along(plot_muts)) {
    v<-vaf.mtx[plot_muts[j],]
    v[v<lowest_VAF_to_show]<-0
    hm[j,]<-col_scheme[as.character(round(v,digits=2))]
  }

  #Plain pdf(), not save_pdf(): the latter is gated by the notebooks' save_plots
  #flag and rescales the canvas, which shrinks base-R text. Figure scripts write
  #unconditionally at the stated size.
  grDevices::pdf(file=paste0(supp_dir,"SuppFig5e.",this_id,"_oocyte_muts.pdf"),width=7,height=2.6)
  tree<-plot_tree(tree=tree,cex.label=0,plot_axis=TRUE,vspace.reserve=1.1,
                  title=paste0(this_id," (",panel_e_donors$tissue[i],", n = ",
                               length(tree$tip.label)," samples)"))
  add_mito_mut_heatmap(tree=tree,heatmap=hm,border="gray",
                       heatmap_bar_height=0.1,cex.label=0.25)
  dev.off()
  cat("Supp Fig 5e:",this_id,"-",length(plot_muts),"oocyte mutation(s) written\n")
}

#-----------------------------------------------------------------------------------#
# PANEL f | mtDNA-defined clone inference in two of the panel e individuals
#
# Clone assignments were produced by the Seurat SNN clustering described in
# Supplementary Note 13 and are cached in data/mito_mut_clones/nonblood/; this
# section draws them onto the phylogenies.
#-----------------------------------------------------------------------------------#

clone_dir<-paste0(root_dir,"/data/mito_mut_clones/nonblood/")
panel_f_donors<-panel_e_donors%>%filter(exp_ID%in%c("PD5182","PD41857"))

cluster_cols<-c("lightgray","#1f77b4","#d62728","#2ca02c","#ff7f0e","#9467bd","#8c564b",
                "#e377c2","#7f7f7f","#bcbd22","#17becf","#ad494a","#e7ba52","#8ca252",
                "#756bb1","#636363","#aec7e8",RColorBrewer::brewer.pal(12,"Paired"))

for(i in 1:nrow(panel_f_donors)) {
  this_id<-panel_f_donors$exp_ID[i]
  f<-paste0(clone_dir,this_id,"_mtdna_clone_assignment.txt")
  if(!file.exists(f)) {cat("Supp Fig 5f:",this_id,"- no clone assignment - skipped\n"); next}

  list<-readRDS(paste0(root_dir,"/data/nonblood/mito_mutation_data_",
                       panel_f_donors$dataset[i],".RDS"))[[this_id]]
  exp_clones<-read.delim(f)
  n_clones<-length(unique(exp_clones$cluster_id))
  exp_cluster_cols<-cluster_cols[1:n_clones]
  names(exp_cluster_cols)<-unique(exp_clones$cluster_id)

  #The largest clone is drawn light grey, so swap it with whichever id is 0
  biggest<-names(table(exp_clones$cluster_id))[which.max(table(exp_clones$cluster_id))]
  was_0<-which(exp_clones$cluster_id==0)
  is_biggest<-which(as.character(exp_clones$cluster_id)==biggest)
  exp_clones$cluster_id[is_biggest]<-0
  exp_clones$cluster_id[was_0]<-as.integer(biggest)

  grDevices::pdf(file=paste0(supp_dir,"SuppFig5f.Inferred_clones_",this_id,".pdf"),width=7,height=3)
  par(mfrow=c(1,1))
  tree<-plot_tree(tree=list$tree,cex.label=0)
  clone_hm<-matrix(NA,nrow=1,ncol=length(tree$tip.label),
                   dimnames=list("Clones",tree$tip.label))
  for(k in 1:nrow(exp_clones)) {
    if(exp_clones$sample_id[k]%in%tree$tip.label)
      clone_hm[1,exp_clones$sample_id[k]]<-exp_cluster_cols[as.character(exp_clones$cluster_id[k])]
  }
  add_heatmap(tree=tree,heatmap=clone_hm,border="gray",cex.label=1)
  legend("topleft",inset=c(.01,.01),title="Clone no.",names(exp_cluster_cols),
         fill=exp_cluster_cols,horiz=FALSE,cex=0.7,ncol=1)
  dev.off()
  cat("Supp Fig 5f:",this_id,"-",n_clones,"clones written\n")
}

cat("\nSupplementary Fig. 5 panels a-f written to",supp_dir,"\n")
