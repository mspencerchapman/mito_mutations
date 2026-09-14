#-----------------------------------------------------------------------------------#
# Generate_ExtData_Fig3.R
#
# Extended Data Fig. 3 | Mitochondrial mutation signature decomposition by VAF.
#
#   a  The 192-profile mutational signatures extracted by a hierarchical
#      Dirichlet process. Mutations are aggregated by individual and VAF bin.
#      Positive bars are mutations on the heavy strand, negative bars the light
#      strand, coloured by substitution type. Signature N1 is the genuine one.
#   b  Absolute signature contributions to the mutation set in each VAF bin in
#      each individual, showing large artefactual contributions below ~0.5% VAF
#      and almost exclusively N1 above ~1%.
#   c  Mutational profile of the calls in the comparator (cross-tissue) datasets,
#      divided by strand.
#   d  As c, but restricted to lower VAF mutations.
#
# Panels a and b follow full_analysis_scripts/Mutational_signature_extraction_by_VAF.R.
# Panels c and d follow mtDNA_mutations_comparator_tissues.Rmd, chunks
# trinuc_profiles1 and trinuc_profiles2.
#
# The hdp extraction itself is a farm job (see
# full_analysis_scripts/Mutational_signature_extraction_bash/). This script does
# NOT re-run it: it reads the extracted components and exposures that the
# extraction wrote to data/mutational_signatures/.
#-----------------------------------------------------------------------------------#

cran_packages=c("dplyr","tidyr","tibble","ggplot2","stringr","RColorBrewer","ggh4x")
#trinucleotide_plot() additionally needs BSgenome/MutationalPatterns, loaded by the function itself
for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

options(stringsAsFactors = F)

source(here::here("config.R")) #sets root_dir, genomeFile, project paths and my_theme
plots_dir=paste0(root_dir,"/plots/")
dir.create(paste0(plots_dir,"Extended_Data_Figure_03"),showWarnings=FALSE,recursive=TRUE)
ed3_dir=paste0(plots_dir,"Extended_Data_Figure_03/")

source(paste0(root_dir,"/data/mito_mutations_blood_functions.R")) #supplies mitochondrial_extracted_signature_plot

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))+
  theme(legend.key.height=unit(3,"mm"),legend.title = element_text(size=8))

#-----------------------------------------------------------------------------------#
# Extracted signatures and their exposures
#
# components: one row per extracted signature, one column per 192-profile channel
#             (trinucleotide context x substitution, separately for the heavy (_H)
#             and light (_L) replication strand)
# exposures:  one column per <individual>_<VAF bin>, giving each signature's
#             proportional contribution to that mutation set
#-----------------------------------------------------------------------------------#

sig_dir<-paste0(root_dir,"/data/mutational_signatures/")
exposures<-read.csv(paste0(sig_dir,"exposures.csv"))
components<-read.csv(paste0(sig_dir,"components.csv"))
mutation_profiles_mat<-as.matrix(read.table(paste0(sig_dir,"trinuc_mut_mat.txt")))

#The csv files carry no row names and read.csv mangles the column names
#("C>A,A-A_H" becomes "C.A.A.A_H", and sample names starting with a digit gain an
#"X"). Rebuild them exactly as Mutational_signature_extraction_by_VAF.R does, so
#that the signature labels read N0..N6 and the mutation contexts read "C>A,A-A".
sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
colnames(components)<-c(paste(full_vec,"H",sep="_"),paste(full_vec,"L",sep="_"))
rownames(components)<-paste0("N",0:(nrow(components)-1))
rownames(exposures)<-paste0("N",0:(nrow(exposures)-1))
colnames(exposures)<-gsub("^X","",colnames(exposures))
colnames(exposures)<-gsub("\\.pcw","pcw",colnames(exposures))

#Signature colours, matched to panel b: that panel draws Signature with
#scale_fill_brewer("Paired") over this level order, so reproducing the order here
#keeps the two panels consistent.
sig_levels<-c(rownames(exposures)[nrow(exposures):3],"N0","N1")
sig_cols<-RColorBrewer::brewer.pal(max(3,length(sig_levels)),"Paired")[seq_along(sig_levels)]
names(sig_cols)<-sig_levels

cat("Signatures extracted:",nrow(components),"| mutation sets:",ncol(exposures),"\n")

#-----------------------------------------------------------------------------------#
# Fig ED3a | The extracted 192-profile signatures
#
# Split across two panels purely for legibility. N1 (the genuine signature) is
# also drawn on its own.
#-----------------------------------------------------------------------------------#

#' Add signature-coloured y-axis strips to a signature plot
#'
#' The facet strips on the right name the signature; colouring them with the same
#' scheme panel b uses makes the two panels readable together. ggh4x is used
#' because ggplot2 alone cannot style individual strips.
colour_signature_strips<-function(p,sigs) {
  cols<-unname(sig_cols[sigs]); cols[is.na(cols)]<-"grey85"
  p+ggh4x::facet_grid2(Signature~substitution,scales="free_y",
                       strip=ggh4x::strip_themed(
                         background_y=ggh4x::elem_list_rect(fill=cols),
                         text_y=ggh4x::elem_list_text(colour="black")))
}

#Panel height scales with the number of signatures so the profiles are not
#squashed vertically
sig_panel_height<-function(n) 1.25*n+0.9

sigs_1<-rownames(components)[1:min(4,nrow(components))]
pdf(paste0(ed3_dir,"ExtDataFig3a.extracted_signatures_1to4.pdf"),width=7,height=sig_panel_height(length(sigs_1)))
print(colour_signature_strips(mitochondrial_extracted_signature_plot(components[sigs_1,]),sigs_1))
dev.off()

if(nrow(components)>4) {
  sigs_2<-rownames(components)[5:nrow(components)]
  pdf(paste0(ed3_dir,"ExtDataFig3a.extracted_signatures_5to7.pdf"),width=7,height=sig_panel_height(length(sigs_2)))
  print(colour_signature_strips(mitochondrial_extracted_signature_plot(components[sigs_2,]),sigs_2))
  dev.off()
}

#N1 is the genuine signature
pdf(paste0(ed3_dir,"ExtDataFig3a.genuine_signature_N1.pdf"),width=7,height=sig_panel_height(1))
print(colour_signature_strips(mitochondrial_extracted_signature_plot(components["N1",,drop=FALSE]),"N1"))
dev.off()
cat("ED3a: written\n")

#-----------------------------------------------------------------------------------#
# Fig ED3b | Signature contributions by VAF bin and individual
#
# Exposures are proportions, so they are multiplied by the number of mutations in
# each set to give absolute counts.
#-----------------------------------------------------------------------------------#

VAF_groups=unique(sapply(stringr::str_split(colnames(exposures),pattern="_"),
                         function(vec) paste(vec[2:3],collapse="_")))
new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%",
                 "3.1-6.2%","6.2-12.5%","12.5-25%","25-50%","50-100%")[seq_along(VAF_groups)]
names(new_VAF_groups)<-VAF_groups

mut_numbers=data.frame(SampleID=rownames(mutation_profiles_mat),nmuts=rowSums(mutation_profiles_mat))

#Bins with too few mutations were not fitted by the hdp run. These are the high
#VAF bins, which are essentially pure N1, so they are filled in as such rather
#than left as gaps in the figure.
excluded_cats<-rownames(mutation_profiles_mat)[!rownames(mutation_profiles_mat)%in%colnames(exposures)]
all_N1=c(0,1,rep(0,nrow(exposures)-2))
excluded_mat<-matrix(rep(all_N1,times=length(excluded_cats)),nrow=length(all_N1),
                     dimnames=list(rownames(exposures),excluded_cats))

donor_levels<-c("8pcw","18pcw","CB001","CB002","KX001","KX002","SX001","AX001",
                "KX007","KX008","KX003","KX004")

mito_sigs_plot<-t(cbind(exposures,excluded_mat))%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="SampleID")%>%
  left_join(mut_numbers,by="SampleID")%>%
  mutate(exp_ID=stringr::str_split(SampleID,pattern="_",simplify=T)[,1])%>%
  mutate(exp_ID=factor(exp_ID,levels=donor_levels))%>%
  mutate(VAF_range=sapply(stringr::str_split(SampleID,pattern="_"),function(vec) paste(vec[2:3],collapse="_")))%>%
  mutate(VAF_range=factor(new_VAF_groups[VAF_range],levels=new_VAF_groups))%>%
  dplyr::select(-SampleID)%>%
  tidyr::gather(-VAF_range,-exp_ID,-nmuts,key="Signature",value="Exposure")%>%
  mutate(abs_muts=Exposure*nmuts)%>%
  mutate(Signature=factor(Signature,levels=sig_levels))%>%
  filter(!is.na(exp_ID),!is.na(VAF_range))%>%
  ggplot(aes(x=VAF_range,y=abs_muts,fill=Signature))+
  geom_bar(stat="identity",position="stack",col="black",linewidth=0.25)+
  theme_bw()+
  scale_fill_brewer(palette="Paired")+
  my_theme+
  theme(axis.text.x=element_text(size=5,angle=90),
        strip.text.x=element_text(size=7,margin=unit(c(1,0,1,0),"mm")))+
  facet_wrap(~exp_ID,ncol=6)+
  labs(x="VAF range",y="Number of mutations assigned")

ggsave(filename=paste0(ed3_dir,"ExtDataFig3b.signature_contributions_by_VAF.pdf"),mito_sigs_plot,width=8,height=4)
cat("ED3b: written\n")

#-----------------------------------------------------------------------------------#
# Fig ED3c and ED3d | Mutational profiles in the comparator datasets
#
# Observed/expected 192-profiles across the cross-tissue cohorts: c uses every
# retained call, d only the lower-VAF ones, where artefactual signatures dominate.
#
# The notebook chunk filters vaf<0.05, but its output file is named "under0.1"
# and the published caption says "<10%" - and 0.1 is what the published panel
# actually shows, so 0.1 is used here.
#-----------------------------------------------------------------------------------#

#trinucleotide_plot() reads the trinucleotide context for each call straight from
#the reference FASTA, which it picks up from this variable rather than taking it
#as an argument. Set it to your own copy before running.
#genomeFile is set in config.R
genomeFile=path.expand(genomeFile)

#It also expects the mtDNA trinucleotide frequency reference under this name,
#used to convert observed counts into observed/expected ratios
mtdna_trinuc_freq<-readRDS(paste0(root_dir,"/data/mtDNA_trinuc_freqs_coding_dloop_heavy_light.Rds"))

low_vaf_cutoff<-0.1

all_cohorts<-c("KY","HL","SO","PR","LM","NW","lymph","blood")
exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_456_C_T","MT_567_A_C","MT_574_A_C",
               "MT_8270_C_T","MT_16170_A_C","MT_16181_A_C","MT_16182_A_C",
               "MT_16183_A_C","MT_16189_T_C")
mutCN_cutoff=25
vaf_cut_off<-0.03

mito_cn=read.csv(paste0(root_dir,"/data/whole_genome_coverage_pileup_and_bedtools_annotated.csv"),header=T)

#Tidy table of retained calls across the comparator cohorts. Heteroplasmic oocyte
#mutations and copy-number-correlating artefacts are removed, as in the notebook.
all_df_tidy<-dplyr::bind_rows(lapply(all_cohorts,function(dataset) {
  f<-paste0(root_dir,"/data/nonblood/mito_mutation_data_",dataset,".RDS")
  if(!file.exists(f)) {cat("  no data for",dataset,"- skipped\n"); return(NULL)}
  dataset_mito_data<-readRDS(f)
  cn_muts_ds<-dataset_mito_data$CN_correlating_muts

  df<-dplyr::bind_rows(Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID) {
    if(is.null(list$matrices)) return(NULL)
    CN_removal<-list$matrices$CN_correlating_mut_removal_mat
    if(is.null(CN_removal)) return(NULL)
    implied_tidy<-list$matrices$implied_mutCN%>%as.data.frame()%>%
      tibble::rownames_to_column(var="mut_ref")%>%
      tidyr::gather(key="Sample",value="implied_mut_CN",-mut_ref)
    (list$matrices$vaf*list$matrices$SW*CN_removal)%>%
      as.data.frame()%>%
      tibble::rownames_to_column(var="mut_ref")%>%
      tidyr::gather(key="Sample",value="vaf",-mut_ref)%>%
      left_join(implied_tidy,by=c("Sample","mut_ref"))%>%
      dplyr::filter(!grepl("DEL|INS",mut_ref) & !mut_ref%in%exclude_muts &
                      !mut_ref%in%list$het_oocyte_muts)%>%
      dplyr::filter(vaf>=vaf_cut_off)%>%
      mutate(exp_ID=exp_ID)
  }))
  if(!nrow(df)) return(NULL)
  df%>%
    left_join(mito_cn%>%dplyr::select(Sample,bedtools_mtDNA_genomes)%>%
                dplyr::distinct(Sample,.keep_all=TRUE),by="Sample")%>% #one row per sample
    mutate(implied_mut_CN=vaf*bedtools_mtDNA_genomes)%>%
    dplyr::filter(!(mut_ref%in%cn_muts_ds & implied_mut_CN<=mutCN_cutoff))%>%
    mutate(dataset=dataset)
}))
cat("ED3c/d: ",nrow(all_df_tidy)," calls across ",
    length(unique(all_df_tidy$dataset))," comparator cohorts\n",sep="")

prep_for_profile<-function(df) {
  df%>%
    filter(Sample!="global" & !mut_ref%in%exclude_muts)%>%
    tidyr::separate(mut_ref,into=c("chr","pos","ref","mut"))%>%
    dplyr::rename("donor"=exp_ID)%>%
    mutate(pos=as.numeric(pos))
}


if(nrow(all_df_tidy) && file.exists(genomeFile)) {
  trinucleotide_plot(mutations=prep_for_profile(all_df_tidy),
                     analysis_region="all_mtDNA",analysis_type="obs_exp",
                     file_name=paste0(ed3_dir,"ExtDataFig3c.mutational_profile_comparator_all.pdf"))
  cat("ED3c: written\n")

  trinucleotide_plot(mutations=prep_for_profile(all_df_tidy%>%filter(vaf<low_vaf_cutoff)),
                     analysis_region="all_mtDNA",analysis_type="obs_exp",
                     file_name=paste0(ed3_dir,"ExtDataFig3d.mutational_profile_comparator_low_vaf.pdf"))
  cat("ED3d: written\n")
} else {cat("ED3c/d: no comparator data - skipped\n")}

cat("\nExtended Data Fig. 3 panels written to",ed3_dir,"\n")
