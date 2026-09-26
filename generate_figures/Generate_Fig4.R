#-----------------------------------------------------------------------------------#
# --------Load packages (and install if they are not installed yet)-------------------
#-----------------------------------------------------------------------------------#
cran_packages=c("devtools","ggridges","ape","stringr","dplyr","tidyr","ggplot2","gridExtra","phylosignal")
bioconductor_packages=c("MutationalPatterns")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if (!require("BiocManager", quietly = T, warn.conflicts = F))
  install.packages("BiocManager")
for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    BiocManager::install(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

if(!require("dndscv", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("im3sanger/dndscv")
  library("dndscv",character.only=T,quietly = T, warn.conflicts = F)
}

#-----------------------------------------------------------------------------------#
# ----------------------------------Set paths and import files------------------------
#-----------------------------------------------------------------------------------#

options(stringsAsFactors = F)

#Set these file paths before running the script
#genomeFile is set in config.R
source(here::here("config.R")) #sets root_dir, genomeFile, project paths and my_theme
source(paste0(root_dir,"/data/mito_mutations_blood_functions.R"))

#Set the key file paths using the root dir
tree_file_paths = list.files(paste0(root_dir,"/data/tree_files"),pattern=".tree",full.names = T)
ref_file=paste0(root_dir,"/data/Samples_metadata_ref.csv")
#plots_dir, rebuttal_figs_dir and my_theme all come from config.R


#Create the figure output directories if they do not already exist
for(d in c("Extended_Data_Figure_06")) dir.create(paste0(plots_dir,d),showWarnings=FALSE,recursive=TRUE)
for(d in c("Figure_04")) dir.create(paste0(plots_dir,d),showWarnings=FALSE,recursive=TRUE)

#Read in the mitochondrial copy number data
mito_cn=read.csv(paste0(root_dir,"/data/whole_genome_coverage_pileup_and_bedtools_annotated.csv"),header=T)

#Combine the sample level metadata info for the adult and foetal blood samples
sample_level_metadata_EM<-read.csv(paste0(root_dir,"/data/EM_sample_level_metadata.csv"))%>%dplyr::rename("exp_ID"="donor_id","Sample"="PDID","Cell_type"="cell_type","coverage"="mean_depth")
sample_level_metadata_foetal<-read.csv(paste0(root_dir,"/data/foetal_sample_level_metadata.csv"))%>%mutate(Sample=paste(Donor_ID,"hum",sep="_"))%>%dplyr::select(-Donor_ID,-Percentage)
sample_level_metadata<-dplyr::bind_rows(sample_level_metadata_EM,sample_level_metadata_foetal)

phenotype_data_EM<-read.csv(paste0(root_dir,"/data/Summary_pheno_pdid.csv"),stringsAsFactors = F)%>%
  dplyr::rename("exp_ID"=Donor_ID,"Sample"=PDID)

#Import individual level metadata for the adult and foetal blood samples
ref_df=read.csv(ref_file)%>%filter(Dataset!="Lymphocyte")
Individual_cols=RColorBrewer::brewer.pal(12,"Paired")
names(Individual_cols)<-ref_df$Sample[order(ref_df$Age)]

mito_cn<-mito_cn%>%
  left_join(sample_level_metadata,by="Sample",relationship="many-to-many")%>%
  left_join(ref_df%>%dplyr::rename("exp_ID"=Sample))
mito_cn$exp_ID<-factor(mito_cn$exp_ID,levels=ref_df$Sample[order(ref_df$Age)]) #Make the exp_ID a factor, with levels increasing by individual age
mito_cn<-left_join(mito_cn,phenotype_data_EM,relationship="many-to-many")

#Now import the mitochondrial mutation data
mito_data_file=paste0(root_dir,"/data/mito_data.Rds")
mito_data<-readRDS(mito_data_file)
CN_correlating_muts<-readRDS(paste0(root_dir,"/data/CN_correlation.RDS"))

#Define the 'black listed' mutation set - those with recurent artefacts despite the Shearwater filtering
exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_567_A_C","MT_574_A_C","MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C")

#-----------------------------------------------------------------------------------#
## Additional mutational signature info for this analysis----
#-----------------------------------------------------------------------------------#

exposures<-read.csv(paste0(root_dir,"/data/mutational_signatures/exposures.csv"))
components<-read.csv(paste0(root_dir,"/data/mutational_signatures/components.csv"))
mutation_profiles_mat=read.table(paste0(root_dir,"/data/mutational_signatures/trinuc_mut_mat.txt"))
sig_ref<-readRDS(file=paste0(root_dir,"/data/mutational_signatures/sig_ref_file.Rds"))

#Clean up the row and column names
sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G"); ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
colnames(components)<-c(paste(full_vec,"H",sep="_"),paste(full_vec,"L",sep="_"))
rownames(components)<-paste0("N",0:(nrow(components)-1))
rownames(exposures)<-paste0("N",0:(nrow(exposures)-1))
colnames(exposures)<-gsub("^X","",colnames(exposures))
colnames(exposures)<-gsub("\\.pcw"," pcw",colnames(exposures))

VAF_groups=unique(sapply(stringr::str_split(colnames(exposures),pattern="_"),function(vec) paste(vec[2:3],collapse="_")))
new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")
names(new_VAF_groups)<-VAF_groups

mut_numbers=data.frame(SampleID=rownames(mutation_profiles_mat),nmuts=rowSums(mutation_profiles_mat))

#Fill in the rare bins with too low mutation numbers - assume all are N1, as these are high VAF mutations
# (in samples with adequate mutations for extraction, high VAF mutations are all N1)
excluded_cats<-rownames(mutation_profiles_mat)[!rownames(mutation_profiles_mat)%in%colnames(exposures)]
all_N1=c(0,1,rep(0,nrow(exposures)-2))
excluded_mat<-matrix(rep(all_N1,times=length(excluded_cats)),nrow=length(all_N1),dimnames = list(rownames(exposures),excluded_cats))

real_muts_dist_df<-t(cbind(exposures,excluded_mat))%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="SampleID")%>%
  left_join(mut_numbers)%>%
  mutate(exp_ID=stringr::str_split(SampleID,pattern="_",simplify=T)[,1])%>%
  mutate(VAF_range=sapply(stringr::str_split(SampleID,pattern="_"),function(vec) paste(vec[2:3],collapse="_")))%>%
  mutate(VAF_range=factor(new_VAF_groups[VAF_range],levels=new_VAF_groups))%>%
  dplyr::select(-SampleID)%>%
  gather(-VAF_range,-exp_ID,-nmuts,key="Signature",value="Exposure")%>%
  mutate(abs_muts=Exposure*nmuts)%>%
  mutate(Signature=factor(Signature,levels=c(rownames(exposures)[nrow(exposures):3],"N0","N1")))%>%
  dplyr::filter(Signature=="N1")%>%
  mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))

real_mut_distribution_plot<-real_muts_dist_df%>%
  ggplot(aes(x=VAF_range,y=abs_muts))+
  geom_bar(stat="identity",fill="#FDBF6F",col="black",linewidth=0.25)+
  geom_line(aes(group=exp_ID))+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(size=5,angle=90),strip.text.x = element_text(size=6))+
  #facet_grid(rows=vars(factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)])),scales="free")+
  facet_wrap(~factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)]),scales="fixed",dir="v",ncol=6)+
  labs(x="VAF range",y="Number of N1 mutations detected")

#-----------------------------------------------------------------------------------#
### Generate FIG. 4A ---------
#-----------------------------------------------------------------------------------#

#Normalize mutation numbers by the number of samples included
nsamp_df<-Map(list=mito_data,Exp_ID=names(mito_data),function(list,Exp_ID) {
  data.frame(exp_ID=Exp_ID,n_samp=length(colnames(list$matrices$SW)))
})%>%dplyr::bind_rows()%>%
  mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))

real_mut_distribution_plot_normalized<-left_join(real_muts_dist_df,nsamp_df)%>%
  mutate(muts_per_samp=abs_muts/n_samp)%>%
  ggplot(aes(x=VAF_range,y=muts_per_samp))+
  geom_bar(stat="identity",fill="#FDBF6F",col="black",linewidth=0.25)+
  geom_line(aes(group=exp_ID))+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(size=5,angle=90),strip.text.x = element_text(size=6))+
  #facet_grid(rows=vars(factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)])),scales="free")+
  facet_wrap(~factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)]),scales="fixed",dir="v",ncol=6)+
  labs(x="VAF range",y="Number of N1 mutations detected per sample")
ggsave(filename=paste0(plots_dir,"Figure_04/Fig4a.real_mut_distribution_plot.pdf"),real_mut_distribution_plot_normalized,width=7,height=2.5)

#-----------------------------------------------------------------------------------#
##----------------------VISUALIZE WRIGHT-FISHER MODELS OF DRIFT----------------------
#-----------------------------------------------------------------------------------#

#This code is to plot the shifting distribution of VAFs over time - single run of simulation but storing information after each generation
muts_per_mitochondria_per_generation=5e-4
print(muts_per_mitochondria_per_generation)
mito_copy_number=600 #This is the
number_of_cells_in_simulation=1000 #Need enough 'cells' to have adequate mutation numbers to define the distribution
mutation_introductions_per_generation=muts_per_mitochondria_per_generation*number_of_cells_in_simulation*mito_copy_number

#Define the VAF 'bins' that will be used for comparing VAF distributions
new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")
VAF_groups=data.frame(
  labels=new_VAF_groups,
  lower_limit=c(0,2^(-10:-1)),
  upper_limit=2^(-10:0)
)

#The Wright-Fisher simulation below is stochastic (fisher_wright_drift), as is
#the binomial resampling of its output, so fix the seed here - before the loop -
#to make Fig 4b and Extended Data Fig 7a exactly reproducible.
set.seed(42)

#Record the distribution of VAFs every 10 generations, though the generation of mutation acquisition is recorded exactly
gens_to_include=seq(10,1500,10)
gens_record=vector(mode="list",length = length(gens_to_include))
vaf_by_gen_list=vector(mode="list")
last_gen=0
curr_vafs=c()
for(ngen in gens_to_include) {
  print(ngen)
  if(length(vaf_by_gen_list)>0) {
    drifted_old_muts=lapply(vaf_by_gen_list,function(gen_vafs) {
      new_gen_vafs<-sapply(gen_vafs,function(vaf) fisher_wright_drift(vaf,population_size = mito_copy_number,1,ngen-last_gen))
      return(new_gen_vafs[new_gen_vafs>0])
    })
  } else {
    drifted_old_muts<-NULL
  }
  
  new_muts_gen_list=vector(mode="list",length=ngen-last_gen)
  names(new_muts_gen_list)<-(last_gen+1):ngen
  for(j in 1:(ngen-last_gen)){
    gen_mut_vafs<-vector(length=mutation_introductions_per_generation)
    for(i in 1:mutation_introductions_per_generation) {
      this_mut_final_vaf<- fisher_wright_drift(1/mito_copy_number,population_size = mito_copy_number,1,j)
      gen_mut_vafs[i]<-this_mut_final_vaf
    }
    new_muts_gen_list[[j]]<-gen_mut_vafs[gen_mut_vafs>0]
  }
  vaf_by_gen_list<-c(drifted_old_muts,new_muts_gen_list)
  gen_summary=Map(vafs=vaf_by_gen_list,gen=names(vaf_by_gen_list),function(vafs,gen) if(length(vafs)>0){data.frame(gen=gen,vaf=vafs)}else {NULL})%>%dplyr::bind_rows()%>%mutate(total_gens=ngen)
  gens_record[[which(gens_to_include==ngen)]]<-gen_summary
  last_gen<-ngen
}

#Here, 7.5 is the 'haploid sequencing coverage' in this simulated experiment
#This approximates the ~15X diploid coverage of most of the experiments
haploid_coverage=7.5
gens_record<-dplyr::bind_rows(gens_record)
gens_record$observed_vaf=sapply(gens_record$vaf,function(vaf) rbinom(n=1,size=(haploid_coverage*mito_copy_number),prob=vaf)/(haploid_coverage*mito_copy_number))

gens_record$VAF_group=sapply(1:nrow(gens_record),function(i) {
  if(gens_record$observed_vaf[i]==0) {
    return("Absent")
  } else {
    return(VAF_groups$labels[VAF_groups$lower_limit<gens_record$observed_vaf[i] & VAF_groups$upper_limit>=gens_record$observed_vaf[i]])
  }
})

#-----------------------------------------------------------------------------------#
### Generate FIG. 4B ---------
#-----------------------------------------------------------------------------------#

#Visualize the distribution at a selection of "generations" to get a sense of the evolution of the distribution through time
modelled_VAF_dist<-gens_record%>%
  filter(total_gens%in%c(100,200,400,800,1200))%>%
  mutate(total_gens=factor(paste(total_gens,"generations"),levels=paste(c(100,200,400,800,1200),"generations")))%>%
  mutate(VAF_group=factor(VAF_group,levels=new_VAF_groups),gen=factor(gen,levels=1:1200))%>%
  ggplot(aes(x=VAF_group,fill=gen))+
  geom_bar(stat="count")+
  facet_wrap(~total_gens,ncol=5)+
  theme_bw()+
  theme(legend.position = "none",axis.text.x = element_text(size=5,angle=90),strip.text.x = element_text(size=6))+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8,"Spectral"))(1200))+
  labs(x="VAF group",y="Count")+
  my_theme

ggsave(filename=paste0(plots_dir,"Figure_04/Fig4b.Modelled_VAF_dist.pdf"),width = 6,height=2)

#-----------------------------------------------------------------------------------#
### Generate EXTENDED DATA FIG. 6A ---------
#-----------------------------------------------------------------------------------#

#One notable feature is that these simulations suggest that high VAF mutations were acquired early in life
#whereas low VAF mutations were invariably acquired very recently.
acquisition_time_by_VAF_ridges<-gens_record%>%
  mutate(VAF_group=factor(VAF_group,levels=VAF_groups$labels))%>%
  filter(total_gens==1250)%>%
  ggplot(aes(x=as.numeric(gen),y=VAF_group,fill = factor(stat(quantile),levels=1:4))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,quantiles = 4, quantile_lines = TRUE,scale=4) +
  scale_fill_viridis_d(name = "Quartiles")+
  scale_x_continuous(limits=c(0,1350),breaks=seq(0,1250,250))+
  theme_classic()+
  my_theme+
  labs(x="Time of mutation acquisition\n(WF Generations)",y="VAF level")

ggsave(filename = paste0(plots_dir,"Extended_Data_Figure_06/ExtDataFig6a.acquisition_time_by_VAF_ridges.pdf"),acquisition_time_by_VAF_ridges,width=3.3,height=2.5)

#-----------------------------------------------------------------------------------#
### Generate animation of evolving VAF distribution ---------
#-----------------------------------------------------------------------------------#

## Can also be informative to visualize the evolving distribution over time as an animation
gens_record_summary<-gens_record%>%
  group_by(total_gens,gen,VAF_group)%>%
  summarise(n=n())
gens_record_summary$VAF_group=factor(gens_record_summary$VAF_group,levels=VAF_groups$labels)
gens_record_summary$gen=as.numeric(gens_record_summary$gen)

library(gganimate)
drift_with_age_animation<-gens_record_summary%>%
  filter(!is.na(VAF_group))%>%
  arrange(gen)%>%
  ggplot(aes(x=VAF_group,y=n,fill=gen))+
  geom_bar(stat="identity")+
  transition_manual(frames=total_gens)+
  theme_bw()+
  theme(title = element_text(size=15),legend.position = "none",axis.text.x = element_text(angle=90,size=15),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12))+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(8,"Spectral"))+
  labs(title='{current_frame} generations',x="Observed VAF",y="Count")
#Render at an explicit size: gganimate's default (480x480) depends on the
#package version, so without this the animation silently changes dimensions
#between environments. 960x768 matches the committed gif.
anim_save(filename=paste0(rebuttal_figs_dir,"driftwithage.gif"),
          animation = animate(drift_with_age_animation, width = 960, height = 768,
                              renderer = gifski_renderer()))

#Plot the colour scale for these figures
max_gen=1500
lut=colorRampPalette(RColorBrewer::brewer.pal(8,"Spectral"))(max_gen)
scale = (length(lut)-1)
pdf(paste0(rebuttal_figs_dir,"mito_muts_by_gen_scale.pdf"),height=6,width = 3)
plot(c(0,10), c(0,1), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
axis(side=2,pos=1,las=1,at=seq(0,max_gen,200)/scale,labels=seq(0,max_gen,200))
for (i in 1:(length(lut)-1)) {
  y = (i-1)/scale
  rect(1,y,2,y+1/scale, col=lut[i], border=NA)
}
dev.off()

#-----------------------------------------------------------------------------------#
#----------------------PERFORM THE ABC----------------------
#-----------------------------------------------------------------------------------#

#This section uses summary statistics collected from many simulations and compares the results to those from the
#data statistics for each individual (as collected above).

#Get the summary stats for the data into a dataframe
sumstats.data<-left_join(real_muts_dist_df,nsamp_df)%>%
  mutate(muts_per_samp=abs_muts/n_samp)%>%
  dplyr::select(exp_ID,VAF_range,muts_per_samp)%>%
  pivot_wider(names_from="VAF_range",values_from="muts_per_samp")%>%
  replace(is.na(.), 0)

#The actual simulations are run in the separate script "Mito_VAF_distribution_simulations.R"
library(abc)
all_sumstats=readRDS(file=paste0(root_dir,"/data/Drift_ABC_VAF_distribution/VAF_distribution_ABC_simulation_sumstats_combined.Rds"))

new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")

#A VAF bin containing no mutations means zero mutations in that bin, but
#pivot_wider() in the simulation script simply omits the column, so bind_rows()
#filled it with NA. abc() then returns an NA distance for that simulation, which
#can never be among the closest - so those simulations are silently unreachable.
#The NAs are concentrated at low generation counts (98.7% of simulations under
#100 generations have at least one), which is exactly the region the youngest
#individuals need, so the reference table was effectively truncated from below.
#Restoring the zeros returns 234 usable simulations under 100 generations (was 3).
all_sumstats[,new_VAF_groups]<-lapply(all_sumstats[,new_VAF_groups],
                                      function(x) ifelse(is.na(x),0,x))

#Exclude the lowest VAF categories from abc that may be unreliably captured by sequencing
VAF_categories_to_include_in_ABC=new_VAF_groups[3:11]

mtDNA_CN=600
abc_res_df<-lapply(1:nrow(sumstats.data),function(i) {
  Exp_ID=sumstats.data[i,1]
  res<-abc(target = as.numeric(sumstats.data[i,VAF_categories_to_include_in_ABC]),
           param = all_sumstats[,1:2],
           sumstat = all_sumstats[,VAF_categories_to_include_in_ABC],
           tol = 0.02,
           transf = c("log","log"),
           method = "neuralnet")
  stats=apply(res$unadj.values,2,quantile,c(0.025,0.5,0.975))
  stats=as.data.frame(stats)%>%tibble::rownames_to_column(var="quantile")%>%pivot_wider(names_from = "quantile",values_from=c("muts_per_mitochondria_per_generation","total_generations"))
  return(cbind(data.frame(exp_ID=Exp_ID),stats))
})%>%dplyr::bind_rows()%>%
  left_join(ref_df,by=c("exp_ID"="Sample"))

#-----------------------------------------------------------------------------------#
### Generate FIG. 4E ---------
#-----------------------------------------------------------------------------------#

abc_res_plot_mut_rate<-abc_res_df%>%
  #Exclude only the foetal samples. The cord bloods are retained: once the
  #empty VAF bins are restored to zero above, the reference table resolves the
  #low generation counts they need, and their posteriors are no longer pinned
  #against the reachable floor.
  filter(!exp_ID%in%c("8pcw","18pcw"))%>%
  ggplot(aes(x=forcats::fct_reorder(exp_ID,Age),y=`muts_per_mitochondria_per_generation_50%`,
             ymin=`muts_per_mitochondria_per_generation_2.5%`,
             ymax=`muts_per_mitochondria_per_generation_97.5%`))+
  geom_point(alpha=0.75,size=0.5)+
  geom_errorbar(width=0.3,alpha=0.5)+
  theme_classic()+
  scale_y_continuous(limits=c(1e-4,7e-4))+
  labs(x="Individual",y="Mutations per\nmtDNA genome\nper generation")+
  my_theme+
  #coord_flip() swaps the panels but not the theme elements: axis.title.x still
  #refers to the bottom axis (the rate), axis.title.y to the donor axis. Blank
  #the latter, keeping labs(x=) set so no default expression can leak through.
  theme(axis.title.y = element_blank())+
  coord_flip()

ggsave(filename=paste0(plots_dir,"Figure_04/Fig4e.abc_res_plot_mut_rate.pdf"),abc_res_plot_mut_rate,width=1.4,height=2)

#-----------------------------------------------------------------------------------#
### Generate FIG. 4C ---------
#-----------------------------------------------------------------------------------#

##Now do linear mixed effects modelling using all values of the posterior distribution and using individual as a random effect
abc_res_df<-lapply(1:nrow(sumstats.data),function(i) {
  Exp_ID=sumstats.data[i,1]
  res<-abc(target = as.numeric(sumstats.data[i,VAF_categories_to_include_in_ABC]),
           param = all_sumstats[,1:2],
           sumstat = all_sumstats[,VAF_categories_to_include_in_ABC],
           tol = 0.05,
           transf = c("log","log"),
           method = "neuralnet")
  
  #Use the unadjusted (rejection) posterior, as in panel e. The neural-network
  #regression adjustment produced implausibly narrow posteriors for some
  #individuals (e.g. KX003) and extrapolated below the range of the reference
  #table for the cord bloods, so the rejection posterior is the more
  #conservative choice and matches the published Methods.
  return(as.data.frame(res$unadj.values)%>%mutate(exp_ID=Exp_ID$exp_ID))
})%>%dplyr::bind_rows()%>%
  left_join(ref_df%>%dplyr::select(Sample,Age),by=c("exp_ID"="Sample"))%>%
  #Exclude the foetal samples only. The cord bloods are retained: with the empty
  #VAF bins restored to zero above, the reference table resolves the low
  #generation counts they require, and they anchor the regression at age 0,
  #which roughly halves the width of the interval on the fitted drift rate.
  dplyr::filter(!exp_ID%in%c("8pcw","18pcw"))

#Perform the LMER
lme.gens_by_age<-lme4::lmer(total_generations~Age+(1|exp_ID),data=abc_res_df)
summary(lme.gens_by_age)
mtDNA_CN*365/lme.gens_by_age@beta[2]
mtDNA_CN*365/confint(lme.gens_by_age)["Age",]

#Plot these results
#95% CI of the fitted trend - see lmer_confidence_band() for what the band covers
gens_by_age_band<-lmer_confidence_band(lme.gens_by_age,abc_res_df,predictor="Age")

abc_res_plot_total_generations<-abc_res_df%>%
  dplyr::filter(!exp_ID%in%c("8pcw","18pcw"))%>%
  mutate(exp_ID=factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)]))%>%
  ggplot(aes(x=Age,y=total_generations))+
  geom_ribbon(data=gens_by_age_band,aes(x=Age,ymin=lower,ymax=upper),
              inherit.aes=FALSE,fill="grey50",alpha=0.25)+
  geom_point(aes(col=exp_ID),alpha=0.05,size=0.25)+
  scale_y_continuous(limits=c(0,1700))+
  geom_abline(slope=lme.gens_by_age@beta[2],intercept = lme.gens_by_age@beta[1],linetype=1)+
  scale_color_manual(values=Individual_cols[-c(1:2)])+
  theme_classic()+
  guides(colour=guide_legend(override.aes = list(alpha=1)))+
  labs(x="Age",y="Total WF generations\n(Posterior distribution from ABC)",col="")+
  my_theme+
  theme(legend.box.spacing = unit(0,"mm"),
        legend.key.size = unit(0.5,"mm"),
        legend.spacing = unit(0,"mm"),
        legend.box.margin = margin(c(0,0,0,0)),
        legend.text = element_text(margin = margin(t=0)))

ggsave(filename=paste0(plots_dir,"Figure_04/Fig4c.abc_res_plot_total_generations.pdf"),abc_res_plot_total_generations,width=2,height=2)

#-----------------------------------------------------------------------------------#
### Generate FIG. 4D ---------
#-----------------------------------------------------------------------------------#

drift_parameter_ml=mtDNA_CN*365/lme.gens_by_age@beta[2]
drift_parameter_CI<-mtDNA_CN*365/confint(lme.gens_by_age)["Age",]
drift_parameter_lowerCI=drift_parameter_CI[2]
drift_parameter_upperCI=drift_parameter_CI[1]

population_size_estimate_1=740
population_size_estimate_2=population_size_estimate_1/5
xrange=seq(50,10000,10)
#Annotation text is sized in mm, so it does not shrink with the panel. This
#panel is saved at 2x2in to match 4c (it was 3x3in), so the sizes and anchors
#below are scaled to suit: the day labels are left-aligned off the y-axis
#rather than centred on x=60, and the rotated "Pop. size" labels are anchored
#at their top (hjust=1) so they hang downwards inside the panel.
text_size=1.6

pop_vs_gentime_schematic<-data.frame(x=xrange,y=drift_parameter_ml/xrange,ymin=drift_parameter_lowerCI/xrange,ymax=drift_parameter_upperCI/xrange)%>%
  ggplot(aes(x=x,y=y,ymin=ymin,ymax=ymax))+
  geom_line()+
  annotate(geom="text",label="Combinations consistent\nwith drift rate",x=3300,y=9.5,col="black",angle=317,size=text_size,lineheight=0.9)+
  geom_ribbon(alpha=0.2)+
  geom_path(data=data.frame(x=c(0,population_size_estimate_1,population_size_estimate_1),y=c(drift_parameter_ml/population_size_estimate_1,drift_parameter_ml/population_size_estimate_1,0)),aes(x=x,y=y),col="red",linetype=2,inherit.aes = F)+
  annotate(geom="text",label="Pop. size = mtDNA CN",x=population_size_estimate_1+130,y=13,col="red",angle=270,hjust=1,size=text_size)+
  annotate(geom="text",label=paste0(round(drift_parameter_ml/population_size_estimate_1)," days"),x=52,y=(drift_parameter_ml/population_size_estimate_1)*1.25,col="red",hjust=0,size=text_size)+
  geom_path(data=data.frame(x=c(0,population_size_estimate_2,population_size_estimate_2),y=c(drift_parameter_ml/population_size_estimate_2,drift_parameter_ml/population_size_estimate_2,0)),aes(x=x,y=y),col="blue",linetype=2,inherit.aes = F)+
  annotate(geom="text",label="Pop. size = # of mitochondria",x=population_size_estimate_2+26,y=13,col="blue",angle=270,hjust=1,size=text_size)+
  annotate(geom="text",label=paste0(round(drift_parameter_ml/population_size_estimate_2)," days"),x=52,y=(drift_parameter_ml/population_size_estimate_2)*1.25,col="blue",hjust=0,size=text_size)+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  my_theme+
  labs(x="Effective mitochondrial\npopulation size",y="Generation time (days)")

ggsave(filename = paste0(plots_dir,"Figure_04/Fig4d.pop_vs_gentime_schematic.pdf"),pop_vs_gentime_schematic,width=2,height=2)


#-----------------------------------------------------------------------------------#
### Generate FIG. 4F ---------
#-----------------------------------------------------------------------------------#
# Drift during embryogenesis, inferred from the heteroplasmic oocyte mutations
# in the 8 pcw foetus. The per-mutation simulations were run in
# full_analysis_scripts/Heteroplasmic_oocyte_mutation_analysis.R and saved to
# data/Drift_ABC_foetal; this section re-runs only the inference from them.
#
# NB this is a single joint ABC: the summary statistics of all six mutations are
# concatenated into one target vector, fitted against one shared generation-time
# parameter. It is not a product of per-mutation posteriors - the per-mutation
# abc.nn.file.*.RDS outputs are only a check for outliers.
# full_analysis_scripts/Drift_ABC2_analysis.R does combine them by multiplying
# binned posteriors, but that script is orphaned: it reads from a directory that
# does not exist (data/Drift_ABC2) and log10()s a column that is already
# log_generation_time. It is superseded by this section.

het_oocyte_muts<-readRDS(paste0(root_dir,"/data/het_oocyte_muts.Rds"))
foetal_muts<-het_oocyte_muts[["8pcw"]]

#Use the rejection posterior, as in panels c and e. The originally reported
#value came from the neural-network regression-adjusted posterior, which gives a
#roughly five-fold narrower interval here - the same over-tightening seen in
#panel e - so the rejection posterior is the more conservative choice.
foetal_posterior_type<-"unadj.values"

#Observed summary statistics per mutation - the ABC target, recomputed from the
#VAF matrix exactly as in the source script.
foetal_sumstat_fun<-function(v) c(mean(v),min(v),max(v),
                                  quantile(v,0.05),quantile(v,0.5),quantile(v,0.95),sd(v))
foetal_tips<-mito_data$`8pcw`$tree$tip.label
all_targets<-lapply(foetal_muts,function(mut)
  foetal_sumstat_fun(as.numeric(mito_data$`8pcw`$matrices$vaf[mut,foetal_tips])))

all_params_and_sumstats<-lapply(foetal_muts,function(mut)
  readRDS(paste0(root_dir,"/data/Drift_ABC_foetal/params_and_sumstats_",mut,".RDS")))

#Each mutation contributes its own block of summary statistics; the generation
#time parameter is shared, so the simulations are matched across mutations.
all_sumstats_mat<-Reduce(cbind,Map(list=all_params_and_sumstats,mut=foetal_muts,
  f=function(list,mut) {
    colnames(list$sumstats)<-paste(mut,colnames(list$sumstats),sep="_")
    list$sumstats
  }))
all_params<-all_params_and_sumstats[[1]]$params[,"log_generation_time",drop=FALSE]
all_targets_vec<-Reduce(c,all_targets)

abc.comb.nn<-abc::abc(target=all_targets_vec,param=all_params,sumstat=all_sumstats_mat,
                      tol=0.1,transf=c("none"),method="neuralnet")

foetal_gentime_posterior<-10^quantile(abc.comb.nn[[foetal_posterior_type]],c(0.025,0.5,0.975))

#Population size is fixed within these simulations, so the drift parameter -
#population size x generation time - is directly comparable with panel d.
foetal_pop_size<-unique(all_params_and_sumstats[[1]]$params[,"population_size"])
stopifnot(length(foetal_pop_size)==1)
foetal_drift_parameter<-foetal_pop_size*foetal_gentime_posterior

cat("\nFig 4f - embryonic drift (",foetal_posterior_type,", ",length(foetal_muts),
    " oocyte mutations, population size ",foetal_pop_size,")\n",sep="")
cat("  generation time: ",signif(foetal_gentime_posterior[2],3)," days (95% PI ",
    signif(foetal_gentime_posterior[1],3),"-",signif(foetal_gentime_posterior[3],3),")\n",sep="")
cat("  drift parameter: ",signif(foetal_drift_parameter[2],3)," mitochondria days (95% PI ",
    signif(foetal_drift_parameter[1],3),"-",signif(foetal_drift_parameter[3],3),")\n",sep="")
cat("  fold faster than the adult estimate: ",
    signif(drift_parameter_ml/foetal_drift_parameter[2],3),"x\n",sep="")

foetal_gentime_posterior_plot<-data.frame(type="posterior",log_gen_time=abc.comb.nn[[foetal_posterior_type]])%>%
  rbind(data.frame(type="prior",log_gen_time=as.numeric(all_params)))%>%
  mutate(gen_time=10^log_gen_time)%>%
  ggplot(aes(x=gen_time,col=type))+
  stat_density(geom="line",position="identity")+
  labs(x="WF generation time (days)",y="Density",col="")+
  theme_bw()+
  scale_x_log10(limits=c(0.1,500))+
  #Headroom above the posterior peak so the summary values can be printed at the
  #top without overlapping the curve. coord_cartesian() rather than a scale
  #limit, so the density is never censored if the peak grows.
  coord_cartesian(ylim=c(0,3))+
  my_theme+
  theme(legend.position="top",legend.box.spacing=unit(0,"mm"))

ggsave(filename = paste0(plots_dir,"Figure_04/Fig4f.foetal_gentime_posterior.pdf"),foetal_gentime_posterior_plot,width=1.75,height=2)
