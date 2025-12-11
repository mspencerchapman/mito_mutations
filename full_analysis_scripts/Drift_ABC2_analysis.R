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
R_scripts_dir = "~/R_work/my_functions/"
treemut_dir="~/R_work/treemut/"
R_scripts=list.files(R_scripts_dir,pattern = ".R",full.names = T)
sapply(R_scripts[-2],source)

#Set the key file paths using the root dir
tree_file_paths = list.files(paste0(root_dir,"data/tree_files"),pattern=".tree",full.names = T)
ref_file=paste0(root_dir,"data/Samples_metadata_ref.csv")
figures_dir=paste0(root_dir,"figures/")

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


#Pulling together the abc files and multiply together the posteriors
library(abc)
library(dplyr)
library(ggplot2)
my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"),
                strip.text.x = element_text(size=6))

abc_output_files<-list.files(path = paste0(root_dir,"data/Drift_ABC2"),pattern="abc.file.",full.names = F)
abc_values<-lapply(abc_output_files,function(file){
  abc.nn<-readRDS(paste0(root_dir,"data/Drift_ABC2/",file))
  exp_ID<-stringr::str_split(file,pattern="\\.",simplify = T)[,3]
  mut<-stringr::str_split(file,pattern="\\.",simplify = T)[,4]
  return(data.frame(exp_ID=exp_ID,mut=mut,adj.values=log10(abc.nn$adj.values[,3]),unadj.values=log10(abc.nn$unadj.values[,3])))
})

#Group posterior results into "bins" of width 0.1 (on a log scale)
bin_width=0.1
bins=seq(-1,2.7-bin_width,bin_width)
all_post<-lapply(abc_values,function(df){
  n_per_bin<-sapply(bins,function(ll) sum(df$adj.values>ll&df$adj.values<ll+bin_width))
  names(n_per_bin)<-bins
  return(n_per_bin)
})%>%
  dplyr::bind_rows()

#Now multiply these posterior distributions across all mutations & normalize
#This gives a final posterior distribution for the drift that it most consistent with ALL individual mutations
comb_post<-apply(all_post,2,prod) #Multiply
comb_post<-comb_post/sum(comb_post) #Normalize

#Now plot results
#First plot the posteriors for each individual mutations
individual_mutation_posteriors<-abc_values%>%
  dplyr::bind_rows()%>%
  mutate(generation_time=10^unadj.values)%>%
  ggplot(aes(x=generation_time))+
  geom_histogram()+
  facet_wrap(~mut,nrow = 2)+
  theme_bw()+
  scale_x_log10(limits=c(0.1,500))+
  labs(x="Generation time (days)",y="Posterior distribution")+
  my_theme

ggsave(filename = paste0(figures_dir,"Supp_Figure_07/Drift_ABC2_individual_posteriors.pdf"),individual_mutation_posteriors,width = 7,height=2.5)

combined_posterior<-data.frame(log_generation_time=as.numeric(names(comb_post)),
                  posterior_probability=comb_post)%>%
  mutate(generation_time=10^(log_generation_time+0.5*bin_width))%>%
  ggplot(aes(x=generation_time,
             y=posterior_probability))+
  geom_bar(col="black",fill="lightblue",stat="identity")+
  geom_point(size=0.3)+
  scale_x_log10()+
  geom_line(size=0.25)+
  theme_bw()+
  labs(x="Generation time (days)",y="Posterior distribution")+
  my_theme

ggsave(filename = paste0(figures_dir,"Supp_Figure_07/Drift_ABC2_combined_posterior.pdf"),combined_posterior,width = 2,height=1.5)

#This cumulative probability plot helps to estimate the 95% posterior interval
posterior_probability_cumsum<-data.frame(log_generation_time=as.numeric(names(comb_post)),
           posterior_probability=comb_post)%>%
  mutate(log_generation_time=log_generation_time+0.5*bin_width)%>%
  mutate(cumsum=cumsum(posterior_probability))%>%
  ggplot(aes(x=log_generation_time,y=cumsum))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = c(0.025,0.975),linetype=2)+
  scale_x_continuous(breaks=seq(-1,2.7,by=0.05))+
  labs(x="Log10 effective generation time (days)",y="Cumulative sum of posterior probability")+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90))

#Using this plot estimate posterior intervals as:
#1.05 - 1.45, max posterior density = 1.25

ggsave(filename = paste0(figures_dir,"Additional_plots/Drift_ABC2_posterior_probability_cumsum.pdf"),posterior_probability_cumsum,width = 2,height=1.5)





