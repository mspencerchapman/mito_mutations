## This is a simple simulation model for drift of the overall heteroplasmy level during clonal expansion from a single cell
#1. A 'starting_vaf' (between 0 and 1) represents the initial heteroplasmy level in the founder cell
#2. During clonal expansion, cell division is imagined to proceed in a co-ordinated fashion for a specified number of generations ('ngen'), such that in each generation, each cell divides once
#3. During each cell division, the mtDNA copy number doubles, with new mtDNA molecules being randomly made from the original mtDNA population according to a binomial distribution.
#4. The mtDNA molecules are then randomly partitioned between the two daughter cells

# setwd("/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/colony_drift_simulation")
# save_dir<-"model_data/"
# plots_dir<-"plots"

setwd("~/R_work/mito_mutations_blood/")
save_dir<-"data/colony_drift_simulations/"
plots_dir<-"rebuttal_plots/"
CORES=5

library(dplyr)
library(ggplot2)
library(parallel)

#SET BASIC PARAMETERS FOR ALL THE MODELS
mito_CN=675
ngen=12 #(clone will grow to 2^12 cells ~4,100)
nsim=200
vafs_for_testing<-c(0.05,0.2,0.5,0.8,0.95)

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

#############################
# Record the actual individual cell heteroplasmy levels
############################
print("STARTING MODEL 1")
model1_path<-paste0(save_dir,"model1_res.Rds")
if(file.exists(model1_path)) {
  model1_res<-readRDS(model1_path)
} else {
  model1_res<-mclapply(X=vafs_for_testing,mc.cores=CORES,FUN=function(starting_vaf) {
    print(starting_vaf)
    final_colony_vaf_list<-lapply(1:nsim,function(i) {
      if(i%%10==0) {cat(i,sep="\n")}
      cell_vec<-starting_vaf
      for(gen in 1:ngen) {
        
        cell_vec<-unlist(lapply(cell_vec, function(vaf) {
          
          if(vaf==0) {stop(return(c(0,0)))}
          
          #existing mtDNA molecules are chosen at random to replicate during doubling of the mtDNA copy number
          new_mut_mitos<-rbinom(1,size=mito_CN,prob=vaf)
          total_mut_mito<-new_mut_mitos+round(vaf*mito_CN)
          
          #Mutant mitoDNA molecules are randomly partitioned between daughter cells
          #Sampling is based on the idea of hypergeometric distribution as mutant mtDNA is sampled WITH REPLACEMENT
          daughter1_mut_mito<-rhyper(nn=1,m=total_mut_mito,n=(2*mito_CN)-total_mut_mito,k=mito_CN)
          daughter2_mut_mito=total_mut_mito-daughter1_mut_mito
          c(daughter1_mut_mito/mito_CN,daughter2_mut_mito/mito_CN)
        }))
      }
      return(cell_vec)
    })
    
    final_colony_vaf_mat<-Reduce(rbind,final_colony_vaf_list)
    return(final_colony_vaf_mat)
  })
  saveRDS(object=model1_res,file=model1_path)
}

## MODEL 2
## As MODEL 1, but now, mtDNA copy number doubling happens according to a birth process
print("STARTING MODEL 2")
model2_path<-paste0(save_dir,"model2_res.Rds")
if(file.exists(model2_path)) {
  model2_res<-readRDS(model2_path)
} else {
  model2_res<-mclapply(X=vafs_for_testing,mc.cores=CORES,FUN=function(starting_vaf) {
    print(starting_vaf)
    final_colony_vaf_list<-lapply(1:nsim,function(i) {
      if(i%%5==0) {cat(i,sep="\n")}
      cell_vec<-starting_vaf
      for(gen in 1:ngen) {
        
        cell_vec<-unlist(lapply(cell_vec, function(vaf) {
          
          if(vaf==0) {stop(return(c(0,0)))}
          
          #existing mtDNA molecules are chosen at random to replicate during doubling of the mtDNA copy number
          initial_mut_mito=round(vaf*mito_CN)
          new_vaf<-vaf
          
          #This time, grow the mitoCN by a birth process i.e.
          for(j in 1:mito_CN) {
            new_mol=rbinom(1,1,prob = new_vaf)
            curr_mito_CN=mito_CN+j
            new_vaf<-(round(new_vaf*(curr_mito_CN-1)) + new_mol)/curr_mito_CN
          }
          total_mut_mito<-round(curr_mito_CN*new_vaf)
          
          #Mutant mitoDNA molecules are randomly partitioned between daughter cells
          #Sampling is based on the idea of hypergeometric distribution as mutant mtDNA is sampled WITH REPLACEMENT
          daughter1_mut_mito<-rhyper(nn=1,m=total_mut_mito,n=curr_mito_CN-total_mut_mito,k=mito_CN)
          daughter2_mut_mito=total_mut_mito-daughter1_mut_mito
          c(daughter1_mut_mito/mito_CN,daughter2_mut_mito/mito_CN)
        }))
      }
      return(cell_vec)
    })
    
    final_colony_vaf_mat<-Reduce(rbind,final_colony_vaf_list)
    return(final_colony_vaf_mat)
  })
  
  saveRDS(object=model2_res,file=model2_path)
}



## MODEL 3
## As MODEL 2, but now cell division at each generation happens by a fisher-wright process, i.e. individual cells may divide multiple times, and others not at all

print("STARTING MODEL 3")
model3_path<-paste0(save_dir,"model3_res.Rds")
if(file.exists(model3_path)) {
  model3_res<-readRDS(model3_path)
} else {
  model3_res<-mclapply(X=vafs_for_testing,mc.cores=CORES,FUN=function(starting_vaf) {
    print(starting_vaf)
    final_colony_vaf_list<-lapply(1:nsim,function(i) {
      if(i%%10==0) {cat(i,sep="\n")}
      cell_vec<-starting_vaf
      for(gen in 1:ngen) {
        
        #Add in a sampling with replacement step (this is the Wright-Fisher drift step)
        cell_vec<-sample(cell_vec,replace = T)
        
        #Proceed as before
        cell_vec<-unlist(lapply(cell_vec, function(vaf) {
          
          if(vaf==0) {stop(return(c(0,0)))}
          
          #existing mtDNA molecules are chosen at random to replicate during doubling of the mtDNA copy number
          initial_mut_mito=round(vaf*mito_CN)
          new_vaf<-vaf
          
          #This time, grow the mitoCN by a birth process i.e.
          for(j in 1:mito_CN) {
            new_mol=rbinom(1,1,prob = new_vaf)
            curr_mito_CN=mito_CN+j
            new_vaf<-(round(new_vaf*(curr_mito_CN-1)) + new_mol)/curr_mito_CN
          }
          total_mut_mito<-round(curr_mito_CN*new_vaf)
          
          #Mutant mitoDNA molecules are randomly partitioned between daughter cells
          #Sampling is based on the idea of hypergeometric distribution as mutant mtDNA is sampled WITH REPLACEMENT
          daughter1_mut_mito<-rhyper(nn=1,m=total_mut_mito,n=curr_mito_CN-total_mut_mito,k=mito_CN)
          daughter2_mut_mito=total_mut_mito-daughter1_mut_mito
          c(daughter1_mut_mito/mito_CN,daughter2_mut_mito/mito_CN)
        }))
      }
      return(cell_vec)
    })
    
    final_colony_vaf_mat<-Reduce(rbind,final_colony_vaf_list)
    return(final_colony_vaf_mat)
  })
  
  saveRDS(object=model3_res,file=model3_path)
}



## MODEL 4
## As MODEL 1, but now cell division at each generation happens by a fisher-wright process, i.e. individual cells may divide multiple times, and others not at all
print("STARTING MODEL 4")
model4_path<-paste0(save_dir,"model4_res.Rds")
if(file.exists(model4_path)) {
  model4_res<-readRDS(model4_path)
} else {
  model4_res<-mclapply(X=vafs_for_testing,mc.cores=CORES,FUN=function(starting_vaf) {
    print(starting_vaf)
    final_colony_vaf_list<-lapply(1:nsim,function(i) {
      if(i%%10==0) {cat(i,sep="\n")}
      cell_vec<-starting_vaf
      for(gen in 1:ngen) {
        
        #Add in a sampling with replacement step
        cell_vec<-sample(cell_vec,replace = T)
        
        cell_vec<-unlist(lapply(cell_vec, function(vaf) {
          
          if(vaf==0) {stop(return(c(0,0)))}
          
          #existing mtDNA molecules are chosen at random to replicate during doubling of the mtDNA copy number
          new_mut_mitos<-rbinom(1,size=mito_CN,prob=vaf)
          total_mut_mito<-new_mut_mitos+round(vaf*mito_CN)
          
          #Mutant mitoDNA molecules are randomly partitioned between daughter cells
          #Sampling is based on the idea of hypergeometric distribution as mutant mtDNA is sampled WITH REPLACEMENT
          daughter1_mut_mito<-rhyper(nn=1,m=total_mut_mito,n=(2*mito_CN)-total_mut_mito,k=mito_CN)
          daughter2_mut_mito=total_mut_mito-daughter1_mut_mito
          c(daughter1_mut_mito/mito_CN,daughter2_mut_mito/mito_CN)
        }))
      }
      return(cell_vec)
    })
    
    final_colony_vaf_mat<-Reduce(rbind,final_colony_vaf_list)
    return(final_colony_vaf_mat)
  })
  saveRDS(object=model4_res,file=model4_path)
}


################################
# Generate individual plots
################################

plot1<-Map(vaf_level=vafs_for_testing,cell_vafs=model1_res,function(vaf_level,cell_vafs) {
  mean_vafs=apply(cell_vafs,1,mean)
  data.frame(model="model1",vaf_level=vaf_level,stat="Mean heteroplasmy\nacross colony",res=mean_vafs)%>%
    dplyr::bind_rows(data.frame(model="model1",vaf_level=vaf_level,stat="Individual cell\nheteroplasmy distribution",res=sample(as.numeric(cell_vafs),size=4*nsim,replace=F)))
})%>%dplyr::bind_rows()%>%
  ggplot(aes(x=res,fill = stat))+
  geom_density(alpha=0.25)+
  geom_vline(aes(xintercept=vaf),linetype=2,data = data.frame(vaf_level=vafs_for_testing,vaf=vafs_for_testing))+
  facet_grid(rows=vars(vaf_level))+
  theme_classic()+
  my_theme+
  theme(panel.border = element_rect(fill = NA),strip.text.y=element_text(size=7,angle=0),legend.position="none")+
  labs(x="Mutant mtDNA heteroplasmy",y="Density",fill="")

plot2<-Map(vaf_level=vafs_for_testing,cell_vafs=model2_res,function(vaf_level,cell_vafs) {
  mean_vafs=apply(cell_vafs,1,mean)
  data.frame(model="model2",vaf_level=vaf_level,stat="Mean heteroplasmy\nacross colony",res=mean_vafs)%>%
    dplyr::bind_rows(data.frame(model="model2",vaf_level=vaf_level,stat="Individual cell\nheteroplasmy distribution",res=sample(as.numeric(cell_vafs),size=4*nsim,replace=F)))
})%>%dplyr::bind_rows()%>%
  ggplot(aes(x=res,fill = stat))+
  geom_density(alpha=0.25)+
  geom_vline(aes(xintercept=vaf),linetype=2,data = data.frame(vaf_level=vafs_for_testing,vaf=vafs_for_testing))+
  facet_grid(rows=vars(vaf_level))+
  theme_classic()+
  my_theme+
  theme(panel.border = element_rect(fill = NA),strip.text.y=element_text(size=7,angle=0),legend.position="none")+
  labs(x="Mutant mtDNA heteroplasmy",y="Density",fill="")


plot3<-Map(vaf_level=vafs_for_testing,cell_vafs=model4_res,function(vaf_level,cell_vafs) {
  mean_vafs=apply(cell_vafs,1,mean)
  data.frame(model="model4",vaf_level=vaf_level,stat="Mean heteroplasmy\nacross colony",res=mean_vafs)%>%
    dplyr::bind_rows(data.frame(model="model4",vaf_level=vaf_level,stat="Individual cell\nheteroplasmy distribution",res=sample(as.numeric(cell_vafs),size=4*nsim,replace=F)))
})%>%dplyr::bind_rows()%>%
  ggplot(aes(x=res,fill = stat))+
  #geom_histogram(binwidth=0.01)+
  geom_density(alpha=0.25)+
  geom_vline(aes(xintercept=vaf),linetype=2,data = data.frame(vaf_level=vafs_for_testing,vaf=vafs_for_testing))+
  facet_grid(rows=vars(vaf_level))+
  theme_classic()+
  my_theme+
  theme(panel.border = element_rect(fill = NA),strip.text.y=element_text(size=7,angle=0),legend.position="none")+
  labs(x="Mutant mtDNA heteroplasmy",y="Density",fill="")

plot4<-Map(vaf_level=vafs_for_testing,cell_vafs=model3_res,function(vaf_level,cell_vafs) {
  mean_vafs=apply(cell_vafs,1,mean)
  data.frame(model="model3",vaf_level=vaf_level,stat="Mean heteroplasmy\nacross colony",res=mean_vafs)%>%
    dplyr::bind_rows(data.frame(model="model3",vaf_level=vaf_level,stat="Individual cell\nheteroplasmy distribution",res=sample(as.numeric(cell_vafs),size=4*nsim,replace=F)))
})%>%dplyr::bind_rows()%>%
  ggplot(aes(x=res,fill = stat))+
  #geom_histogram(binwidth=0.01)+
  geom_density(alpha=0.25)+
  geom_vline(aes(xintercept=vaf),linetype=2,data = data.frame(vaf_level=vafs_for_testing,vaf=vafs_for_testing))+
  facet_grid(rows=vars(vaf_level))+
  theme_classic()+
  my_theme+
  theme(panel.border = element_rect(fill = NA),strip.text.y=element_text(size=7,angle=0),legend.position="right")+
  labs(x="Mutant mtDNA heteroplasmy",y="Density",fill="")

################################
# Combined visualization
################################
comb_plot<-Map(model=list(model1=model1_res,model2=model2_res,model3=model4_res,model4=model3_res),name=paste("Model",1:4),function(model,name) {
  Map(vaf_level=vafs_for_testing,cell_vafs=model,function(vaf_level,cell_vafs) {
    mean_vafs=apply(cell_vafs,1,mean)
    data.frame(vaf_level=vaf_level,stat="Mean heteroplasmy\nacross colony",res=mean_vafs)%>%
      dplyr::bind_rows(data.frame(vaf_level=vaf_level,stat="Individual cell\nheteroplasmy distribution",res=sample(as.numeric(cell_vafs),size=4*nsim,replace=F)))
  })%>%dplyr::bind_rows()%>%
    mutate(model=name,.before=1)
})%>%dplyr::bind_rows()%>%
  mutate(vaf_level=factor(paste0(100*vaf_level,"%"),levels=paste0(100*vafs_for_testing,"%")),res=100*res)%>%
  ggplot(aes(x=res,fill = stat))+
  geom_density(alpha=0.25,size=0.2)+
  geom_vline(aes(xintercept=vaf),linetype=2,linewidth=0.5,data = data.frame(vaf_level=factor(paste0(100*vafs_for_testing,"%"),levels=paste0(100*vafs_for_testing,"%")),vaf=100*vafs_for_testing))+
  facet_grid(rows=vars(vaf_level),cols=vars(model))+
  theme_classic()+
  my_theme+
  theme(panel.border = element_rect(fill = NA),
        strip.text.y=element_text(size=8,angle=0),
        legend.position="top",
        legend.text=element_text(size=8))+
  labs(x="Mutant mtDNA heteroplasmy (%)",y="Density",fill="")

ggsave(filename=paste0(plots_dir,"colony_drift_plots.pdf"),plot = comb_plot,width = 7,height=4)


model_summary_table<-Map(model=list(model1=model1_res,model2=model2_res,model3=model4_res,model4=model3_res),name=paste("Model",1:4),function(model,name) {
  Map(vaf=vafs_for_testing,res_mat=model,function(vaf,res_mat) {
    mean_colony_vafs<-apply(res_mat,1,mean)
    
    as.data.frame(t(quantile(mean_colony_vafs,c(0.025,0.5,0.975))))%>%
      dplyr::mutate(initial_cell_vaf=vaf,.before=1)
  })%>%dplyr::bind_rows()%>%
    mutate(model=name,.before=1)
})%>%dplyr::bind_rows()

model_summary_plot<-model_summary_table%>%
  mutate_if(is.numeric,function(x) {100*x})%>%
  ggplot(aes(x=initial_cell_vaf,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+
  geom_abline(intercept=0,slope=1,linetype = 2,col="red")+
  geom_point(size=0.9)+
  geom_errorbar(width=2)+
  theme_classic()+
  scale_x_continuous(limits=c(0,100))+
  scale_y_continuous(limits=c(0,100))+
  facet_grid(~model)+
  my_theme+
  labs(x="True mtDNA heteroplasmy\nin founder cell (%)",y="Mean mtDNA heteroplasmy\nacross colony cells (%)")

model_summary_table%>%filter(initial_cell_vaf==0.5)
model_summary_table%>%filter(initial_cell_vaf==0.05)

ggsave(filename=paste0(plots_dir,"model_summary_plot.pdf"),plot = model_summary_plot,width = 7,height=2.7)


####RECORDING ONLY THE 95% CONFIDENCE INTERVALS

# model1_res<-lapply(vafs_for_testing,function(starting_vaf) {
#   print(starting_vaf)
#   final_colony_vaf<-sapply(1:nsim,function(i) {
#     if(i%%10==0) {cat(i,sep="\n")}
#     cell_vec<-starting_vaf
#     for(gen in 1:ngen) {
#       
#       cell_vec<-unlist(lapply(cell_vec, function(vaf) {
#         
#         if(vaf==0) {stop(return(c(0,0)))}
#         
#         #existing mtDNA molecules are chosen at random to replicate during doubling of the mtDNA copy number
#         new_mut_mitos<-rbinom(1,size=mito_CN,prob=vaf)
#         total_mut_mito<-new_mut_mitos+round(vaf*mito_CN)
#         
#         #Mutant mitoDNA molecules are randomly partitioned between daughter cells
#         #Sampling is based on the idea of hypergeometric distribution as mutant mtDNA is sampled WITH REPLACEMENT
#         daughter1_mut_mito<-rhyper(nn=1,m=total_mut_mito,n=(2*mito_CN)-total_mut_mito,k=mito_CN)
#         daughter2_mut_mito=total_mut_mito-daughter1_mut_mito
#         c(daughter1_mut_mito/mito_CN,daughter2_mut_mito/mito_CN)
#       }))
#     }
#     return(mean(cell_vec))
#   })
#   
#   as.data.frame(t(quantile(final_colony_vaf,c(0.025,0.5,0.975))))%>%
#     dplyr::mutate(initial_cell_vaf=starting_vaf,.before=1)
# })%>%dplyr::bind_rows()
# 
# model1_res%>%
#   ggplot(aes(x=initial_cell_vaf,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+
#   geom_point()+
#   geom_errorbar(width=0.02)+
#   theme_classic()+
#   scale_x_continuous(limits=c(0,1))+
#   scale_y_continuous(limits=c(0,1))+
#   geom_abline(intercept=0,slope=1,linetype=2)+
#   theme()
# 
# ## MODEL 2
# ## As MODEL 1, but now, mtDNA copy number doubling happens according to a birth process
# 
# model2_res<-lapply(vafs_for_testing,function(starting_vaf) {
#   print(starting_vaf)
#   final_colony_vaf<-sapply(1:nsim,function(i) {
#     if(i%%10==0) {cat(i,sep="\n")}
#     cell_vec<-starting_vaf
#     for(gen in 1:ngen) {
#       
#       cell_vec<-unlist(lapply(cell_vec, function(vaf) {
#         
#         if(vaf==0) {stop(return(c(0,0)))}
#         
#         #existing mtDNA molecules are chosen at random to replicate during doubling of the mtDNA copy number
#         initial_mut_mito=round(vaf*mito_CN)
#         new_vaf<-vaf
#         
#         #This time, grow the mitoCN by a birth process i.e.
#         for(j in 1:mito_CN) {
#           new_mol=rbinom(1,1,prob = new_vaf)
#           curr_mito_CN=mito_CN+j
#           new_vaf<-(round(new_vaf*(curr_mito_CN-1)) + new_mol)/curr_mito_CN
#         }
#         total_mut_mito<-round(curr_mito_CN*new_vaf)
#         
#         #Mutant mitoDNA molecules are randomly partitioned between daughter cells
#         #Sampling is based on the idea of hypergeometric distribution as mutant mtDNA is sampled WITH REPLACEMENT
#         daughter1_mut_mito<-rhyper(nn=1,m=total_mut_mito,n=curr_mito_CN-total_mut_mito,k=mito_CN)
#         daughter2_mut_mito=total_mut_mito-daughter1_mut_mito
#         c(daughter1_mut_mito/mito_CN,daughter2_mut_mito/mito_CN)
#       }))
#     }
#     return(mean(cell_vec))
#   })
#   
#   as.data.frame(t(quantile(final_colony_vaf,c(0.025,0.5,0.975))))%>%
#     dplyr::mutate(initial_cell_vaf=starting_vaf,.before=1)
# })%>%dplyr::bind_rows()
# 
# model2_res%>%
#   ggplot(aes(x=initial_cell_vaf,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+
#   geom_point()+
#   geom_errorbar(width=0)+
#   theme_classic()+
#   scale_x_continuous(limits=c(0,1))+
#   scale_y_continuous(limits=c(0,1))+
#   geom_abline(intercept=0,slope=1,linetype=2)+
#   theme()
# 
# saveRDS(object=list(model1=model1_res,model2=model2_res),"colony_drift.Rds")
# 
# 
# ## MODEL 3
# ## As MODEL 2, but now cell division at each generation happens by a fisher-wright process, i.e. individual cells may divide multiple times, and others not at all
# 
# mito_CN=675
# ngen=14
# nsim=200
# 
# model3_res<-lapply(vafs_for_testing,function(starting_vaf) {
#   print(starting_vaf)
#   final_colony_vaf<-sapply(1:nsim,function(i) {
#     if(i%%10==0) {cat(i,sep="\n")}
#     cell_vec<-starting_vaf
#     for(gen in 1:ngen) {
#       
#       #Add in a sampling with replacement step
#       cell_vec<-sample(cell_vec,replace = T)
#       
#       #Proceed as before
#       cell_vec<-unlist(lapply(cell_vec, function(vaf) {
#         
#         if(vaf==0) {stop(return(c(0,0)))}
#         
#         #existing mtDNA molecules are chosen at random to replicate during doubling of the mtDNA copy number
#         initial_mut_mito=round(vaf*mito_CN)
#         new_vaf<-vaf
#         
#         #This time, grow the mitoCN by a birth process i.e.
#         for(j in 1:mito_CN) {
#           new_mol=rbinom(1,1,prob = new_vaf)
#           curr_mito_CN=mito_CN+j
#           new_vaf<-(round(new_vaf*(curr_mito_CN-1)) + new_mol)/curr_mito_CN
#         }
#         total_mut_mito<-round(curr_mito_CN*new_vaf)
#         
#         #Mutant mitoDNA molecules are randomly partitioned between daughter cells
#         #Sampling is based on the idea of hypergeometric distribution as mutant mtDNA is sampled WITH REPLACEMENT
#         daughter1_mut_mito<-rhyper(nn=1,m=total_mut_mito,n=curr_mito_CN-total_mut_mito,k=mito_CN)
#         daughter2_mut_mito=total_mut_mito-daughter1_mut_mito
#         c(daughter1_mut_mito/mito_CN,daughter2_mut_mito/mito_CN)
#       }))
#     }
#     return(mean(cell_vec))
#   })
#   
#   as.data.frame(t(quantile(final_colony_vaf,c(0.025,0.5,0.975))))%>%
#     dplyr::mutate(initial_cell_vaf=starting_vaf,.before=1)
# })%>%dplyr::bind_rows()
# 
# model3_res%>%
#   ggplot(aes(x=initial_cell_vaf,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+
#   geom_point()+
#   geom_errorbar()+
#   theme_classic()+
#   scale_x_continuous(limits=c(0,1))+
#   scale_y_continuous(limits=c(0,1))+
#   geom_abline(intercept=0,slope=1)+
#   theme()
# 
# 
# ## MODEL 4
# ## As MODEL 1, but now cell division at each generation happens by a fisher-wright process, i.e. individual cells may divide multiple times, and others not at all
# 
# mito_CN=675
# ngen=14
# nsim=200
# 
# 
# for(gen in 1:ngen) {
#   
#   #Add in a sampling with replacement step
#   cell_vec<-sample(cell_vec,replace = T)
#   
#   cell_vec<-unlist(lapply(cell_vec, function(vaf) {
#     
#     if(vaf==0) {stop(return(c(0,0)))}
#     
#     #existing mtDNA molecules are chosen at random to replicate during doubling of the mtDNA copy number
#     new_mut_mitos<-rbinom(1,size=mito_CN,prob=vaf)
#     total_mut_mito<-new_mut_mitos+round(vaf*mito_CN)
#     
#     #Mutant mitoDNA molecules are randomly partitioned between daughter cells
#     #Sampling is based on the idea of hypergeometric distribution as mutant mtDNA is sampled WITH REPLACEMENT
#     daughter1_mut_mito<-rhyper(nn=1,m=total_mut_mito,n=(2*mito_CN)-total_mut_mito,k=mito_CN)
#     daughter2_mut_mito=total_mut_mito-daughter1_mut_mito
#     c(daughter1_mut_mito/mito_CN,daughter2_mut_mito/mito_CN)
#   }))
# }
# 
# 
# 
# 
# model4_res<-lapply(vafs_for_testing,function(starting_vaf) {
#   print(starting_vaf)
#   final_colony_vaf<-sapply(1:nsim,function(i) {
#     if(i%%10==0) {cat(i,sep="\n")}
#     cell_vec<-starting_vaf
#     for(gen in 1:ngen) {
#       
#       #Add in a sampling with replacement step
#       cell_vec<-sample(cell_vec,replace = T)
#       
#       cell_vec<-unlist(lapply(cell_vec, function(vaf) {
#         
#         if(vaf==0) {stop(return(c(0,0)))}
#         
#         #existing mtDNA molecules are chosen at random to replicate during doubling of the mtDNA copy number
#         new_mut_mitos<-rbinom(1,size=mito_CN,prob=vaf)
#         total_mut_mito<-new_mut_mitos+round(vaf*mito_CN)
#         
#         #Mutant mitoDNA molecules are randomly partitioned between daughter cells
#         #Sampling is based on the idea of hypergeometric distribution as mutant mtDNA is sampled WITH REPLACEMENT
#         daughter1_mut_mito<-rhyper(nn=1,m=total_mut_mito,n=(2*mito_CN)-total_mut_mito,k=mito_CN)
#         daughter2_mut_mito=total_mut_mito-daughter1_mut_mito
#         c(daughter1_mut_mito/mito_CN,daughter2_mut_mito/mito_CN)
#       }))
#     }
#     return(mean(cell_vec))
#   })
#   
#   as.data.frame(t(quantile(final_colony_vaf,c(0.025,0.5,0.975))))%>%
#     dplyr::mutate(initial_cell_vaf=starting_vaf,.before=1)
# })%>%dplyr::bind_rows()
# 
# ###Combined visualization
# 
# dplyr::bind_rows(model1_res%>%mutate(model="model_1",.before=1),
#                  model2_res%>%mutate(model="model_2",.before=1),
#                  model4_res%>%mutate(model="model_4",.before=3)
# )%>%ggplot(aes(x=initial_cell_vaf,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+
#   geom_point()+
#   geom_errorbar(width=0.01)+
#   theme_classic()+
#   scale_x_continuous(limits=c(0,1))+
#   scale_y_continuous(limits=c(0,1))+
#   geom_abline(intercept=0,slope=1)+
#   facet_grid(~model)+
#   theme()
# 

