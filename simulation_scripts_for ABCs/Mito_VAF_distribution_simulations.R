library(dplyr)
library(tidyr)
library(stringr)

#Define function for simulating allelic drift according to a F-W model
#This is parameterised by the population size and number of generations
fisher_wright_drift=function(starting_vaf,
                             population_size,
                             generation_time,
                             total_time){
  curr_vaf<-starting_vaf
  total_generations=round(total_time/generation_time)
  for(i in 1:total_generations) {
    curr_vaf<-rbinom(1,size = population_size,prob = curr_vaf)/population_size #This is equivalent to sampling with replacement
    if(curr_vaf==0|curr_vaf==1){break} #if the VAF becomes 0 or 1, no point in continuing - mutation is lost or fixed
  }
  return(curr_vaf)
}

#Define the VAF 'bins' that will be used for comparing VAF distributions
new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")
VAF_groups=data.frame(
  labels=new_VAF_groups,
  lower_limit=c(0,2^(-10:-1)),
  upper_limit=2^(-10:0)
)

#START THE ACTUAL SIMULATIONS
nsim=10000
all_sumstats=lapply(1:nsim,function(k) {
  if(k%%10==0) {print(k)}
  
  #Define the fixed parameters - these remain constant in all simulations
  mito_copy_number=600 #This is approximately the average copy number from the adult individuals
  number_of_cells_in_simulation=1000 #Need enough 'cells' to have adequate mutation numbers to define the distribution
  
  #Select the varying parameters from their prior distributions
  muts_per_mitochondria_per_generation=runif(1,1e-4,1e-3)
  ngen=round(runif(1,1,3000))
  
  #Use these parameters to introduce mutations and allow them to drift up until the time of sampling
  #NB. Mutations are introduced each generation & allowed to drift until the end of simulation
  # Therefore, if the total number of generations is 1000, some mutations will be introduced in the first
  #generation and will undergo 1000 generations of drift.  Some mutations will be introduced later (e.g. during
  #the 600th generation) and will have fewer remaining generations of drift (e.g. 400 remaining
  #generations).  The final VAF distribution is the aggregate of the final VAFs of all the mutations
  #introduced over the lifetime.  Mutations are introduced at the same rate in each
  #generation, but with some poisson variation in absolute number.
  
  gen_list=list()
  #Iterate through each generation defining the
  for(j in 1:ngen){
    mutation_introductions_this_generation=sum(rpois(number_of_cells_in_simulation*mito_copy_number,lambda=muts_per_mitochondria_per_generation))
    gen_mut_vafs<-vector(length=mutation_introductions_this_generation)
    for(i in 1:mutation_introductions_this_generation) {
      this_mut_final_vaf<- fisher_wright_drift(starting_vaf=1/mito_copy_number,
                                               population_size = mito_copy_number,
                                               generation_time=1,
                                               total_time=ngen-j+1)
      gen_mut_vafs[i]<-this_mut_final_vaf
    }
    gen_list[[j]]<-gen_mut_vafs
  }
  
  final_vafs=Map(gen_muts=gen_list,gen=1:ngen,function(gen_muts,gen) {
    if(length(gen_muts[gen_muts>0])>0) {
      return(data.frame(total_generations=ngen,mut_rate=muts_per_mitochondria_per_generation,gen=gen,gen_mut_vafs=gen_muts[gen_muts>0]))
    } else {
      return(NULL)
    }
    
  })%>%dplyr::bind_rows()
  
  #The observed VAF will not always perfectly reflect the true VAF (by chance) and depends on the sequencing coverage
  #Here, 7.5 is the 'haploid sequencing coverage' in this simulated experiment. This approximates the ~15X diploid coverage of most of the experiments
  haploid_coverage=7.5
  final_vafs$observed_vaf=sapply(final_vafs$gen_mut_vafs,function(vaf) rbinom(n=1,size=(haploid_coverage*mito_copy_number),prob=vaf)/(haploid_coverage*mito_copy_number))
  
  #Assign each observed VAF into its 'VAF group' bin
  final_vafs$VAF_group=sapply(1:nrow(final_vafs),function(i) {
    if(final_vafs$observed_vaf[i]==0) {
      return("Absent")
    } else {
      return(VAF_groups$labels[VAF_groups$lower_limit<final_vafs$observed_vaf[i] & VAF_groups$upper_limit>=final_vafs$observed_vaf[i]])
    }
  })
  
  #Summarise these into numbers of mutations in each bin per 'simulated cell' in the simulation
  sumstats<-final_vafs%>%
    group_by(total_generations,VAF_group)%>%
    summarise(n=n())%>%
    filter(VAF_group!="Absent")%>%
    mutate(VAF_group=factor(VAF_group,levels=VAF_groups$labels))%>%
    mutate(muts_per_cell=n/number_of_cells_in_simulation)%>%
    dplyr::select(-n)%>%
    tidyr::pivot_wider(names_from="VAF_group",values_from="muts_per_cell")%>%
    mutate(muts_per_mitochondria_per_generation=muts_per_mitochondria_per_generation,.before=total_generations)
  
  return(sumstats)
})%>%dplyr::bind_rows()

#Save the output for doing the ABC
saveRDS(all_sumstats,file = "VAF_distribution_ABC_simulation_sumstats_",ids::random_id(),".Rds")