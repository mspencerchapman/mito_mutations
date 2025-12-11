library(optparse)
library(rsimpop)
library(ggplot2)
library(dplyr)

option_list = list(
  make_option(c("-n", "--n_of_job_run"), action="store", default=1, type='numeric', help="index of simulation so that results are distinguishable"),
  make_option(c("-s", "--sims_per_job"), action="store", default=1, type='numeric', help="number of simulations per script"),
  make_option(c("-t", "--test_run"), action="store_true", default=FALSE, type='logical', help="option to run test parameters rather than pick from the priors"),
  make_option(c("-p", "--pop_size"), action="store", default=NULL, type='numeric', help="option to set the log10 cell population size within the job"),
  make_option(c("-o", "--output_dir"), action="store", default='.', type='character', help="output directory")
)

#Parse the inputted options
opt = parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))
print(opt)

runs_per_job=opt$s
log_total_cell_pop=opt$p
params_file=paste0(opt$o,"/parameters_",opt$n,".txt")
sumstats_file=paste0(opt$o,"/summarystats_",opt$n,".txt")

#Set up function for mitochondrial mutation and drift between CELL GENERATIONS
between_cell_generation_drift=function(mut_vec,
                                       mito_fitness_vec,
                                       ngen,
                                       mutation_rate,
                                       non_synonymous_prob=0.5, #What is the probability of having a 'neutral' vs 'non-neutral' event
                                       prop_of_ns_affecting_function=0.2,
                                       mito_fitness_func=fitnessGammaFn) {
  
  mito_CN<-length(mut_vec)
  for(i in 1:ngen) {
    new_muts=rbinom(n=1,size=mito_CN,prob=mutation_rate)
    if(new_muts>0) {
      
      #Decide if synonymous or non-synonymous
      ns_muts=rbinom(n=new_muts,size=1,prob=non_synonymous_prob)
      function_altering_muts=rbinom(n=length(ns_muts),size=1,prob=prop_of_ns_affecting_function)
      
      #Generate mutation ids for the new mutations and assign fitness coefficients
      mut_IDs=ids::random_id(n=new_muts,bytes=4)
      mut_fitness=1+sapply(1:new_muts,function(i) mito_fitness_func() )
      mut_IDs<-paste(ifelse(ns_muts==1,"n","s"),mut_IDs,sep = "_")
      names(mut_fitness)<-mut_IDs
      
      #However, for synonymous mutations, reset fitness to 1
      mut_fitness[ns_muts==0|function_altering_muts==0]<-1
      
      #Add these coefficients to the vector recording the fitness of each mutation
      mito_fitness_vec<-c(mito_fitness_vec,mut_fitness)
      
      #Decide which mitochondria will be mutated, and replace the id
      for(j in 1:new_muts) {
        mut_select=sample(1:mito_CN,size = 1)
        original_mut_id=mut_vec[mut_select]
        
        if(mut_vec[mut_select]=="WT"){
          #If original ID is simply "WT" (wild type), update it with the mutation id
          mut_vec[mut_select]<-mut_IDs[j]
        } else {
          #If molecule already has a mutation, create combined mutation id
          combined_mut_ID<-paste(mut_vec[mut_select],mut_IDs[j],sep="-")
          mut_vec[mut_select]<-combined_mut_ID
          
          mito_fitness_vec[mut_IDs[j]]<-mito_fitness_vec[mut_IDs[j]]*mito_fitness_vec[original_mut_id]
          names(mito_fitness_vec)[which(names(mito_fitness_vec)==mut_IDs[j])] <- combined_mut_ID
        }
      }
    }
    new_mut_vec=sample(mut_vec,prob = mito_fitness_vec[mut_vec],replace = T)
    mut_vec<-new_mut_vec
  }
  return(list(mut_vec=mut_vec,mito_fitness_vec=mito_fitness_vec))
}

perform_cell_generation=function(total_cell_pop_list,
                                 mito_fitness_vec,
                                 n_mito_gen,
                                 mito_mutation_rate,
                                 vaf_to_fitness_func,
                                 mito_fitness_func,
                                 non_synonymous_prob,
                                 prop_of_ns_affecting_function,
                                 perform_resample=T) {
  
  #Perform the mitochondrial drift
  post_drift_res=lapply(total_cell_pop_list,function(mut_vec) {
    res=between_cell_generation_drift(mut_vec = mut_vec,
                                      mito_fitness_vec = mito_fitness_vec,
                                      ngen=n_mito_gen,
                                      mutation_rate=mito_mutation_rate,
                                      mito_fitness_func=mito_fitness_func,
                                      prop_of_ns_affecting_function,
                                      non_synonymous_prob=non_synonymous_prob)
    
    #Remove mutations from the fitness vector that are not in the population
    res$mito_fitness_vec<-res$mito_fitness_vec[unique(res$mut_vec)]
    
    return(res)
  })
  
  #Extract the new population list
  new_total_cell_pop_list=lapply(post_drift_res,function(list) list$mut_vec)
  
  #Extract the new mitochondrial fitness vector
  temp=unlist(lapply(post_drift_res,function(list) list$mito_fitness_vec))
  new_mito_fitness_vec<-temp[unique(names(temp))]
  
  #Now, from the cell population, use the mitochondrial mutations to work out the CELL FITNESS
  #This will depend on the proportion of mitochondria non-synonymous mutations
  
  if(perform_resample) {
    #Work out the proportion of non-synonymous mutations
    VAF_ns=sapply(new_total_cell_pop_list,function(mut_vec) {
      return(sum(table(mut_vec)[new_mito_fitness_vec[names(table(mut_vec))]>1])/mito_CN)
    })
    
    #Convert the VAFs to a fitness
    cell_fitness=sapply(VAF_ns,vaf_to_fitness_func)
    new_idx=sample(1:total_cell_pop,replace = T,prob = cell_fitness)
    resampled_total_cell_pop_list=new_total_cell_pop_list[new_idx]
    
    #Discard IDs that aren't in the resample cell population
    new_mito_fitness_vec=new_mito_fitness_vec[unique(unlist(resampled_total_cell_pop_list))]
    
    return(list(total_cell_pop_list=resampled_total_cell_pop_list,
                mito_fitness_vec=new_mito_fitness_vec))
  } else {
    return(list(total_cell_pop_list=new_total_cell_pop_list,
                mito_fitness_vec=new_mito_fitness_vec))
  }
  
}

mature_compartment=function(mito_CN,
                            total_cell_pop,
                            n_cell_gen,
                            n_mito_gen = 10,
                            mito_mutation_rate=1e-5,
                            non_synonymous_prob,
                            prop_of_ns_affecting_function,
                            vaf_to_fitness_func=dec_sigmoid,
                            mito_fitness_func,
                            perform_resample,
                            verbose=F) { #verbose mode will also record the mean sum of vaf measure for each generation in a 'ledger'
  
  #Set up initial population list from the simulation parameters
  total_cell_pop_list=list()
  for(k in 1:total_cell_pop) {total_cell_pop_list[[k]]<-rep("WT",mito_CN)}
  
  #Set up the mitochondrial fitness vector - this records the fitness of each mitochondrial "clone"
  mito_fitness_vec=1
  names(mito_fitness_vec)="WT"
  mean_sov=NULL
  
  if(verbose) {
    cat(paste("Cell compartment contains",length(total_cell_pop_list),"cells"),sep="\n")
    cat(paste("Performing",n_cell_gen,"cell generations, with",n_mito_gen,"mitochondrial generation between each cell generation."),sep="\n")
    cat(paste("After each cell generation, cells will",ifelse(perform_resample,"","not"),"be sampled with replacement."),sep="\n")
  }
  for(l in 1:n_cell_gen) {
    cat(paste("Starting cell generation",l),sep="\n")
    this_res=perform_cell_generation(total_cell_pop_list = total_cell_pop_list,
                                     mito_fitness_vec=mito_fitness_vec,
                                     n_mito_gen = n_mito_gen,
                                     mito_mutation_rate=mito_mutation_rate,
                                     vaf_to_fitness_func=vaf_to_fitness_func,
                                     mito_fitness_func=mito_fitness_func,
                                     non_synonymous_prob=non_synonymous_prob,
                                     prop_of_ns_affecting_function,
                                     perform_resample=perform_resample)
    total_cell_pop_list=this_res$total_cell_pop_list
    mito_fitness_vec=this_res$mito_fitness_vec
    
    if(verbose) {
      cat(paste("Total mtDNA mutations in compartment =",length(mito_fitness_vec)),sep="\n")
      fully_WT_cells<-sum(sapply(total_cell_pop_list,function(x) all(x=="WT")))
      cat(paste(fully_WT_cells,"cells have no mtDNA mutations"),sep = "\n")
      sov=sapply(total_cell_pop_list,function(mut_vec) {
        geno<-unlist(sapply(mut_vec,stringr::str_split,"-"))
        geno<-geno[geno!="WT"]
        return(length(geno)/length(mut_vec))
      })
      mean_sov=c(mean_sov,mean(sov))
      cat(paste("The mean 'sum of VAF' amongst cells is",round(mean(sov),digits = 4)),sep = "\n")
    }
  }
  if(verbose) {
    return(list(total_cell_pop_list=total_cell_pop_list,
                mito_fitness_vec=mito_fitness_vec,
                mean_sum_of_vaf_ledger=mean_sov))
  } else {
    return(list(total_cell_pop_list=total_cell_pop_list,
                mito_fitness_vec=mito_fitness_vec))
  }
  
}

##SET A RUNID
params_list=list()
sumstats_list=list()

for(x in 1:runs_per_job) {
  	
	cat(paste("Simulation number",x),sep="\n")
	run_id=ids::random_id(n = 1, bytes = 16, use_openssl = TRUE)
  
  if(opt$t) { #for test runs, choose quick parameters
    cell_fitness_param_Th=0.05 #Value at which
    cell_fitness_param_h=0.6 #Half point Valuei.e. at which the selective coeffecient is 0.5 
    cell_fitness_param_n=3 #Steepness of the drop
    
    #2. For the mitochondrial fitness distribution
    fitness_threshold=0.0001
    gamma_shape=0.5#0.5
    gamma_rate=200 #200
    
    #3. For the population dynamics
    mito_CN=500 #700
    total_mitochondrial_generations=1400#1400
    n_cell_gen=70 #70 #Number of cell cycles
    n_mito_gen = round(total_mitochondrial_generations/n_cell_gen) #Number of mitochondrial generations between cell cycles
    
    log_total_cell_pop=4
    total_cell_pop=round(10^log_total_cell_pop)
    log_mito_mutation_rate=-3.3
    mito_mutation_rate=10^log_mito_mutation_rate #This is the mutation rate measure as the number of new mutations per mitochondrial replication
    prop_of_ns_under_selection=0.25
    non_synonymous_prob=0.5
    
    #Create a list of the parameters to save
    params=list(run_id=run_id,
                cell_fitness_param_Th=cell_fitness_param_Th,
                cell_fitness_param_h=cell_fitness_param_h,
                cell_fitness_param_n=cell_fitness_param_n,
                fitness_threshold=fitness_threshold,
                gamma_shape=gamma_shape,
                gamma_rate=gamma_rate,
                mito_CN=mito_CN,
                total_mitochondrial_generations=total_mitochondrial_generations,
                n_cell_gen=n_cell_gen,
                n_mito_gen=n_mito_gen,
                log_total_cell_pop=log_total_cell_pop,
                log_mito_mutation_rate=log_mito_mutation_rate,
                prop_of_ns_under_selection=prop_of_ns_under_selection,
                non_synonymous_prob=non_synonymous_prob)
    
    params_df=as.data.frame(t(unlist(params)))
    
  } else {
    #Draw parameters from prior distributions
    #1. For the cell fitness distribution
    cell_fitness_param_Th=runif(1,min=0,max=0.6) #Value at which
    cell_fitness_param_h=runif(1,min=0.1,max=1) #Half point Valuei.e. at which the selective coeffecient is 0.5 
    cell_fitness_param_n=runif(1,min=0,max=10) #Steepness of the drop
    
    #2. For the mitochondrial fitness distribution
    fitness_threshold=0.0001
    gamma_shape=runif(1,0.1,1)#0.5
    gamma_rate=runif(1,min=10,1000) #200
    
    #3. For the population dynamics
    mito_CN=round(runif(1,500,1200)) #700
    total_mitochondrial_generations=round(runif(1,min=1000,2000)) #1400
    n_cell_gen=round(runif(1,min=30,max=200)) #70 #Number of cell cycles
    n_mito_gen = round(total_mitochondrial_generations/n_cell_gen) #Number of mitochondrial generations between cell cycles
    
    if(is.null(log_total_cell_pop)) {
      log_total_cell_pop=runif(1,min=3,max=4.9)#2e3
    }
    total_cell_pop=round(10^log_total_cell_pop)
    log_mito_mutation_rate=runif(1,min=-3.7,max=-3.2)
    mito_mutation_rate=10^log_mito_mutation_rate #This is the mutation rate measure as the number of new mutations per mitochondrial replication
    prop_of_ns_under_selection=runif(1,min=0,max=1)
    non_synonymous_prob=0.7
    
    #Create a list of the parameters to save
    params=list(run_id=run_id,
                cell_fitness_param_Th=cell_fitness_param_Th,
                cell_fitness_param_h=cell_fitness_param_h,
                cell_fitness_param_n=cell_fitness_param_n,
                fitness_threshold=fitness_threshold,
                gamma_shape=gamma_shape,
                gamma_rate=gamma_rate,
                mito_CN=mito_CN,
                total_mitochondrial_generations=total_mitochondrial_generations,
                n_cell_gen=n_cell_gen,
                n_mito_gen=n_mito_gen,
                log_total_cell_pop=log_total_cell_pop,
                log_mito_mutation_rate=log_mito_mutation_rate,
                prop_of_ns_under_selection=prop_of_ns_under_selection,
                non_synonymous_prob=non_synonymous_prob)
    
    params_df=as.data.frame(t(unlist(params)))
    
  }
  
  print(params_df)
  
  ##Now define a function to do a cell generation of drift
  #The 'Th' (Threshold) shifts the curve such that the effect on fitness only starts above that level
  dec_sigmoid=function(X,
                       Th=cell_fitness_param_Th,
                       m=1,
                       h=cell_fitness_param_h,
                       n=cell_fitness_param_n) {
    if(X<Th) {
      return(1)
    } else {
      m * h^n / (h^n + (X-Th)^n)
    }
  }
  #plot(sapply(seq(0,1,0.01),dec_sigmoid)~seq(0,1,0.01),ylim=c(0,1),xlim=c(0,1))
  
  #Define function for generating mitochondrial fitness
  genGammaFitness=function(fitness_threshold,shape,rate){
    function() rtrunc(n=1,a=fitness_threshold, b=Inf,"gamma",shape=shape,rate=rate)
  }
  fitnessGammaFn=genGammaFitness(fitness_threshold=fitness_threshold,shape = gamma_shape, rate=gamma_rate)
  
  #Visualize the non-synonymous mitochondrial fitness function
  #plot(density(sapply(1:1e4,function(i) fitnessGammaFn())),xlim=c(0,0.05))
  
  res=mature_compartment(mito_CN=mito_CN,
                         n_cell_gen=n_cell_gen,#n_cell_gen, #Number of cell cycles
                         total_cell_pop=total_cell_pop,
                         n_mito_gen = n_mito_gen, #Number of mitochondrial generations between cell cycles
                         mito_mutation_rate=mito_mutation_rate, #This is the mutation rate measure as the number of new mutations per mitochondrial replication
                         non_synonymous_prob=non_synonymous_prob,
                         prop_of_ns_affecting_function=prop_of_ns_under_selection,
                         vaf_to_fitness_func=dec_sigmoid, #links the proportion of mtDNA in a cell with a non-synonymous mutation with the overall cellular fitness
                         mito_fitness_func=fitnessGammaFn,
                         perform_resample=T) #The default is to resample the cells (based on cell fitness) after each generation. Use 'FALSE' to suppress this behvaiour
  
  ##############################################################
  #-----CALCULATE THE SUMMARY STATISTICS FOLLOWING SIMULATION--------
  ##############################################################
  
  cat("Calculating the summary statistics",sep = "\n")
  
  bins=c(0.01,0.05,0.1,0.2,0.5,1)
  bin_order=c("<1%","1-5%","5-10%","10-20%","20-50%",">50%")
  
  #Perform a subsampling step to reflect experimental protocol
  n_cells_sampled=300
  sub_res=res
  sub_res$total_cell_pop_list=res$total_cell_pop_list[sort(sample(1:length(res$total_cell_pop_list),size = n_cells_sampled,replace = F))]
  sub_res$mito_fitness_vec=res$mito_fitness_vec[unique(unlist(sub_res$total_cell_pop_list))]
  
  #Now extract the sumstats
  unique_muts<-unique(unlist(stringr::str_split(names(sub_res$mito_fitness_vec),pattern = "-")))[-1]
  mut_vaf_df=lapply(unique_muts, function(mut) {
    vafs=sapply(sub_res$total_cell_pop_list,function(mut_vec) {
      sum(grepl(pattern = mut,mut_vec))/length(mut_vec)
    })
    data.frame(mut_ID=mut,mut_type=substr(mut,1,1),ncell=sum(vafs>0),max_vaf=max(vafs))
  })%>%dplyr::bind_rows()
  
  # ggplot(mut_vaf_df,aes(x=max_vaf,fill=mut_type))+
  #   scale_x_log10()+
  #   facet_grid(~mut_type)+
  #   geom_histogram(alpha=1)
  
  dNdS_sumstats<-mut_vaf_df%>%
    mutate(bin=ifelse(max_vaf<0.01,"<1%",ifelse(max_vaf>=0.01&max_vaf<0.05,"1-5%",ifelse(max_vaf>=0.05&max_vaf<0.1,"5-10%",ifelse(max_vaf>=0.1&max_vaf<0.2,"10-20%",ifelse(max_vaf>=0.2&max_vaf<0.5,"20-50%",">50%"))))))%>%
    group_by(bin)%>%
    summarise(syn=sum(mut_type=="s"),non_syn=sum(mut_type=="n"),total=n())%>%
    mutate(total_per_cell=total/n_cells_sampled)%>%
    mutate(dNdS=(non_syn/syn)/(non_synonymous_prob/(1-non_synonymous_prob)),bin=factor(bin,levels=bin_order))%>%
    arrange(bin)
  
  # ggplot(dNdS_sumstats,aes(x=bin,y=dNdS))+
  #   geom_point()+
  #   geom_hline(yintercept=1,linetype=2)+
  #   theme_classic()  
  
  #Now get the mutation burden info
  sum_of_vaf=sapply(res$total_cell_pop_list,function(mut_vec) {
    mut_ids<-unlist(stringr::str_split(mut_vec,pattern = "-"))
    mut_ids<-mut_ids[mut_ids!="WT"]
    return(length(mut_ids)/length(mut_vec))
  })
  
  mut_burden_sumstats=quantile(sum_of_vaf,seq(0,1,0.1))
  mut_burden_sumstats_df<-t(as.data.frame(mut_burden_sumstats))
  colnames(mut_burden_sumstats_df)<-paste("mut_burden_quant",names(mut_burden_sumstats),sep = "_")
  dNdS_sumstats_df<-dNdS_sumstats%>%dplyr::select(-total,-syn,-non_syn)%>%tidyr::pivot_wider(names_from = "bin",values_from = c("total_per_cell","dNdS"))
  all_sumstats<-cbind(data.frame(run_id=run_id),mut_burden_sumstats_df,dNdS_sumstats_df)%>%tibble::remove_rownames()
  
  sumstats_list[[x]]<-all_sumstats
  params_list[[x]]<-params_df
}

cat("Completed simulations",sep="\n")
#Condense the lists into dataframes
params_comb<-dplyr::bind_rows(params_list)
sumstats_comb<-dplyr::bind_rows(sumstats_list)

#Write the dataframes into tables
cat("Writing out the parameters and summary stats",sep="\n")
write.table(params_comb,file = params_file,quote = F,row.names = F,col.names = T,sep = "\t")
write.table(sumstats_comb,file = sumstats_file,quote = F,row.names = F,col.names = T,sep = "\t")

cat("Completed.",sep="\n")
