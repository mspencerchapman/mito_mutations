##### Drift through bottleneck

mito_CN_ss=675
mito_CN_bn=10
ngen_pre_bn=50
ngen_post_bn=5
starting_vaf=0.1

print(starting_vaf)
post_ss_vaf_list<-lapply(1:nsim,function(i) {
  if(i%%10==0) {cat(i,sep="\n")}
  cell_vec<-starting_vaf
  for(gen in 1:ngen_pre_bn) {
    
    cell_vec<-unlist(lapply(cell_vec, function(vaf) {
      
      if(vaf==0) {stop(return(0))}
      
      #existing mtDNA molecules are chosen at random to replicate during doubling of the mtDNA copy number
      new_mut_mitos<-rbinom(1,size=mito_CN_ss,prob=vaf)
      total_mut_mito<-new_mut_mitos+round(vaf*mito_CN_ss)
      
      #Mutant mitoDNA molecules are randomly partitioned between daughter cells
      #Sampling is based on the idea of hypergeometric distribution as mutant mtDNA is sampled WITH REPLACEMENT
      daughter1_mut_mito<-rhyper(nn=1,m=total_mut_mito,n=(2*mito_CN_ss)-total_mut_mito,k=mito_CN_ss)
      daughter2_mut_mito=total_mut_mito-daughter1_mut_mito
      #c(daughter1_mut_mito/mito_CN_ss,daughter2_mut_mito/mito_CN_ss)
      c(daughter1_mut_mito/mito_CN_ss)
    }))
  }
  return(cell_vec)
})

post_bn_cell_vec<-lapply(post_ss_vaf_list, function(vaf) {
  if(vaf==0) {stop(return(0))}
  
  #existing mtDNA molecules are chosen at random to replicate during doubling of the mtDNA copy number
  new_mut_mitos<-rbinom(1,size=mito_CN_bn,prob=vaf)
  c(new_mut_mitos/mito_CN_bn)
})

final_cell_vec<-lapply(post_bn_cell_vec, function(vaf) {
  
  if(vaf==0) {stop(return(c(0,0)))}
  
  #existing mtDNA molecules are chosen at random to replicate during doubling of the mtDNA copy number
  new_mut_mitos<-rbinom(1,size=mito_CN_ss,prob=vaf)
  total_mut_mito<-new_mut_mitos+round(vaf*mito_CN_ss)
  
  #Mutant mitoDNA molecules are randomly partitioned between daughter cells
  #Sampling is based on the idea of hypergeometric distribution as mutant mtDNA is sampled WITH REPLACEMENT
  daughter1_mut_mito<-rhyper(nn=1,m=total_mut_mito,n=(2*mito_CN_ss)-total_mut_mito,k=mito_CN_ss)
  daughter2_mut_mito=total_mut_mito-daughter1_mut_mito
  #c(daughter1_mut_mito/mito_CN_ss,daughter2_mut_mito/mito_CN_ss)
  c(daughter1_mut_mito/mito_CN_ss)
})

mean(unlist(final_cell_vec))



