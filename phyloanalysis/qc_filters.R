
do_hard_minvaf_filter=function(inf){
  clones=get_clones(inf)
  minvaf=inf$params$hard_minvaf_limit
  for(x in c("snv","indel")){
    # Sums of MTR and DEP at mutant genotype sites
    MTR=rowSums(inf[[x]]$mtr[,clones]*(inf[[x]]$geno>0),na.rm=TRUE)
    DEP=rowSums(inf[[x]]$dep[,clones]*(inf[[x]]$geno>0),na.rm=TRUE)
    VAFLIMIT=rowSums(inf[[x]]$dep[,clones]*(inf[[x]]$geno>0)*inf[[x]]$minmutvaf[,clones],na.rm=TRUE)/DEP
    ##  This accounts for samples specific LOH
    VAFLIMIT=(minvaf/0.5)*VAFLIMIT
    inf[[x]]$details$meanvaf_g=MTR/DEP
    ##  meanvaf_g will be NA if there are no variants that have a mutant genotype (e.g. will also fail count)
    inf[[x]]$details=inf[[x]]$details %>% mutate(filter_hard_minvaf=ifelse(is.na(meanvaf_g) | meanvaf_g<VAFLIMIT,1,0))
  }
  inf
}


do_simple_vaf_filter=function(inf){
  clones=get_clones(inf)
  
  for(x in c("snv","indel")){
    # Sums of MTR and DEP at mutant genotype sites
    MTR=rowSums(inf[[x]]$mtr[,clones]*(inf[[x]]$geno>0),na.rm=TRUE)
    DEP=rowSums(inf[[x]]$dep[,clones]*(inf[[x]]$geno>0),na.rm=TRUE)
    # Calculate the sample aggregate VAF limit taking into account copy number events via minmutvaf (also includes MINACF)
    VAFLIMIT=rowSums(inf[[x]]$dep[,clones]*(inf[[x]]$geno>0)*inf[[x]]$minmutvaf[,clones],na.rm=TRUE)/DEP
    # This p-value assesses whether VAF is lower than is consistent with the VAFLIMIT
    p=pbinom(MTR,DEP,p=VAFLIMIT)
    ##Similar stats at WT sites
    MTR2=rowSums(inf[[x]]$mtr[,clones]*(inf[[x]]$geno==0),na.rm=TRUE)
    DEP2=rowSums(inf[[x]]$dep[,clones]*(inf[[x]]$geno==0),na.rm=TRUE)
    ## Reject if "error rate" is significantly>0.01 # for mouse try 0.05
    p2=pbinom(MTR2-1,DEP2,p=0.05,lower.tail=FALSE) # equivalent to binom.test(MTR2,DEP2,"greater")
    ## Similar to above except rather than using mutant genotype require at least 2 mutant reads for assessment
    MTR3=rowSums(inf[[x]]$mtr[,clones]*(inf[[x]]$mtr[,clones] >= inf$params$simple_vaf_minm),na.rm=TRUE)
    DEP3=rowSums(inf[[x]]$dep[,clones]*(inf[[x]]$mtr[,clones] >= inf$params$simple_vaf_minm),na.rm=TRUE)
    VAFLIMIT3=rowSums(inf[[x]]$dep[,clones]*(inf[[x]]$mtr[,clones]>1)*inf[[x]]$minmutvaf[,clones],na.rm=TRUE)/DEP3
    p3=pbinom(MTR3,DEP3,p=VAFLIMIT3) # equivalent to binom.test(MTR,DEP,"less")
    inf[[x]]$details$meanvaf_g=MTR/DEP
    
    # VAF too low based on VAF at genotype=mutant sites
    inf[[x]]$details$filter_vaf_too_low_s=ifelse(p<0.001,1,0)
    inf[[x]]$details$filter_vaf_too_low_m2=ifelse(p3<0.001,1,0)
    inf[[x]]$details$meanvaf_m2=MTR3/DEP3
    
    inf[[x]]$details$vaflimit_s=VAFLIMIT
    inf[[x]]$details$vaflimit_m2=VAFLIMIT3
    inf[[x]]$details$mtr_s=MTR
    inf[[x]]$details$dep_s=DEP
    inf[[x]]$details$mtr_m2=MTR3
    inf[[x]]$details$dep_m2=DEP3
    
    
    inf[[x]]$details$vaflimit_m2=VAFLIMIT3
    
    ##  VAF to high at zero genotyped sites.
    inf[[x]]$details$filter_vaf_zg_too_noisy=ifelse(p2<0.001,1,0)
    
    ##Record stats/p-values for posterity
    ##Aggregate VAF at wild type sites
    inf[[x]]$details$error_at_wt_sites=MTR2/DEP2
    inf[[x]]$details$p_error_at_wt_sites_lt_0_01=p2
    inf[[x]]$details$p_vaf_at_mt_sites_gt_0_45=p
    inf[[x]]$details$p_vaf_at_mtr_gt_1_sites_gt_0_45=p3
  }
  inf
}

do_beta_binomial_filter=function(inf){
  inf=apply_rho(inf,"snv")
  inf$snv$details$filter_beta_binomial=ifelse(inf$snv$details$rho>0.1,0,1)
  inf=apply_rho(inf,"indel")
  inf$indel$details$filter_beta_binomial=ifelse(inf$indel$details$rho>0.15,0,1)
  inf
}

##Sets the filter_near_indel flag. This removes SNVs that are  within 10 bp of neighbouring indels.
# Following few functions adapted from https://github.com/HLee-Six/HSC_phylodynamics/blob/master/filtering_subs/filter_mpileup_subs_pre_shearwater.R 
do_filter_near_indels=function(inf){
  if(dim(inf$indel$details)[1]==0){
    inf$snv$details$filter_near_indel=0
    return(inf)
  }
  vv=inf$snv$details
  inds=inf$indel$details[,1:4]
  inds$type <- NA
  ## TODO review the equality case below
  inds$type[(nchar(inds$Ref))>=(nchar(inds$Alt))] = "del"
  inds$type[(nchar(inds$Ref))<(nchar(inds$Alt))] = "ins"
  inds$Start <- NA
  inds$End <- NA
  inds$Start[inds$type=="ins"] <- inds$Pos[inds$type=="ins"]-10
  inds$End[inds$type=="ins"] <- inds$Pos[inds$type=="ins"]+10
  inds$Start[inds$type=="del"] <- inds$Pos[inds$type=="del"]-10
  inds$End[inds$type=="del"] <- inds$Pos[inds$type=="del"] + nchar(inds$Ref[inds$type=="del"]) + 10
  grinds <- GRanges(inds$Chrom, IRanges(start = inds$Start, end=inds$End))
  grsubs <- GRanges(vv$Chrom, IRanges(start=as.numeric(vv$Pos), width=1))
  overlaps <- findOverlaps(grsubs, grinds)
  nearindels <- unique(queryHits(overlaps))
  inf$snv$details$filter_near_indel=0
  inf$snv$details$filter_near_indel[nearindels]=1
  inf
}


get_too_close_idx=function(df,minbp=10){
  tt=df
  cc <- tt[with(tt, order(Chrom, Pos)),]
  uu <- data.frame()
  for (i in unique(cc$Chrom)) {
    thischr <- cc[cc$Chrom==i,]
    thischr$toprev <- NA
    # get the distance to the previous mutation
    if(dim(thischr)[1]>1){
      for (n in 2:nrow(thischr)) {
        #cat(i,"\t",n,"\n")
        dd <- thischr[n, "Pos"] - thischr[n-1, "Pos"]
        thischr$toprev[n] <- dd
      }
      thischr$tonext <- c(thischr$toprev[2:nrow(thischr)],minbp+1)
    }else{
      thischr$tonext = minbp+1
    }
    thischr$toprev[1] <- minbp+1 # no need to lose the first line.
    # put in the distance to the next one too
    #thischr$tonext <- c(thischr$toprev[2:nrow(thischr)],minbp+1)
    uu <- as.data.frame(rbind(uu, thischr))
  }
  idx.too.close=which(!(as.numeric(uu$toprev)>minbp & as.numeric(uu$tonext) > minbp))
}

##Sets "filter_too_close" filter
# Flags SNVs that are within 10 bp of another SNV
# Flags Indels that are within 10 bp of another Indel
do_filter_too_close=function(inf,minbp=10){
  idx.too.close.snv=get_too_close_idx(inf$snv$details,minbp=minbp)
  inf$snv$details$filter_too_close=0
  inf$snv$details$filter_too_close[idx.too.close.snv]=1
  idx.too.close.indel=get_too_close_idx(inf$indel$details,minbp=minbp)
  if(dim(inf$indel$details)[1]==0){
    inf$indel$details$filter_too_close=c()
  }else{
    inf$indel$details$filter_too_close=0
    inf$indel$details$filter_too_close[idx.too.close.indel]=1
  }
  inf
}
## Creat copy of functions with the correct name (keep the old for backward compatibility)
do_near_indels_filter=do_filter_near_indels
do_too_close_filter=do_filter_too_close


## Sets "filter_max_miss" 
## Flags sites 1/3rd of the samples are allowed to have depth<=6 
## and "filter_max_gmiss"
## Flags sites 1/3rd of the samples have missing genotypes (see add_counts)
do_missing_filter=function(inf,maxmiss=floor(length(inf$meta$clones)/3)){
  #max(5,floor(length(inf$meta$clones)/5)
  cat("applying max missing filter #NMISS<=%s",maxmiss,"\n")
  mtrs=inf$meta$mtrs
  deps=inf$meta$deps
  
  clones=get_clones(inf)#inf$meta$clones
  
  inf$snv$details$nmiss=apply(inf$snv$vaf_d[,clones],1,function(row) length(which(is.na(row))))
  inf$indel$details$nmiss=apply(inf$indel$vaf_d[,clones],1,function(row) length(which(is.na(row))))
  inf$snv$details$filter_max_miss=ifelse(inf$snv$details$nmiss>maxmiss,1,0)
  inf$indel$details$filter_max_miss=ifelse(inf$indel$details$nmiss>maxmiss,1,0)
  inf$snv$details$filter_max_gmiss=ifelse(inf$snv$details$g_nmiss>maxmiss,1,0)
  inf$indel$details$filter_max_gmiss=ifelse(inf$indel$details$g_nmiss>maxmiss,1,0)
  
  
  
  inf
}

## and "filter_max_gmiss"
## Flags sites 1/3rd of the samples have missing genotypes (see add_counts)
do_gmiss_filter=function(inf,maxmiss=floor(length(inf$meta$clones)/3)){
  #max(5,floor(length(inf$meta$clones)/5)
  cat("applying max missing filter #NMISS<=%s",maxmiss,"\n")
  mtrs=inf$meta$mtrs
  deps=inf$meta$deps
  clones=get_clones(inf)#inf$meta$clones
  inf$snv$details$filter_max_gmiss=ifelse(inf$snv$details$g_nmiss>maxmiss,1,0)
  inf$indel$details$filter_max_gmiss=ifelse(inf$indel$details$g_nmiss>maxmiss,1,0)
  inf
}


## Set "filter_count"
## Flags any site where genotype status indicates all colonies are mutant or no colonies are mutant.
do_count_filter=function(inf){
  clones=get_clones(inf)
  n_clone=length(clones)
  for(x in c("snv","indel")){
    idx=which((inf[[x]]$details$g_counts %in% 1:(n_clone-1)) &  ##At least one non-zero genotype##non-zero mutant count in some but not all clones
                inf[[x]]$details$ref_counts>0   ##At least one of the clones has no mutant reads..
    )
    if(!is.na(get_normal(inf))){
      ## Sometimes there is a clonal trunk - here we use the normal sample to identify likely somatic variants.
      idx=which(inf[[x]]$details$g_counts>0)
    }
    if(dim(inf[[x]]$details)[1]==0){
      inf[[x]]$details$filter_count=c() ##Handle empy indel file
    }else{
      inf[[x]]$details$filter_count=1
      inf[[x]]$details$filter_count[idx]=0
    }
  }
  inf
}


##Calculates GLOD 
get_glod=function(mut,ref,e,f=0.05){
  lmNonGermline=ref*log10((f*e/3+(1-f)*(1-e)))+mut*log10(f*(1-e)+(1-f)*e/3)
  lm05=(mut+ref)*log10(0.5*(e/3+(1-e)))
  LOD=lmNonGermline-lm05
  LOD
}

##Sets filter_bglod 
# Flags sites that are likely germline variants.
# In the absence of a normal sample all sites where aggregate mutant read count is not significantly less that MINACF*0.5 
# are assumed to be germline. Additionally 
do_bglod_filter=function(inf){
  if(is.na(inf$meta$normal)){
    warning("bglod filter unsupported for missing normal sample: using basic_germline instead")
    clones=get_clones(inf)
    for(x in c("snv","indel")){
      DEP=rowSums(inf[[x]]$dep[,clones])
      MTR=rowSums(inf[[x]]$mtr[,clones])
      ##  minmutvaf2 reflects the expected minimum vaf of germline variants.
      VAFLIMIT=rowSums(inf[[x]]$dep[,clones]*inf[[x]]$minmutvaf2[,clones],na.rm=TRUE)/DEP
      # one sided binimial test
      p.germline=pbinom(MTR,DEP,p=VAFLIMIT)
      inf[[x]]$details$filter_bglod=ifelse(p.germline<1e-3,0,1)
      inf[[x]]$details$allcolony_meanvaf=MTR/DEP
    }
    return(inf)
  }
  mtr=rowSums(inf$snv$mtr[,inf$meta$clones_short])
  dep=rowSums(inf$snv$dep[,inf$meta$clones_short])
  totvaf=mtr/dep
  #p=p.adjust(pbinom(mtr,dep,p=0.45,lower.tail = TRUE))
  mtri=rowSums(inf$indel$mtr[,inf$meta$clones_short])
  depi=rowSums(inf$indel$dep[,inf$meta$clones_short])
  totvafi=mtri/depi
  normal=get_normal(inf)
  normal_acf=inf$meta$NORMAL_ACF
  #Note that BGLOD is already calculated in the input file as part of the pileups process
  #The following makes the BGLOD filter less aggressive by effectively reinstating loci with sufficiently few mutant reads in normal.
  #The rationale is that these variants will end up on the trunk or not fitting the tree well and can be subsequently filtered.
  adjminmutvaf=apply(inf$snv$minmutvaf2,1,min,na.rm=TRUE)
  adjmaxmutvaf=apply(inf$snv$minmutvaf2,1,max,na.rm=TRUE)
  germlinevaf=0.95*apply(inf$snv$germlinevaf,1,max,na.rm=TRUE)## The 0.95 accounts for ref bias and also accomodates homozygous variants.
  mtrg=inf$snv$mtr[,normal]
  depg=inf$snv$dep[,normal]
  p=pbinom(mtrg,depg,germlinevaf,lower.tail = TRUE)
  # Calculate GLOD value under specified contamination and  maximum vaf of variant.
  inf$snv$details$GLOD=get_glod(inf$snv$mtr[,normal],inf$snv$dep[,normal]-inf$snv$mtr[,normal],0.001,normal_acf*adjmaxmutvaf)
  bglod=with(inf$snv$details,ifelse(ifelse(GSNP==0,GLOD>1.35,GLOD>5),0,1))
  # If we reject the hypothesis that the variant is heterozygous (VAF=0.5) (p<0.001) then the variant passes the filter.
  # If we are unable to reject the hypothesis e.g. because the read count is quite high then we use the result of the GLOD filter which asseses
  # the probability of the read counts resulting from 10% contamination rather than being a SNP.
  inf$snv$details$p_het_germline=p
  inf$snv$details$BGLOD=bglod
  inf$snv$details$BGLOD_ORIG=bglod
  # If the sample fails glod then we should reinstate if it has aggregate VAF<30%..
  idx=which(totvaf<0.6*adjminmutvaf)
  # This might let a few SNPs back in - but they will be subject to further significance based clonal VAF testing. 
  if(length(idx)>0){
    inf$snv$details$BGLOD[idx]=0
  }
  inf$snv$details$filter_bglod=inf$snv$details$BGLOD #rescued_bglod
  ##Indels
  adjminmutvaf=apply(inf$indel$minmutvaf2,1,min,na.rm=TRUE)
  adjmaxmutvaf=apply(inf$indel$minmutvaf2,1,max,na.rm=TRUE)
  ## Use higher error rate. Assume same relative prior for indels as snps. 
  inf$indel$details$GLOD=get_glod(inf$indel$mtr[,normal],inf$indel$dep[,normal]-inf$indel$mtr[,normal],0.01,normal_acf*adjmaxmutvaf)
  bglod=with(inf$indel$details,ifelse(ifelse(is.na(GSNP) | GSNP==0,GLOD>1.35,GLOD>5),0,1))
  ##
  inf$indel$details$BGLOD_ORIG=bglod
  inf$indel$details$BGLOD=bglod
  idx=which(totvafi<0.6*adjminmutvaf)
  # This might let a few SNPs back in - but they will be subject to further significance based clonal VAF testing. 
  if(length(idx)>0){
    inf$indel$details$BGLOD[idx]=0
  }
  # Consider replacing with BGLOD  - althought the following reflects the idea that we can be more lax with indels since they don't inform tree topology.
  inf$indel$details$filter_bglod=ifelse(is.na(inf$indel$vaf_d[,normal]) | (inf$indel$vaf_d[,normal]<0.2),0,inf$indel$details$BGLOD)
  inf
}

do_19_to_38_change_arm_filter=function(inf){
  if(!dir.exists("liftover_resources")){
    cat("Not generating informational xfilter_change_arm field as liftover not available")
    return(inf)
  }
  for(x in c("snv","indel")){
    details=get_cytoband_hg19_38(inf[[x]]$details)
    inf[[x]]$details$xfilter_change_arm=ifelse(substr(details$cyto_hg19,1,1)!=substr(details$cyto_hg38,1,1),1,0)
    inf[[x]]$details=cbind(inf[[x]]$details,details[,c("cyto_hg19","cyto_hg38")])
  }
  inf
}
##  Identifies sites that do have an aggregate VAF that is consistent with being germline
do_basic_germline_filter=function(inf,x){
  clones=get_clones(inf)
  for(x in c("snv","indel")){
    DEP=rowSums(inf[[x]]$dep[,clones])
    MTR=rowSums(inf[[x]]$mtr[,clones])
    ##  Note minmutvaf is a ploidy specific matrix that has previously been setup.
    VAFLIMIT=rowSums(inf[[x]]$dep[,clones]*inf[[x]]$minmutvaf[,clones],na.rm=TRUE)/DEP
    # one sided binimial test
    p.germline=pbinom(MTR,DEP,p=VAFLIMIT)
    inf[[x]]$details$filter_basic_germline=ifelse(p.germline<1e-3,0,1)
    inf[[x]]$details$allcolony_meanvaf=MTR/DEP
  }
  inf
}

