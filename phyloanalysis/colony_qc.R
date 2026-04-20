#
library(GenomicRanges)
#library(ape)
#library(phangorn)
library(yaml)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(parallel)
source("project.R")
source("qc_filters.R")
# The code here supports the following QC process

if(system("uname",intern = TRUE)=="Darwin"){
  MC_CORES=4
  LUSTRE="/Users/nw14/volumes/nw14_network/LUSTRE"
  NFS="/Users/nw14/volumes/nw14_network/NFS"
}else{
  MC_CORES=as.integer(Sys.getenv("LSB_MAX_NUM_PROCESSORS"))
  if(is.na(MC_CORES)){
    MC_CORES=4
  }
  LUSTRE="/lustre"
  NFS="/nfs"
}

#MINACF=0.9 -For clonal samples
#MINACF=0.7

##Adds meta data from the yaml file
add_meta_yaml=function(inf,prefix){
  #meta=list()
  #meta$SEX=inf$SEX

  cfg=inf$cfg
  LOH=lapply(inf$yamlf$LOH,function(loh){

    idx=which(sapply(inf$yamlf$SAMPLES,function(x){loh$LABEL %in% x$LOH}));
    samples=sapply(inf$yamlf$SAMPLES[idx],function(y) y$SHORT_LABEL,USE.NAMES=FALSE)
    loh$samples=samples
    loh
  }
  )
  CNA=lapply(inf$yamlf$CNA,function(cna){

    idx=which(sapply(inf$yamlf$SAMPLES,function(x){cna$LABEL %in% x$CNA}));
    samples=sapply(inf$yamlf$SAMPLES[idx],function(y) y$SHORT_LABEL,USE.NAMES=FALSE)
    cna$samples=samples
    cna
  }
  )
  label_l=cfg$LABEL[which(cfg$TYPE=="C")]
  label_s=cfg$SHORT_LABEL[which(cfg$TYPE=="C")]
  label_o=cfg$LABEL[which(cfg$TYPE!="C")]
  normal=cfg$LABEL[which(cfg$TYPE=="N")][1]
  for(field in c("LABEL","SHORT_LABEL")){
    if(length(unique(cfg[[field]]))!=length(cfg[[field]])){
      stop(sprintf("%s should contain distinct entries!",field))
    }
  }

  inf$meta=list(SEX=inf$yamlf$SEX,LOH=LOH,CNA=CNA,clones=label_l,clones_short=label_s,samples=cfg$LABEL,normal=normal,prefix=prefix)
  inf
}


load_data=function(config_dir,data_dir,PREFIX,exclude=NULL){
  yamlf_path=sprintf("%s/%s.yaml",config_dir,PREFIX)
  yamlf=yaml.load_file(yamlf_path)
  fnames=c("LABEL", "TYPE","PROJECT", "LCM","SHORT_LABEL") ##Removed VERSION as required field
  if(!all(fnames %in% names(yamlf$SAMPLES[[1]]))){
    stop(sprintf("%s does not contain all required fields:%s",yamlf_path,paste(fnames,collapse=",")))
  }
  cfg=do.call("rbind",
              lapply(yamlf$SAMPLES,
                     function(x) as.data.frame(
                       lapply(fnames,function(z) x[[z]]),col.names=fnames,stringsAsFactors=FALSE
                      )
              )
      )
  if(length(unique(cfg$SHORT_LABEL))<length(cfg$SHORT_LABEL)){
    stop("Duplicates in cfg$SHORT_LABEL")
  }

  if(!is.null(exclude)){
    if(any(is.na(match(exclude,cfg$SHORT_LABEL)))){
      stop("Can't find excluded all excluded samples")
    }
    cfg$TYPE[cfg$SHORT_LABEL %in% exclude]="E"
    cat("Excluding:",paste(exclude,collapse=","),"\n")
  }
  clones=cfg$LABEL[which(cfg$TYPE=="C")]
  snv=read.csv(sprintf("%s/%s_snp_pileup.tsv",data_dir,PREFIX), sep="\t", header=T, stringsAsFactors = F)
  colnames(snv)[1:4] = c("Chrom", "Pos", "Ref", "Alt")
  snv$Chrom=gsub("^chr","",snv$Chrom)
  #
  indel_file=sprintf("%s/%s_indel_pileup.tsv",data_dir,PREFIX)
  if(file.exists(indel_file)){
    indel=read.csv(indel_file, sep="\t", header=T, stringsAsFactors = F)
  }else{
    indel=snv[which(snv$Pos<0),]##Just get an empty indel file
  }
  colnames(indel)[1:4] <- c("Chrom", "Pos", "Ref", "Alt")
  indel$Chrom=gsub("^chr","",indel$Chrom)
  loci=sprintf("%s:%s:%s:%s",snv$Chrom,snv$Pos,snv$Ref,snv$Alt)


  projects=paste(unique(cfg$PROJECT[which(cfg$TYPE %in% c("C","E"))]),collapse=",")
  #Disable lane multiplexing QC - only sets mtr=1 to 0 - the downstream process is robust to this
  #lanes=load_lanes(config_dir,PREFIX,projects)
  inf=list(cfg=cfg,lanes=NULL,yamlf=yamlf)
  inf=add_meta_yaml(inf,PREFIX)
  #inf$lanes=inf$lanes[which(inf$lanes$NAME %in% inf$meta$clones),]
  samples=inf$cfg$LABEL
  
  
  mtrs=match(sprintf("%s_MTR",samples),colnames(snv))
  if( any(is.na(mtrs))){
    stop("Missing data for ",paste(samples[which(is.na(mtrs))],collapse=","))
  }
  wtrs=match(sprintf("%s_WTR",samples),colnames(snv))
  deps=match(sprintf("%s_DEP",samples),colnames(snv))
  ofs=match(sprintf("%s_OFS",samples),colnames(snv))
  mtrf=match(sprintf("%s_MTRF",samples),colnames(snv))

  driver_file=sprintf("%s/%s_%s_drivers.bed",data_dir,PREFIX,SPECIES_BUILD)
  if(file.exists(driver_file)){
    drivers=read.csv(driver_file, sep="\t", header=F, stringsAsFactors = F)
    colnames(drivers) = c("Chrom", "Pos", "Ref", "Alt")
    drivers$Chrom=gsub("^chr","",drivers$Chrom)
  
  }else{
    cat("WARNING:No driver file:",driver_file,"provided..\n")
    drivers=snv[which(snv$Pos<0),c("Chrom", "Pos", "Ref", "Alt")]#
    
  }
  drivers$key=with(drivers,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
  ##Typically non-existent now
  embfile=sprintf("%s/embryonic_mutations.txt",data_dir)
  
  if(file.exists(embfile)){
    cat("Reading predefined embryonic mutations from",embfile,"\n")
    embmuts=read.csv(embfile, sep="\t", header=T, stringsAsFactors = F)
    colnames(embmuts) = c("Chrom", "Pos", "Ref", "Alt","Patient")
    embmuts=embmuts[which(embmuts$Patient==PREFIX),]
    embmuts$key=with(embmuts,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
  }else{
    cat("Not using list of predefined embryonic mutatons\n")
    embmuts=list(key=c())
  }


  #SNV
  inf$snv$details=snv[,1:17]
  key=with(inf$snv$details,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
  inf$snv$details$keep_driver=ifelse(key %in% drivers$key,1,0)
  inf$snv$details$keep_embryonic=ifelse(key %in% embmuts$key,1,0)
  inf$snv$mtr=as.matrix(snv[,mtrs])
  inf$snv$wtr=as.matrix(snv[,wtrs])
  inf$snv$dep=as.matrix(snv[,deps])
  inf$snv$ofs=as.matrix(snv[,ofs])
  if(!any(is.na(mtrf))){
	inf$snv$mtrf=as.matrix(snv[,mtrf])
	}
  #Indel
  mtrs=match(sprintf("%s_MTR",samples),colnames(indel))
  wtrs=match(sprintf("%s_WTR",samples),colnames(indel))
  deps=match(sprintf("%s_DEP",samples),colnames(indel))
  ofs=match(sprintf("%s_OFS",samples),colnames(indel))
  mdr=match(sprintf("%s_MTRF",samples),colnames(indel))
  inf$indel$details=indel[,1:17]
  inf$indel$mtr=as.matrix(indel[,mtrs])
  inf$indel$wtr=as.matrix(indel[,wtrs])
  inf$indel$dep=as.matrix(indel[,deps])
  inf$indel$ofs=as.matrix(indel[,ofs])
  if(!any(is.na(mdr))){
	inf$indel$mdr=as.matrix(indel[,mdr])
  }
  key=with(inf$indel$details,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
  inf$indel$details$keep_driver=ifelse(key %in% drivers$key,1,0)
  inf$indel$details$keep_embryonic=ifelse(key %in% embmuts$key,1,0)
  for(x in c("indel","snv")){
    for(y in c("mtr","wtr","dep","ofs","mdr","mtrf")){
      #We use short labels to identify colonies - easier on plots
      if(!is.null(inf[[x]][[y]])){
        colnames(inf[[x]][[y]])=inf$cfg$SHORT_LABEL
      }
    }
  }
  ##Refilter to get rid of sites  only present in excluded samples
  if(length(exclude)>0){
    for(x in c("indel","snv")){
      #The "ofs" field stands for originally filter status. If it is 0 then the variant was called and passed filters (including LCM) otherwise not.
      idx.only.in.failed=which(rowSums(inf[[x]]$ofs[,inf$meta$clones_short]==0)==0  & inf[[x]]$details$keep_driver==0)
      if(length(idx.only.in.failed)>0){
        cat("Deleting",length(idx.only.in.failed),"of type",x," that are only present in excluded samples:",exclude,"\n")
        for(y in c("mtr","wtr","dep","ofs","details","mdr","mtrf")){
          if(!is.null(inf[[x]][[y]])){
            inf[[x]][[y]]=inf[[x]][[y]][-idx.only.in.failed,]
          }
        }
      }
    }
  }

  inf$patient=PREFIX
  ##Do a locus specific min variant VAF.
  ##Minimum mutataion VAF for clonal variant (minimum applies after event..).

  for(type in c("snv","indel")){
    dims=dim(inf[[type]]$mtr)
    inf[[type]]$minmutvaf=inf[[type]]$mtr
    inf[[type]]$minmutvaf[,]=0.5
    if(inf$meta$SEX=="M"){
      ##Won't work for build38..
      ##Identify XnonPAR and Y as ploidy=1.. (If necessary this is over-ridden by subsequent LOH and CNA processing)
      ## mouse GRCm38 169,969,759	170,931,299
      nonpar=read.table(sprintf("non_par_X_%s.bed",SPECIES_BUILD),head=F,stringsAsFactors = FALSE,sep="\t")
      #idx=with(inf[[type]]$details,which((Chrom=="X" & Pos>2699520 & Pos<154931044) | Chrom=="Y"))
      idx=with(inf[[type]]$details,which((Chrom=="X" & Pos>nonpar$V2 & Pos<nonpar$V3) | Chrom=="Y"))
      if(length(idx)>0){
        inf[[type]]$minmutvaf[idx,]=1
      }
    }
    inf[[type]]$germlinevaf=inf[[type]]$minmutvaf
    #CNA and LOH here are each lists with elements that are also a list which contain "chr","start","end","ploidy","samples" where the "samples"
    #field contains the list of samples sharing the copy number event.  This structure is manually curated via the YAML file.
    ##
    cnafile=sprintf("%s/%s.mscna",config_dir,PREFIX)
    lohfile=sprintf("%s/%s.msloh",config_dir,PREFIX)
    if(file.exists(cnafile)){
      cna=read.table(cnafile,head=TRUE,sep="\t")
      #cna$major=sapply(cna$major,function(x) as.integer(unlist(strsplit(x,split="-"))[2]))
      #cna$minor=sapply(cna$minor,function(x) as.integer(unlist(strsplit(x,split="-"))[2]))
      cna$ploidy=cna$major+cna$minor
      inf$meta$CNA=list()
      for(i in 1:nrow(cna)){
        row=cna[i,]
        inf$meta$CNA[[i]]=list(LABEL=row$label,
                               chr=row$chr,
                               start=row$start,
                               end=row$end,
                               ploidy=row$ploidy,
                               samples=unlist(strsplit(row$sample,",")))
      }
    }
    if(file.exists(lohfile)){
      cna=read.table(lohfile,head=TRUE,sep="\t")
      #cna$major=sapply(cna$major,function(x) as.integer(unlist(strsplit(x,split="-"))[2]))
      #cna$minor=sapply(cna$minor,function(x) as.integer(unlist(strsplit(x,split="-"))[2]))
      cna$ploidy=cna$major#+cna$minor
      inf$meta$LOH=list()
      for(i in 1:nrow(cna)){
        row=cna[i,]
        inf$meta$LOH[[i]]=list(LABEL=row$label,
                               chr=row$chr,
                               start=row$start,
                               end=row$end,
                               ploidy=row$ploidy,
                               samples=unlist(strsplit(row$sample,",")))
      }
    }
    for( x in c(inf$meta$LOH,inf$meta$CNA)){
      minmutvaf=1/x$ploidy
      idx=with(inf[[type]]$details,which(Chrom==x$chr & Pos>=x$start & Pos<x$end))
      samples=intersect(x$samples,inf$cfg$SHORT_LABEL)
      if(length(samples)>0  & length(idx)>0){
          inf[[type]]$minmutvaf[idx,samples]=ifelse(inf[[type]]$minmutvaf[idx,samples]>minmutvaf,minmutvaf,inf[[type]]$minmutvaf[idx,samples])
      }
    }
    #For downstream processing we keep the true minmutvaf for a clonal sample as minmutvaf2
    inf[[type]]$minmutvaf2=inf[[type]]$minmutvaf
    #and change minmutvaf to allow for an aberant cell fraction of 0.6 (typical VAF=0.3)
    inf[[type]]$minmutvaf=MINACF*inf[[type]]$minmutvaf
  }
  inf
}


##Sets mtr_d to NA where depth is strictly less than mindepth.
do_apply_depth_limit=function(inf,mindepth=6,b_verbose=FALSE){
  inf$snv$mtr_d=depth_limit(inf,inf$snv$details,inf$snv$mtr,inf$snv$dep,"snv",inf$meta$SEX,b_verbose,mindepth=mindepth)
  inf$indel$mtr_d=depth_limit(inf,inf$indel$details,inf$indel$mtr,inf$indel$dep,"indel",inf$meta$SEX,b_verbose,mindepth=mindepth)
  inf$snv$vaf_d=inf$snv$mtr_d/inf$snv$dep
  inf$indel$vaf_d=inf$indel$mtr_d/inf$indel$dep
  inf
}


depth_limit=function(inf,details,mtr,dep,label,sex,b_verbose,mindepth){
  ### Now filter on depth. NB treat autosomes and sex chromosomes differently.
  clones=get_clones(inf)
  n_sample=length(clones)
  vv2=mtr[,clones]
  dep=dep[,clones]
  #
  #MINDEPTH=mindepth
  if(sex=="M"){

    depthlimit=matrix(rep(ifelse(details$Chrom %in% c("X", "Y"),ceiling(mindepth/2),mindepth),n_sample),ncol=n_sample)
    idx.na=which(dep<depthlimit)
    nmut=dim(dep)[1]
    naut=length(which(!(details$Chrom[((idx.na-1) %% nmut)+1] %in% c("X", "Y"))))
    xaut=length(which((details$Chrom[((idx.na-1) %% nmut)+1] %in% c("X", "Y"))))
    if(b_verbose){
      cat("PROGRESS:",label,"Setting ",naut," autosomal sites to NA. REASON: Depth to low < ",mindepth,"\n")
      cat("PROGRESS:",label,"Setting ",xaut," X & Y sites to NA. REASON: Depth to low < ",ceiling(mindepth/2),"\n")
    }
  }else{
    depthlimit=mindepth
    idx.na=which(dep<depthlimit)
    if(b_verbose){
    cat("PROGRESS:",label,"Setting ",length(idx.na)," autosomal+X sites to NA for female sample. REASON: Depth to low < ",mindepth,"\n")
    }
  }
  vv2[idx.na]=NA
  mtr[,clones]=vv2
  mtr
}

get_clones=function(inf){
  if(length(intersect(inf$meta$clones_short,colnames(inf$snv$mtr)))>0){
    clones=inf$meta$clones_short
  }else{
    clones=inf$meta$clones
  }
}

get_normal=function(inf){
  if(length(intersect(inf$meta$clones_short,colnames(inf$snv$mtr)))>0){
    inf$cfg$SHORT_LABEL[match(inf$meta$normal,inf$cfg$LABEL)]
  }else{
    inf$meta$normal
  }
}


## Adds per locus summary cross-sample count information that are used by downstream filters.
##  geno    : binomial classifier based genotyping
##  g_counts : number of samples with mutant genotype
## g_nsamps  : number of samoples with non missing gentoypes
## g_nmiss   : number of samples with missing genotype
## counts_orig : number of samples with non-zero mutant reads
## ref_counts  : number of samples with zero mutant reads that have sufficient depth (encoded in missing status of vaf_d)
add_counts=function(inf){
  clones=get_clones(inf)
  for(x in c("snv","indel")){
    inf[[x]]$details$counts <- apply(inf[[x]]$vaf_d[,clones], 1, function(row) length(which(row>0)))
    geno=get_basic_binomial_genotypes(inf[[x]],clones)$geno
    inf[[x]]$details$g_counts=rowSums(geno==1,na.rm=TRUE)
    inf[[x]]$details$g_nsamps=rowSums(!is.na(geno),na.rm=TRUE)
    inf[[x]]$details$g_nmiss=rowSums(is.na(geno),na.rm=TRUE)
    inf[[x]]$geno=geno
    inf[[x]]$details$counts_orig <- apply(inf[[x]]$mtr[,clones], 1, function(row) length(which(row>0)))
    inf[[x]]$details$ref_counts  <- apply(inf[[x]]$vaf_d[,clones], 1, function(row) length(which(row==0))) ## Number with zero VAF
  }
  inf
}





##The following calls the functions that generate the "filter_<NAMED FILTER>" fields in "details"
run_filters=function(inf,b_basic_germ=TRUE){
  
  ##Applies mask to mtr count if depth<6 then set to NA.  masked mtr are stored in mtr_d which is used to calculate vaf_d.
  inf=do_apply_depth_limit(inf,mindepth=inf$params$min_cov)
  ##Excludes sites that have depth<6 (or min_cov) for more than 1/3 of the samples
  ##counts_orig is the number colonies with unfiltered mtr>0
  ##counts is the number colonies with depth-filtered mtr>0
  ##Ref count is the number of depth-filtered mtr that have mtr==0.
  inf=add_counts(inf)
  for(filter_name in inf$params$filters){
    cat("Running:",filter_name,":\n")
    fn=get(sprintf("do_%s_filter",filter_name))
    inf=fn(inf)
  }
  
  ##Filters variants within 10 bp of indel
  #inf=do_filter_near_indels(inf)
  
  
  ##Filters variants that are within 10 bp of each other
  #inf=do_filter_too_close(inf)
  
  ##Limited correction for lane hopping.
  ##Checks where samples a and b,c..d are sequenced in the same lane and where a has depth>8 and mtr=1 (1.8% chance if variant is genuine) and b,c,d has max(mtr)>5 then set mtr[a]=0.
  ##inf=do_apply_lane_correction(inf)
  
  ##Filters sites where counts==0 or counts==number of colonies or where no colony exhibits depth-filtered mtr>0
  ##Also filters based on genotype
  #inf=do_missing_filter(inf,maxmiss = max(5,floor(length(inf$meta$clones)/5)))
  
  ##Also does basic genotyping..
  #inf=do_count_filter(inf)
  
  ##Optional filter. Identifies sites are on different chromosome arms when lifted over from hg19 to hg38. This is useful identifying standalone LOH snps that are outside their region (9q)
  ##inf=do_19_to_38_change_arm_filter(inf)
  ##Germline filteres - identifies sites where the normal mtr count is much more likely due to a germline SNP than tumour contamination (assumed at 10%)
  ##For indels - classify sites as germline if VAF>20%.
  
  #inf=do_bglod_filter(inf)
  ##Complex VAF filters
  #inf=do_simple_vaf_filter(inf)
  #inf=do_beta_binomial_filter(inf)
  
  #if(is.na(inf$meta$normal) && b_basic_germ){
  #  inf=do_basic_germline_filter(inf)
  #}
  ##Identify CNA and LOH regions with localx flag..
  inf=do_cna_loh_exclusions(inf)
  inf
}

##The details matrix now contains the various filters "filter_<NAMED FILTER>"
##The following combines the filters to make the definitive FILTER variable
##Also recovers any embryonic mutations that might be lost via the bglod filter.
##Also checks for mutations that need kept for further analysis e.g. drivers.
##Any germline muts that creep through will be assigned to the germline root branch of the tree - so no worries.
do_filter=function(inf){
  #inf$filtered=list()
  if(is.na(inf$meta$normal)){
    exflag=c("filter_bglod")
  }else{
    exflag=c("filter_bglod","filter_basic_germline") ##,"filter_count")
  }
  EXCLUDE_ID=get_manual_exclusions()
  for(x in c("snv","indel")){
    if("filter_beta_binomial" %in% colnames(inf[[x]]$details)){
      # Form a temporary filter that is only fails when beta binomial fails and germline flags fail.
      inf[[x]]$details$filter_beta_binomial_ng=(apply(inf[[x]]$details[,exflag,drop=FALSE],1,sum,na.rm=TRUE)==0)*inf[[x]]$details$filter_beta_binomial
      exflag=union(exflag,"filter_beta_binomial")
    }
    
    filter_colidx=grep("^filter",colnames(inf[[x]]$details))
    filtered=apply(inf[[x]]$details[,filter_colidx],1,sum,na.rm=TRUE)>0
    #Find the filters that are not related to germline status
    filter_colidx_ex_glod=setdiff(filter_colidx,match(exflag,colnames(inf[[x]]$details)))
    inf[[x]]$details$FILTER_EXCEPT_GLOD=apply(inf[[x]]$details[,filter_colidx_ex_glod],1,sum,na.rm=TRUE)>0
    inf[[x]]$details$filter_beta_binomial_ng=NULL
    ##The following makes sure potential driver mutations/predefined embryonic mutations are retained
    if(x=="snv" && !is.na(inf$meta$normal)){
      ##Find embryonic muts..
      idxemb=identify_embryonic(inf)
      inf[[x]]$details$keep_embryonic=0
      if(length(idxemb)>0){
      inf[[x]]$details$keep_embryonic[idxemb]=1
      cat("To manually exclude write loci in the form <chr>:<pos> to",get_manual_exclude_file())
      }
    }

    keep=apply(inf[[x]]$details[,grep("^keep_",colnames(inf[[x]]$details))],1,sum,na.rm=TRUE)>0
    manual_exclude=with(inf[[x]]$details,ifelse(sprintf("%s:%s",Chrom,Pos) %in% EXCLUDE_ID,1,0))
    inf[[x]]$details$filter_manual_exclude=manual_exclude
    if(sum(manual_exclude,na.rm = TRUE)>0){
      cat("Manually excluding:\n")
      print(inf[[x]]$details[which(manual_exclude==1),])
    }
    #Recall: 
    #filtered is 1 if any filter_* are 1
    #keep is TRUE if any keep_* are 1
    inf[[x]]$details$FILTER=ifelse((filtered & !keep) | manual_exclude==1,1,0)
    inf[[x]]$filtered=list()
    idx=which(inf[[x]]$details$FILTER==0)
    cat(x,": Filtering down to ",length(idx),"variants\n")
    for(y in c("details","mtr","wtr","dep","mtr_d","vaf_d","ofs","minmutvaf","minmutvaf2")){
      inf[[x]]$filtered[[y]]=inf[[x]][[y]][idx,]
      if(y != "details"){
        colnames(inf[[x]]$filtered[[y]])=inf$cfg$SHORT_LABEL
      }
    }
  }
  inf
}


do_genotyping_basic_binomial=function(inf,factor=20,p.error=0.01){
  clones=inf$meta$clones_short
  for(x in c("snv","indel")){
    res=get_basic_binomial_genotypes(inf[[x]]$filtered,clones,factor=factor,p.error=p.error)
    inf[[x]]$filtered$geno_colony=res$geno #geno
    inf[[x]]$filtered$gp=res$gp
    colnames(inf[[x]]$filtered$geno_colony)=clones
  }
  inf
}

get_basic_binomial_genotypes=function(inf,clones,factor=20,p.error=0.01){
  #clones=inf$meta$clones_short
  #pmin=p.error
  #p.error is probability of erroneous specified mutant read for reference - so epsilon/3 = p.error  
  epsilon=3*p.error
  mtr=inf$mtr[,clones]
  ref=inf$dep[,clones]-mtr
  expected_mutvaf=inf$minmutvaf2[,clones]
  #bf=mtr*log(0.5)+ref*log(0.5)-(mtr*log(pmin)+ref*log(1-pmin))
  #p(mut)=p(mut|mutant read)p(mutant read)+p(mut|ref read)p(ref read) \]
  #p(mut)=(1-epsilon)mu+(epsilon/3)*(1-mu)
  #p(mut)=(1-(4*epsilon)/3)*mu + epsilon/3
  pmm=(1-4*epsilon/3)*expected_mutvaf+epsilon/3
  pmr=epsilon/3
  bf=mtr*log(pmm)+ref*(log(1-pmm))-(mtr*log(pmr)+ref*log(1-pmr))
  #bf=mtr*log(expected_mutvaf-pmin)+ref*log(1-expected_mutvaf+pmin)-(mtr*log(pmin)+ref*log(1-pmin))
  ##Require classifier is <factor> times as likely as other option to get non-missing genotype.
  geno=as.matrix(ifelse(abs(bf)>log(factor),ifelse(bf>0,1,0),NA),ncol=dim(bf)[2])
  gp=exp(bf)/(1+exp(bf))
  list(geno=geno,gp=gp)
}

#Ensure we set variants to NA if the geno is 0 as could just be lost
do_loh=function(inf){
  for(loh in inf$meta$LOH){
    samps=intersect(loh$samples,inf$meta$clones_short)
    for(x in c("snv","indel")){
      changed=0
      flag=with(data = inf[[x]]$filtered$details,(Chrom==loh$chr) & Pos>loh$start & Pos<loh$end)
      ##Exclude potential germline variants that are now homozygous reference.
      for(samp in samps){
        idx=which(flag & (inf[[x]]$filtered$geno_colony[,samp]==0 | inf[[x]]$filtered$mtr[,samp]==0))
        inf[[x]]$filtered$geno_colony[idx,samp]=NA
        changed=changed+length(idx)
        cat("LOH: NA'd ",length(idx)," homozygous reference genotypes; type = ",x," in ",samp,"\n")
      }
      cat("LOH: NA'd ",changed," homozygous reference genotypes; type = ",x,"\n")
    }
  }
  
  inf
}
#Set exclusion flags
do_cna_loh_exclusions=function(inf){
  for(x in c("snv","indel")){
    if(nrow(inf[[x]]$details)>0){
      inf[[x]]$details$is_localx_excluded=0
      for(event in c(inf$meta$LOH,inf$meta$CNA)){
        samps=intersect(event$samples,inf$meta$clones_short)
        if(length(samps)>0){
          idx=with(data = inf[[x]]$details,(Chrom==event$chr) & Pos>=event$start & Pos<=event$end)
          if(length(idx)>0){
            inf[[x]]$details$is_localx_excluded[idx]=1
          }
        }
      }
      ##Exclude X for males and females.
      idx=which(inf[[x]]$details$Chrom %in% c("X","Y"))
      if(length(idx)>0){
        inf[[x]]$details$is_localx_excluded[idx]=1
      }
    }else{
      inf[[x]]$details$is_localx_excluded=c()
    }
  }
  inf
}



merge_snps_and_indels=function(inf){
  fields=colnames(inf$snv$filtered$details)[which(colnames(inf$snv$filtered$details) %in% colnames(inf$indel$filtered$details))]

  #details=rbind(inf$snv$filtered$details[,1:17],inf$indel$filtered$details[,1:17])
  details=rbind(inf$snv$filtered$details[,fields],inf$indel$filtered$details[,fields])
  chr=match(details$Chrom,c(seq(1,22),"X","Y"))
  idx=order(chr,details$Pos,details$Ref,details$Alt)
  dat=list()
  dat$details=details[idx,]
  for(x in c("mtr","wtr","dep","mtr_d","geno_colony","gp","ofs","minmutvaf")){
    ##check same size
    if(dim(inf$snv$filtered[[x]])[1]+dim(inf$indel$filtered[[x]])[1]!=length(idx)){
      stop("Dimension mimsatch in merge_snps_and_indels!")
    }
    dat[[x]]=rbind(inf$snv$filtered[[x]],inf$indel$filtered[[x]])[idx,]
  }
  dat$details$TYPE=with(dat$details,ifelse(nchar(Ref)+nchar(Alt)==2,"SNV","INDEL"))
  geno=cbind(dat$geno_colony,zeros=0)
  list(dat=dat,geno=geno)
}



plot_vaf_genotyped=function(inf,type){
  clones=inf$meta$clones_short
  with(inf[[type]]$filtered,
            boxplot(lapply(clones,function(x) vaf_d[which(geno_colony[,x]==1),x]),ylim=c(0,1),names=clones,las=2,notch=T)
  )
  abline(h=0.5,col="blue")
  title(sprintf("%s: %s : VAF Distribution at sites where genotype=1",inf$patient,toupper(type)))
}

get_exclude_file=function(){
  sprintf("%s/%s_suggested_exclude.txt",output_dir,prefix)
}

get_tree_based_exclude_file=function(){
  sprintf("%s/%s_treebased_exclude.txt",output_dir,prefix)
}

get_manual_exclude_file=function(){
  sprintf("%s/%s_manual_exclude.txt",output_dir,prefix)
}
get_manual_exclusions=function(){
if(file.exists(get_manual_exclude_file())){
  EXCLUDE_ID=read.table(get_manual_exclude_file(),header = F,stringsAsFactors = FALSE)$V1
}else{
  EXCLUDE_ID=c()
}
  EXCLUDE_ID
}

##Plots VAF but has the side effect of also writing suggested excluded samples to file when type=="snv" and germline==FALSE
## TODO - move the exclusion code into a more appropriately named function!
plot_vaf_ofs=function(inf,type,germline=FALSE,method="box",vaf_threshold=0.5*MINACF,mindepth=6,b_only_called=TRUE){
  if(dim(inf[[type]]$details)[1]==0){
    return(NULL)
  }
  clones=inf$meta$clones_short
  if(germline){
    idx=with(inf[[type]]$details,which(filter_bglod==1 & FILTER_EXCEPT_GLOD==0))
    if(b_only_called){
      vt="Germline [Originally Called]"
    }else{
      vt="Germline [Based on the inferred germline sites]"
    }
  }else{
    idx=with(inf[[type]]$details,which(FILTER==0))
    vt="Somatic (pass all filters)"
    if(!b_only_called){
      stop("Unsupported combo for somatic variants")
    }
  }
  ##The following picks out just those sites/sample that passed in caveman/pindel
  ##Focus on autosomes..
  dat=with(inf[[type]],lapply(clones,function(x){
    if(b_only_called){
      idx2=which(ofs[idx,x]==0 & details[["is_localx_excluded"]][idx] ==0)
    }else{
      idx2=which(details[["Chrom"]][idx] %in% 1:22)
    }
    v=vaf_d[idx[idx2],x];v[which(!is.na(v))]})
  )
  ##For each sample check median depth> threshold and that sample is sufficiently clonal (MINACF*0.5)
  pvaf=with(inf[[type]],sapply(clones,function(x){
    idx2=which(ofs[idx,x]==0 & details[["is_localx_excluded"]][idx] ==0)
    m=mtr[idx[idx2],x]
    d=dep[idx[idx2],x]
    idxx=which(d<3*median(d))
    if(length(idxx)<50 && type=="snv" && !germline){
      #browser()
      cat("Too few variants in",x,"\n")
      return(-1)
    }
    pbinom(sum(m[idxx]),sum(d[idxx]),vaf_threshold)
    })
  )

  idx.vaf.fail=which(pvaf<0.01)

  #idxx=with(inf[[type]]$details,which(filter_bglod==1 & FILTER_EXCEPT_GLOD==0 & is_localx_excluded==0))
  idxx=with(inf[[type]]$details,which(FILTER==0 & is_localx_excluded==0))
  dep=with(inf[[type]],sapply(clones,function(x){
    v=dep[idxx,x];median(v)})
  )
  if(type=="snv" && germline==FALSE && method!="fails"){

    if(length(idx.vaf.fail)>0){
      cat("The following have VAF significantly<",vaf_threshold,":", paste(clones[idx.vaf.fail],collapse=","),"\n")
    }
    idx.depth.fail=which(dep<mindepth)
    if(length(idx.depth.fail)>0){
      #cat("The following have median depth< ",mindepth," at good quality germline sites:", paste(clones[idx.depth.fail],collapse=","),"\n")
      cat("The following have median depth< ",mindepth," at PASS sites:", paste(clones[idx.depth.fail],collapse=","),"\n")
      
    }
    ##Check germline sensitivity
    df=check_sensitivity(inf,type)
    low_sensitivity_colonies=df$sample_short[which(df$sensitivity<inf$meta$MIN_SENSITIVITY)]
    if(length(low_sensitivity_colonies)>0){
      idx.sensitivity.fail=match(low_sensitivity_colonies,clones)
      if(any(is.na(idx.sensitivity.fail))){
        cat("WARNING: Missing lookup of low sensitivity clones!\n")
      }
    }else{
      idx.sensitivity.fail=c()
    }
    idx.colony.fail=sort(unique(c(idx.vaf.fail,idx.depth.fail,idx.sensitivity.fail)))
    #if(length(idx5)>0){
      cat("Suggest the following are removed because of low VAF or depth or senstivity:", paste(clones[idx.colony.fail],collapse=","),"\n")
      newexclude=sort(unique(c(clones[idx.colony.fail],with(inf$cfg,SHORT_LABEL[which(TYPE=="E")]))))
      writeLines(paste(newexclude,collapse=","),con=get_exclude_file())
    #}
      write.table(df,file=paste0(output_dir,"/",prefix,"_sample_sensitivity.txt"),quote = FALSE,row.names = FALSE,col.names = TRUE,sep="\t")
    #write.table(df,)

  }
  if(method=="box"){
    if(length(clones)>50){
      cex.sz=50/length(clones)
    }else{
      cex.sz=1
    }
    boxplot(dat,ylim=c(0,1),names=clones,las=2,notch=T,par=list(ylim=c(0,1)),col=ifelse(pvaf>0.01,"white",ifelse(pvaf<0,"grey","red")),border=ifelse(dep<8,"red","black"),cex.axis=cex.sz)
    #points(rep(1:length(dat),sapply(dat,length))+0.1*(runif(length(unlist(dat)))-0.5),unlist(dat),pch=15,cex=0.5,col="red")
  }else if(method=="vio"){
    dat$names=rep("",length(clones))
    do.call("vioplot",dat)
    axis(side=1,at=1:length(clones),labels=clones,las=2)
  }else if(method=="fails"){
    idx=with(inf[[type]]$details,which(FILTER==1 & filter_bglod==0))
    vt="non-germline FAIL "
    dat=with(inf[[type]],lapply(clones,function(x){
      v=vaf_d[idx[which(ofs[idx,x]==0)],x];v[which(!is.na(v))]})
    )
    boxplot(dat,ylim=c(0,1),names=clones,las=2,notch=T)
  }else{
    stop("unsupported method supplied. try box|vioplot|counts")
  }
  abline(h=0.5,col="blue")
  algo=ifelse(type=="snv","Caveman","Pindel")
  title(sprintf("%s: %s : VAF Distribution at %s sites called in %s",inf$patient,toupper(type),vt,algo))
}

plot_filters=function(inf,type){
  if(dim(inf[[type]]$details)[1]==0){
    return(NULL)
  }
  par(mar=c(8,4,4,2)+0.1)
  fields=colnames(inf[[type]]$details)[grep("^filter_",colnames(inf[[type]]$details))]
  counts=apply(inf[[type]]$details[,fields],2,sum,na.rm=T)
  idx=which(counts>0)
  ym=10**ceiling(log10(max(counts)))
  barplot(counts[idx],names=gsub("filter_","",fields)[idx],las=2,log = "y",ylim=c(1,ym),yaxt="n")
  axis(side = 2,at = 10**seq(0,log10(ym)),labels=sprintf("%d",10**seq(0,log10(ym))),las=2)
  title(sprintf("%s: %s: No. Sites removed by each filter",inf$patient,toupper(type)))
  data.frame(FILTER=fields,COUNT=counts)
}

plot_unique_filters=function(inf,type){
  if(dim(inf[[type]]$details)[1]==0){
    return(NULL)
  }
  par(mar=c(8,4,4,2)+0.1)
  fields=colnames(inf[[type]]$details)[grep("^filter_",colnames(inf[[type]]$details))]
  idx=which(apply(inf[[type]]$details[,fields],1,sum,na.rm=T)==1)
  counts=apply(inf[[type]]$details[idx,fields],2,sum,na.rm=T)
  idx=which(counts>0)
  ym=10**ceiling(log10(max(counts)))
  barplot(counts[idx],names=gsub("filter_","",fields)[idx],las=2,log = "y",yaxt="n")
  axis(side = 2,at = 10**seq(0,log10(ym)),labels=sprintf("%d",10**seq(0,log10(ym))),las=2)
  title(sprintf("%s: %s: No. Sites UNIQUELY removed by each filter",inf$patient,toupper(type)))
  data.frame(FILTER=fields,COUNT=counts)
}


findrho=function(nmuts,depth){
  pseudo=1e-6
  loglik=function(par,nmuts,depth){
    idx=which(!is.na(nmuts/depth))
    sum(VGAM::dbetabinom(x =nmuts[idx],size = depth[idx],prob = par[2],rho = par[1],log = T))
  }
  optres=optim(par=c(0.1,0.1),loglik, gr = NULL,method="L-BFGS-B",lower =c(pseudo,pseudo),upper=c(0.8,0.9),control=list(fnscale=-1),nmuts=nmuts,depth=depth)
  #c(rho=optres$par[1],p=optres$par[2])
  optres$par[1]
}

apply_rho=function(inf,type){
  clones=get_clones(inf)
  if(!is.null(inf$params$ncores)){
    MC_CORES=inf$params$ncores
  }
  cat("PROGRESS: Calculating beta-binomial rho across",length(clones)," colonies at",dim(inf[[type]]$details)[1],"loci using ",MC_CORES," cores [Might take a while..]\n")
  inf[[type]]$details$rho=unlist(mclapply(1:dim(inf[[type]]$details)[1],function(i) {
    #findrho(as.numeric(inf[[type]]$mtr[i, clones]),as.numeric(inf[[type]]$dep[i, clones]))
    estimateRho_gridml(as.numeric(inf[[type]]$mtr[i, clones]),as.numeric(inf[[type]]$dep[i, clones])) 
  },mc.cores=MC_CORES))
  inf
}
estimateRho_gridml = function(NV_vec,NR_vec) {
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(NV_vec)/sum(NR_vec)
  ll = sapply(rhovec, function(rhoj) sum(VGAM::dbetabinom(x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=T)))
  return(rhovec[ll==max(ll)][1])
}

get_filter_columns=function(inf){
  inf$snv$details %>% (function(x) colnames(x)[grep("^filter_",colnames(x))])
}
identify_embryonic=function(inf,exclude=c(),normal_depth_cutoff=20){
  include=setdiff(inf$meta$clones_short,exclude)
  normal=get_normal(inf)#inf$cfg$SHORT_LABEL[which(inf$cfg$TYPE=="N")][1]
  if(is.na(normal)){
    cat("Non normal specified.. No need to find embryonic..")
    return(c())
  }
  ##Exclude LOH and CNA regions..
  ##
  lidx=c()
  for(loh in c(inf$meta$LOH,inf$meta$CNA)){
    idx=with(inf$snv$details,which(Chrom==loh$chr & Pos >= loh$start & Pos <= loh$end))
    lidx=c(lidx,idx)
  }
  lidx=sort(unique(lidx))

  ##We require the embryonic variant to be in the normal sample (BGLOD=1) [if BGLOD=0 it will have been passed filteres anyway] to not be a known SNP (GSNP=0) and to pass other filters..
  # Require a global VAF signicantly < 0.45
  #Then require that the variant exhibit beta-binomial overdispersion >0.05
  #We subsequently manually check whether the variants fit the tree.
  filters=get_filter_columns(inf)
  filters=setdiff(filters,sprintf("filter_%s",c("bglod","count","basic_germline")))  ##  Will
  
  idx=with(inf$snv$details %>% mutate(pass_other_filters=rowSums(inf$snv$details[,filters])==0),
           which(BGLOD==1 & GSNP==0 & pass_other_filters & inf$snv$dep[,normal]>normal_depth_cutoff))

  cat("Found",length(idx),"initial candidates\n")
  idx=setdiff(idx,lidx)
  cat("After eliminating LOH and CNA regions we have ",length(idx)," candidates\n")

  cat("Normal sample is ",normal,"\n")
  mtr=rowSums(inf$snv$mtr[idx,include],na.rm = TRUE)
  dep=rowSums(inf$snv$dep[idx,include],na.rm = TRUE)
  ##Look for VAF significantly less than MINACF*0.5
  pval=pbinom(mtr,dep,MINACF*0.5)

  ##Bonferonni correct..
  idx2=which(pval<0.01)
  if(length(idx2)>0){
    include2=intersect(include,colnames(inf$snv$geno))

    cat("found ",length(idx2),"candidate variants.  Checking beta-binomial measure of lumpiness...\n")
    rho=sapply(idx[idx2],function(i) findrho(inf$snv$mtr[i,include],inf$snv$dep[i,include]))

    idx3=which(rho>0.05)###min(rhocutoff,0.05))
    if(length(idx3)>0){
      cat("Retaining...",length(idx3)," potential embryonic variants")
    }
    ridx=idx[idx2][idx3]
    if(length(ridx)>1){
      profile=apply(inf$snv$geno[ridx,include2],1,function(x) paste(ifelse(is.na(x),"?",x),collapse=""))
    }else{
      if(length(ridx)==1){
      profile=paste(ifelse(is.na(inf$snv$geno[ridx,include2]),"?",inf$snv$geno[ridx,include2]),collapse="")
      }else{
        return(c())
      }
    }
    tmp=cbind(inf$snv$details[ridx,c("Chrom", "Pos","GENE","GLOD",colnames(inf$snv$details)[grep("^filter",colnames(inf$snv$details))])],profile=profile)
    
    print(tmp)
    ##write.table()
    ridx
  }else{
    c()
  }
}


get_cmcalls=function(dat,wild.type.colonies,patient,ymax=NA,type="snv",exclude.germline=T,b.do.barplot=TRUE){
  clones=dat$inf$meta$clones_short
  #browser()
  if(exclude.germline){
    inf=dat$inf[[type]]$filtered
    idx2=which(inf$details$filter_bglod==0 | inf$details$filter_count==0)
  }else{
    inf=dat$inf[[type]]
    idx2=1:dim(inf$details)[1]
  }
  if(dim(inf$details)[1]==0){
    return(NULL)
  }
  
  cnt=apply(1-inf$ofs[idx2,clones],2,sum)
  depth=apply(inf$dep[idx2,clones],2,mean)
  vtype=ifelse(type=="snv","Caveman","Pindel")
  #
  idx=order(depth)
  mt=rep(TRUE,length(cnt))
  mcnt=max(cnt)
  if(is.na(ymax)){
    if(type=="snv"){
      ymax=1000*ceiling(1.3*mcnt/1000)
    }else{
      ymax=100*ceiling(1.3*mcnt/100)
    }
  }
  #browser()
  mt[clones %in% wild.type.colonies]=FALSE
  if(b.do.barplot){
    par(mar=c(5,7,4,5)+0.1)
    xm=barplot(cnt[idx],las=2,col=ifelse(mt[idx],"grey","red"),xlab="Depth",main=sprintf("%s:No. %s Calls vs Sample Type and Depth",patient,vtype),ylim=c(0,ymax))
    MD=5* ((max(depth) %/% 5)+1)
    ym=par("usr")[4]
    lines(xm,ym*depth[idx]/MD,type="b",pch=19,col="black")
    axis(side=4,at=seq(0,MD,5)*par("usr")[4]/MD,labels=seq(0,MD,5),las=2)
    mtext(side=4,line=2,text="Depth",las=2)
    mtext(side=2,line=4,text="#Calls",las=2)
    if(length(wild.type.colonies)>0){
      legend("topleft",c("mt","wt"),col=c("grey","red"),pch=15,pt.cex=3)
    }else{
      legend("topleft",c("#Calls"),col=c("grey"),pch=15,pt.cex=3)
    }
    legend("top",c("Depth"),pch=19,col="black",pt.bg = "white", lty = 1, bty="n")
  }
  df=data.frame(patient=patient,samples=clones[idx],count=cnt[idx],depth=depth[idx],is_mutant=mt[idx])
  
  with(df,plot(count~depth,cex=0.5,pch=19,col="red",main=sprintf("Called %s Count vs Average Depth",toupper(type))))
  with(df,text(depth,count,labels=samples))
  df
}


plot_heatmaps=function(inf,maxrows=1000,minNonNA=10){
  minNonNA=min(minNonNA,length(inf$meta$clones_short)-2)
  s=inf$meta$clones_short
  idx=which(inf$snv$details$FILTER==0)
  vaf=inf$snv$vaf_d[idx,s]
  rownames(vaf)=ifelse(inf$snv$details$keep_driver[idx]==1,inf$snv$details$GENE[idx],"")
  snv_pass=subsample(vaf,maxrows,minNonNA)
  snv_fail_noglod=subsample(inf$snv$vaf_d[which(inf$snv$details$FILTER_EXCEPT_GLOD==1 & inf$snv$details$filter_bglod==0),s],maxrows,minNonNA)
  idx=which(inf$indel$details$FILTER==0)
  vaf=inf$indel$vaf_d[idx,s]
  rownames(vaf)=ifelse(inf$indel$details$keep_driver[idx]==1,inf$indel$details$GENE[idx],"")
  indel_pass=subsample(vaf,maxrows,minNonNA)
  ##
  idx=which(inf$indel$details$FILTER_EXCEPT_GLOD==1 & inf$indel$details$filter_bglod==0)
  vaf=inf$indel$vaf_d[idx,s]
  rownames(vaf)=inf$indel$details$GENE[idx]#,c("GENE","Pos")],1,paste,collapse=":")

  indel_fail_noglod=subsample(vaf,maxrows,minNonNA)
  doheatmap(snv_pass,"SNV PASS")
  doheatmap(snv_fail_noglod,"SNV FAIL (Not Germline)")
  if(dim(inf$indel$details)[1]==0){
    return(NULL)
  }
  doheatmap(indel_pass,"INDEL PASS")

  doheatmap(indel_fail_noglod,"INDEL FAIL (Not Germline)")

}

subsample=function(mat,maxrows,minNonNA){

  mat=mat[which(rowSums(!is.na(mat))>=minNonNA),]
  N=dim(mat)[1]
  idx2=which(rownames(mat)!="")
  if(length(idx2)>20){
    idx2=NULL
  }

  if(N>maxrows){
    mat[unique(c(sample(N,maxrows),idx2)),]
  }else{
    mat
  }


}

doheatmap=function(vaf,title_txt){
  distmat=as.matrix(dist(vaf))
  idx=which(rowSums(is.na(distmat))>0)
  if(length(idx)>0){
    cat("removing",length(idx)," mutations from matrix for heatmap\n")
    vaf=vaf[-idx,]
  }
  cat("VAF dim=",dim(vaf),"\n")


  try(pheatmap(vaf,breaks=c(-0.01,0.1,0.2,0.3,0.8,1.01),
               col=c("white","lightgrey","darkgrey","black","green"),
               na_col="red",main=title_txt,fontsize_row = 5))
}

get_overlap=function(inf,type){
  filters=colnames(inf[[type]]$details)[grep("^filter",colnames(inf[[type]]$details))]
  F=as.matrix(inf[[type]]$details[,filters])
  cnames=gsub("^filter_","",colnames(F))
  colnames(F)=cnames
  print(data.frame(filter=filters,NA_COUNT=colSums(is.na(F))))
  F[is.na(F)]=0
  count_matrix=t(F) %*% F
  total=length(which(rowSums(F)>0))
  prop_matrix=count_matrix/total
  list(count_matrix=count_matrix,prop_matrix=prop_matrix,total=total,pass=length(which(rowSums(F)==0)))
}
plot_table=function(tab,fmt="%s",title="Mutation Overlap",cex.headings=0.8){
  plot(NULL,xlim=c(0,1),ylim=c(0,1),axes = FALSE,xlab="",ylab="",main=title)
  R=dim(tab)[1]
  C=dim(tab)[2]
  dr=1/(R+2)
  dc=1/(C+2)
  sdr=sort(seq(dr,by = dr,length.out = R+1),decreasing = TRUE)
  sdc=seq(dc,by=dc,length.out = C+1)
  text(x=sdc[1],y=sdr[-1],labels=colnames(tab),cex=cex.headings)
  text(x=sdc[-1],y=sdr[1],labels=rownames(tab),cex=cex.headings)
  for(i in 1:R){
    font=rep(1,C)
    font[i]=2
    text(x=sdc[-1],y=sdr[i+1],labels=sprintf(fmt,tab[i,]),font=font)
  }
  arrows(x0=sdc+dc/2,y0=0,y1=1,lend=2,length=0)
  arrows(y0=sdr+dr/2,x0=0,x1=1,lend=2,length=0)
}

plot_too_few_warning=function(type="indels"){
  plot(NULL,xlim=c(0,1),ylim=c(0,1));
  text(0.1,0.5,sprintf("To few %s for plots!",type))
}
do_qc_plots=function(inf,output_dir){
  prefix=inf$meta$prefix
  pdf(sprintf("%s/%s_qc_plots.pdf",output_dir,prefix),w=14,h=10)
  par(mfrow=c(2,1),oma = c(0, 0, 5, 0))
  overlap.snv=get_overlap(inf,"snv")
  overlap.indel=get_overlap(inf,"indel")
  plot_table(100*overlap.snv$prop_matrix,fmt="%3.1f",title=sprintf("SNV:Pairwise Filter Overlap (%% Filtered Variants):Total Fail=%s/Pass=%s",overlap.snv$total,overlap.snv$pass),cex.headings = 0.8)
  plot_table(100*overlap.indel$prop_matrix,fmt="%3.1f",title=sprintf("Indels:Pairwise Filter Overlap (%% Filtered Variants):Total Fail=%s/Pass=%s",overlap.indel$total,overlap.indel$pass))
  mtext(prefix, outer = TRUE, cex = 1.5)


  dfds=get_cmcalls(dat=list(inf=inf),NULL,prefix,type="snv")
  dfdi=get_cmcalls(dat=list(inf=inf),NULL,prefix,type="indel")
  #abline(h=6,col="red")
  #par(mfcol=c(1,1))
  #plot_vaf_fpr(inf,"snv")
  #plot_vaf_fpr(inf,"indel")
  par(mfrow=c(2,1),oma = c(0, 0, 5, 0))
  plot_vaf_ofs(inf,"snv",germline=F,method="box")
  plot_vaf_ofs(inf,"snv",germline=T,method="box",b_only_called = FALSE)
  plot_vaf_ofs(inf,"snv",germline=T,method="box",b_only_called = TRUE)
  plot_vaf_ofs(inf,"snv",germline=F,method="fails")
  plot_vaf_density(inf,"snv")
  if(nrow(inf[["indel"]]$details)>10){
    par(mfrow=c(2,1),oma = c(0, 0, 5, 0))
    plot_vaf_ofs(inf,"indel",germline=F,method="box")
    plot_vaf_ofs(inf,"indel",germline=T,method="box",b_only_called = FALSE)
    plot_vaf_ofs(inf,"indel",germline=T,method="box",b_only_called = TRUE)
    plot_vaf_ofs(inf,"indel",germline=F,method="fails")
    plot_vaf_density(inf,"indel")
  }else{
    plot_too_few_warning()
  }

  #plot_vaf_density
  dfs=plot_filters(inf,"snv")
  dfsu=plot_unique_filters(inf,"snv")
  if(nrow(inf[["indel"]]$details)==0){
    plot_too_few_warning()
    dfi=head(dfs,0)
    dfiu=head(dfsu,0)
  }else{
    dfi=plot_filters(inf,"indel")
    dfiu=plot_unique_filters(inf,"indel")
  }
  #dfe=plot_extra(inf)
  plot_heatmaps(inf)
  dev.off()
  
  write.table(dfs,file =sprintf("%s/%s_filter_snv.txt",output_dir,prefix) ,sep = "\t",row.names=FALSE,col.names = TRUE,
              quote=FALSE)
  write.table(dfsu,file =sprintf("%s/%s_filter_unique_snv.txt",output_dir,prefix) ,sep = "\t",row.names=FALSE,col.names = TRUE,
              quote=FALSE)
  write.table(dfi,file =sprintf("%s/%s_filter_indel.txt",output_dir,prefix) ,sep = "\t",row.names=FALSE,col.names = TRUE,
              quote=FALSE)
  write.table(dfiu,file =sprintf("%s/%s_filter_unique_indel.txt",output_dir,prefix) ,sep = "\t",row.names=FALSE,col.names = TRUE,
              quote=FALSE)
  #write.table(dfe,file =sprintf("%s/%s_filter_extra.txt",output_dir,prefix) ,sep = "\t",row.names=FALSE,col.names = TRUE,
  #            quote=FALSE)
  write.table(dfds,file =sprintf("%s/%s_caveman_counts_vs_depth.txt",output_dir,prefix) ,sep = "\t",row.names=FALSE,col.names = TRUE,
              quote=FALSE)
  write.table(dfdi,file =sprintf("%s/%s_pindel_counts_vs_depth.txt",output_dir,prefix) ,sep = "\t",row.names=FALSE,col.names = TRUE,
              quote=FALSE)
}

check_sensitivity=function(inf,type){
  idx=get_germline_sites(inf,type)
  samples=inf$meta$clones_short
  #browser()
  depth=apply(inf[[type]]$dep[idx,samples,drop=FALSE],2,median,na.rm=TRUE)
  sensitivity=apply(1-inf[[type]]$ofs[idx,samples,drop=FALSE],2,mean,na.rm=TRUE)
  data.frame(patient=inf$meta$prefix,
             sample_short=samples,sample=inf$meta$clones,
             lcm=inf$cfg$LCM[match(samples,inf$cfg$SHORT_LABEL)],
             depth=depth,
             sensitivity=sensitivity
             ,N=length(idx),stringsAsFactors=FALSE)

}

##Functions.
get_germline_sites=function(inf,type="snv"){
  details=inf[[type]]$details
  filters=setdiff(colnames(details)[grep("^filter",colnames(details))],c("filter_bglod","filter_count","filter_basic_germline","filter_beta_binomial"))
  normal=inf$cfg$SHORT_LABEL[match(inf$meta$normal,inf$cfg$LABEL)]
  if(is.na(normal)){
    idx=which(rowSums(details[,filters])==0 & details$filter_bglod==1 & details$Chrom %in% 1:22)
  }else{
    idx=which(rowSums(details[,filters])==0 & details$filter_bglod==1 & inf[[type]]$mtr[,normal]>5  & details$Chrom %in% 1:22)##Need at least 5 mtr reads in normal
  }
  cat("Found",inf$meta$prefix,": Found",length(idx),"germline loci\n")
  ##Now exclude copy number regions
  idx=setdiff(idx,do.call("c",lapply(c(inf$meta$LOH,inf$meta$CNA),
                                     function(x) with(inf[[type]]$details,which(Chrom==x$chr & Pos>=x$start & Pos<x$end))
  )
  ))
  cat("Found",inf$meta$prefix,": Found",length(idx),"germline loci [excluding LOH]")
  idx
}

gather_vaf_unfiltered=function(inf,clones=inf$cfg$SHORT_LABEL[which(inf$cfg$TYPE %in% c("C","E"))],type="snv"){
  clones=clones## Force evaluation of default argument!
  patient=inf$patient
  cfg=inf$cfg
  inf=inf[[type]]
  ddim=dim(inf$mtr)
  df=do.call("rbind",lapply(clones,function(clone){
    #recall is_localx_excluded screens LOH/CNA + Sex Chromosomes
    idx=with(inf,which(ofs[,clone]==0 &  details$filter_bglod==0 &  details$is_localx_excluded==0))
    N=length(idx)
    data.frame(patient=rep(patient,N),
               sample=rep(clone,N),mtr=inf$mtr[idx,clone],depth=inf$dep[idx,clone],filter=inf$details$FILTER[idx],stringsAsFactors = TRUE)
  }))
  df$VAF=df$mtr/df$dep
  df$TYPE=cfg$TYPE[match(df$sample,cfg$SHORT_LABEL)]
  df$TYPE=ifelse(df$TYPE=="C","PASS","FAIL")
  df
}

plot_vaf_density=function(inf,type){
  dfvaf=gather_vaf_unfiltered(inf,type=type)
  typecolScale=RColorBrewer::brewer.pal(3,"Set1")[2:1]
  names(typecolScale)=c("PASS","FAIL")
  typecolScale=scale_color_manual(name="TYPE",values=typecolScale)
  p1=ggplot(data=dfvaf %>% filter(depth>0),aes(x=VAF,col=TYPE))+geom_density()+geom_vline(xintercept=0.5,linetype = "longdash")+xlim(c(0,1))+
    typecolScale+facet_wrap(~sample)+ggtitle(sprintf("%s:Unfiltered VAF distributions for %s: Includes excluded samples",toupper(type),inf$patient))+theme_bw()
  print(p1)
  p1=ggplot(data=dfvaf %>% filter(depth>0,filter==0),aes(x=VAF,col=TYPE))+geom_density()+geom_vline(xintercept=0.5,linetype = "longdash")+xlim(c(0,1))+
    typecolScale+facet_wrap(~sample)+ggtitle(sprintf("%s:Filtered VAF distributions for %s: Includes excluded samples",toupper(type),inf$patient))+theme_bw()
  print(p1)
}

get_cytoband=function(chr,pos,genome_build="hg19"){
  cyto=read.table(sprintf("liftover_resources/cytoBand.%s.txt.gz",genome_build),stringsAsFactors = FALSE)
  out=data.frame(chr=chr,pos=pos,cyto="",cytotype="",stringsAsFactors = FALSE)
  uchr=unique(chr)
  for(chrom in uchr){
    idx=which(chr==chrom)
    cyto2=cyto[which(cyto$V1==sprintf("chr%s",chrom)),]
    intervals=c(cyto2$V2,cyto2$V3[dim(cyto2)[1]])
    idx2=findInterval(pos[idx],intervals,all.inside = TRUE)
    
    out$cyto[idx]=cyto2$V4[idx2]
    out$cytotype[idx]=cyto2$V5[idx2]
  }
  out
}

get_cytoband_hg19_38=function(details){
  cyto1=get_cytoband(details$Chrom,details$Pos,genome_build="hg19")
  coord_hg38=with(details,liftover2(sprintf("chr%s:%s:%s",Chrom,Pos,Alt),Chrom,Pos,"hg19ToHg38"))
  cyto2=get_cytoband(coord_hg38$chr,coord_hg38$position,genome_build="hg19")
  details$cyto_hg19=cyto1$cyto
  details$cyto_hg38=cyto2$cyto
  details$Chrom_hg38=coord_hg38$chr
  details$Pos_hg38=coord_hg38$position
  details
}

liftover2=function(snp,chr,position,hgflag="hg18ToHg19"){
  out=data.frame(chr=sprintf("chr%s",chr),
                 pos1=sprintf("%d",position),pos2=sprintf("%d",position+1) ,
                 snp=snp,stringsAsFactors=FALSE)
  tmp1=tempfile(tmpdir="/lustre/scratch119/casm/team273jn/nw14/tmp")
  tmp2=tempfile(tmpdir="/lustre/scratch119/casm/team273jn/nw14/tmp")
  write.table(out,tmp1,row.names=FALSE,col.names=FALSE,quote=FALSE)
  cmd=sprintf("liftover_resources/liftOver %s liftover_resources/%s.over.chain.gz %s %s.error",
              tmp1,hgflag,tmp2,tmp2)
  cat(cmd,"\n")
  system(cmd,ignore.stderr=FALSE)
  test=read.table(tmp2,stringsAsFactors=FALSE)
  
  idx=match(snp,test$V4)
  #if(length(which(is.na(idx)>0)))
  #{
  #  stop("SNP mismatch")
  #}
  data.frame(snp=snp,chr=chr,position=test$V2[idx],stringsAsFactors=FALSE)
}
