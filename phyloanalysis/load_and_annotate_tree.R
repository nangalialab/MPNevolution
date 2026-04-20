##  This source file contains the code for loading/ annotating and caching trees.
# Includes code to retrieve the trees.
require("dplyr")
require("GenomicRanges")
require("rtreefit")
source("colony_qc.R")
#source("../trees/trees.R")
source("plot_tree_annots.R")
source("define_mutant_nodes.R")
source("background_mutrate_model.R")
source("project.R")
source("utils.R")
source("germline_estimates.R")
#source("vaf_ranking.R")
#source("ABC.R")
##Some Project specific settings

MEAN_AGE_AT_DELIVERY_DAYS=40*7 - 14 #Median gestational age=40 and subrtract 14 for LMP->fertilisation
CACHE="../cache"
MUTCOUNTBIN="../cache/mutcountbin.20210322.RDS"

get_phylo_dat=function(prefix,type,job=NULL,res.dir="../results"){
  if(is.null(job)){
    resdir=res.dir
    cat("looking in:",resdir,"for",prefix,"\n")
    hits=system(sprintf("ls %s/%s_%s*info.Rdat",resdir,prefix,type),intern = TRUE)
    return(hits)
  }else{
    load(job)
    pdx=dat$pdx
    pdx$summary$ml=pdx$summary$edge_ml
    pdx$dat$details$TYPE=with(pdx$dat$details,ifelse(nchar(Ref)+nchar(Alt)==2,"SNV","INDEL"))
    pdx$dat$details$node=pdx$tree_ml$edge[pdx$summary$edge_ml,2]
    #pdx=copy_ofs_to_pdx(pdx,dat$inf)
    pdx$meta=dat$inf$meta
    pdx$cfg=dat$inf$cfg
    ##pdx$tree_ml=root(pdx$tree_ml,"zeros",resolve.root = TRUE)
    list(pdx=pdx,inf=dat$inf)
  }
}
#######
##Get metadata

do_patient_paper_plot=function(PD){
  stop("move me")
  PD$pdx=set_color_by_age(PD$pdx,gtc)
  tree=plot_basic_tree(PD$pdx,"",style = "vertical",cex.label=0,cex.terminal=0.5,legpos=NULL)
  legend("topleft",legend = PD$pdx$age.df$age,col=PD$pdx$age.df$color,pch=19,title = "Age at Sample")
  tree
}

set_color_by_age=function(pdx){
  labels=pdx$cfg$LABEL[match(pdx$tree_ml$tip.label,pdx$cfg$SHORT_LABEL)]
  #browser()
  age.df=data.frame(age=sort(unique(pdx$agedf$age_at_sample_exact[pdx$agedf$age_at_sample_exact>0.1])),stringsAsFactors = FALSE)
  age.df=subset(age.df,!is.na(age.df$age))
  age.df$color=RColorBrewer::brewer.pal(5,"Set1")[1:dim(age.df)[1]]
  pdx$tree_ml$tip.color=age.df$color[match(pdx$agedf$age_at_sample_exact,age.df$age)]
  pdx$age.df=age.df
  pdx
}
get_PDD=function(){
  PDD=lapply(PATIENTS,get_PD)
  names(PDD)=PATIENTS
  PDD
}

get_PD=function(patient,b.force.reload=FALSE,...){
  cfile=sprintf("%s/%s.RDS",CACHE,patient)
  cat("looking for",cfile,"\n")
  if(file.exists(cfile) && b.force.reload==FALSE){
    cat("reading cached data")
    return(readRDS(cfile))
  }else{
    PD=get_pd(patient,...)
    PD=add_adjustment_models(PD)
    cat("caching data")
    PD$nodes=PD$nodes %>% dplyr::select(node,driver,status,driver2,driver3,child_count)
    saveRDS(PD,cfile)
    
    return(PD)
  }
}

get_pd=function(patient,b.raw=FALSE,b.remove.germline=TRUE,res.dir="../results"){
  PD=get_phylo_dat(patient,"mpboot",get_phylo_dat(patient,"mpboot",res.dir=res.dir))
  ##
  PD$patient=patient
  cat("Removing redundant data\n")
  PD=remove_redundant_data(PD,b.remove.germline = b.remove.germline )
  PD$pdx$meta$b.remove.germline=b.remove.germline
  #See define_mutant_nodes.R - this does sample specific processing
  cat("Performing patient specific pre-processing\n")
  PD=add_driver_node_info(PD)
  cat("Adding extra driver annot\n")
  ## Adds additional derived driver annotation
  PD=add_driver2_nodeinfo(PD)
  if(!b.raw){
    PD=add_agedf(PD)
    # Requires cache be installed with ../cache/mutcountbin.20210322.RDS present
    PD=add_excluded_region_flag(PD)
    ##Apply filtering
    PD=apply_manual_local_filter(PD)
    
    tmp=apply_adjustment_germline(PD)
    PD$pdx=tmp$pdx_snv_adj
    PD$pdx$meta=PD$inf$meta
    PD$pdx$meta$b.remove.germline=b.remove.germline
    #PD=add_baitset(PD)
    PD=add_derived_fields(PD)
    PD
  }else{
    PD
  }
}

apply_manual_local_filter=function(PD){
  EXCLUDE_FILE="../data/exclusions.txt"
  if(file.exists(EXCLUDE_FILE)){
    excluded=read.table(EXCLUDE_FILE,head=TRUE,stringsAsFactors = FALSE,sep="\t")
    idx=which(excluded$patient==PD$patient)
    if(length(idx)>0){
      key=with(excluded,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
      okey=with(PD$pdx$dat$details,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
      idx=which(okey %in% key)
      cat("manually apply local filter to:")
      print(PD$pdx$dat$details[idx,1:4])
      PD$pdx$dat$details$is_localx_excluded[idx]=1
    }
  }else{
    cat("apply_manual_local_filter:No manual exclude file provided..\n")
  }
  PD
}


remove_redundant_data=function(PD,b.remove.germline=TRUE){
  ##Get rid of unnecessary data
  PD$inf$snv$filtered=NULL
  PD$inf$indel$filtered=NULL
  for(x in c("vaf_d","mtr_d")){
    for( y in c("snv","indel")){
      PD$inf[[y]][[x]]=NULL
    }
  }
  #Remove unnecesary fields
  PD$pdx$dat$gp=NULL
  PD$pdx$dat$mtr_d=NULL
  PD$pdx$geno=NULL
  ##Remove stray germline variants..
  #if(b.remove.germline){
    root.branches=PD$pdx$tree_ml$edge[get_germline_idx(PD$pdx),2]#PD$pdx$tree_ml$edge[which(PD$pdx$tree_ml$edge[,1]==length(PD$pdx$tree_ml$tip.label)+1),2]
    if(b.remove.germline){
      idx=which(PD$pdx$dat$details$node %in% root.branches & PD$pdx$dat$details$keep_embryonic==0)
      idx.branch=which(PD$pdx$tree_ml$edge[,2] %in% root.branches)
      if(length(idx.branch)==0){
        stop("Can't find germline branches!")
      }
      PD$pdx$tree_ml$edge.length[idx.branch]=0
    }else{
      cat("Removing BGLOD_ORIG=1 variants\n")
      idx=which(PD$pdx$dat$details$node %in% root.branches &  ifelse(PD$pdx$dat$details$TYPE=="SNV",PD$pdx$dat$details$BGLOD_ORIG==1,PD$pdx$dat$details$GLOD<1.35))
      
    }
    if(length(idx)>0){
      cat(PD$patient,":Removing",length(idx)," variants from germline branch..\n")
      saveRDS(PD$pdx$dat$details[idx,],sprintf("../cache/%s_REMOVED.RDS",PD$patient))
      PD$pdx$summary=PD$pdx$summary[-idx,]
      N=dim(PD$pdx$dat$details)[1]
      for(x in names(PD$pdx$dat)){
        if(dim(PD$pdx$dat[[x]])[1]!=N){
          stop("Unexpected dimension pdx$dat")
        }
        PD$pdx$dat[[x]]=PD$pdx$dat[[x]][-idx,]
      }
    }
  #}
  PD
}

add_driver2_nodeinfo=function(PD){
  
  nodes=PD$nodes
  #browser()
  
  tree=PD$pdx$tree_ml
  if(dim(nodes)[1]==0){
    PD$nodes=data.frame(node=integer(), driver=character(),status=numeric(), child_count=integer(),driver2=character(),driver3=character(),driver3b=character())
    return(PD)
  }
  
  df=get_all_tree_drivers(PD$pdx,genes=GENES,cv = CV)
  
  df$label3=gsub(":.*","",df$label2)
  nodes$driver2=sapply(nodes$node,function(node) paste(df$label3[which(df$node==node)],collapse=","))
  driver2b=sapply(nodes$node,function(node) paste(df$label[which(df$node==node)],collapse=","))
  
  ##Let's get the full ancestral history of the node too.
  nodes$driver3=sapply(1:length(nodes$node),function(i){
    node=nodes$node[i]
    pnodes=get_parents(node,tree$edge)
    pnodes=intersect(pnodes,nodes$node)
    idx=match(pnodes,nodes$node)
    paste(nodes$driver2[idx],collapse=":")
  })

  ##Let's get the full ancestral history of the node too.
  nodes$driver3b=sapply(1:length(nodes$node),function(i){
    node=nodes$node[i]
    pnodes=get_parents(node,tree$edge)
    pnodes=intersect(pnodes,nodes$node)
    idx=match(pnodes,nodes$node)
    paste(driver2b[idx],collapse=":")
  })

  PD$nodes=nodes
  PD
}


#' Adds excluded region flag and also calculates genome scaling factors to correct for exclusion zones.
#' Note need to rebuild "../cache/mutcountbins.Rds" (MUTCOUNTBIN) whenever additional patients are added to the cohort.
#' @param PD per patient list structure.
add_excluded_region_flag=function(PD){
  bins=readRDS(MUTCOUNTBIN)
  bins$prob=bins$count/sum(bins$count)
  dfs=list(
    globalx=get_all_excluded(),
    localx=get_all_excluded(list(PD),b_use_cache=FALSE)
  )
  for(x in names(dfs)){
    df=dfs[[x]]
    #browser()
    
    
    
    ##Add in to both inf and pdx
    field=sprintf("is_%s_excluded",x)
    PD$pdx$dat$details[[field]]=0
    if(dim(df)[1]>0){
      for(i in 1:dim(df)[1]){
        idx=with(PD$pdx$dat$details,which(Chrom==df$chr[i] & Pos >= df$start[i] & Pos <= df$end[i]))
        if(length(idx)>0){
          PD$pdx$dat$details[[field]][idx]=1
        }
      }
      ##Update INF as well...
      for(type in c("snv","indel")){
        if(nrow(PD$inf[[type]]$details)>0){
          PD$inf[[type]]$details[[field]]=0
          for(i in 1:dim(df)[1]){
            idx=with(PD$inf[[type]]$details,which(Chrom==df$chr[i] & Pos >= df$start[i] & Pos <= df$end[i]))
            if(length(idx)>0){
              PD$inf[[type]]$details[[field]][idx]=1
            }
          }
        }
      }
    }
  }
  ##Use SNPs to estimate the
  #chrinfo=PD$pdx$dat$details %>% filter() %>% group_by(Chrom) %>% summarise(start=min(Pos),end=max(Pos))
  chrinfo=read.table("../data/global.range",head=T,stringsAsFactors = FALSE)
  globrange=GRanges(chrinfo$Chrom, IRanges(start = chrinfo$start, end=chrinfo$end))
  global.length=sum(as.data.frame(globrange)$width)
  #browser()
  global.overlap=sum(as.data.frame(GenomicRanges::intersect(globrange,
                                                            GRanges(dfs$globalx$chr,IRanges(start=dfs$globalx$start,end=dfs$globalx$end))))$width
  )
  local.overlap=sum(as.data.frame(GenomicRanges::intersect(globrange,
                                                           GRanges(dfs$localx$chr,IRanges(start=dfs$localx$start,end=dfs$localx$end))))$width
  )
  ##Naive (also includes sex chromosomes)
  PD$globalx.correction.naive=global.length/(global.length-global.overlap)
  PD$localx.correction.naive=global.length/(global.length-local.overlap)
  localidx=do.call("c",
                   lapply(1:dim(dfs$localx)[1],function(i) {
                     x=dfs$localx[i,]
                     which(bins$chr==x$chr & bins$start>=(1e5*floor(x$start/1e5)) & bins$end<=(1e5*floor(x$end/1e5)))})
  )
  globalidx=do.call("c",
                    lapply(1:dim(dfs$globalx)[1],function(i) {
                      x=dfs$globalx[i,]
                      which(bins$chr==x$chr & bins$start>=(1e5*floor(x$start/1e5)) & bins$end<=(1e5*floor(x$end/1e5)))})
  )
  PD$localx.correction2=1/(1-sum(bins$prob[unique(localidx)]))
  PD$globalx.correction2=1/(1-sum(bins$prob[unique(globalidx)]))
  PD
}


##Gets excluded regions for a list of samples
get_all_excluded=function(PDD=NULL,fname="loh_cna.exclude.txt",b_use_cache=TRUE){
  if(file.exists(fname) && b_use_cache){
    cat("Reading pre-generated exclude file",fname,"\n")
    return(read.table(fname,sep="\t",head=TRUE,colClasses = c("character","integer","integer"),stringsAsFactors = FALSE))
  }
  if(is.null(PDD)){
    stop(sprintf("need to specify PDD to generate %s",fname))
  }
  allf=do.call("rbind",lapply(PDD,function(x) get_cna_loh(x$inf)))
  allf=rbind(allf,data.frame(chr=c("X","Y"),start=c(0,0),end=c(3e8,3e8)))
  grinds=GRanges(allf$chr, IRanges(start = allf$start, end=allf$end))
  chrinfo=get_chr_info()
  allgrinds=GRanges(chrinfo$chr, IRanges(start = chrinfo$start, end=chrinfo$end))
  df=as.data.frame(GenomicRanges::intersect(GenomicRanges::union(grinds,grinds),allgrinds))#Must be a more elegant way of merging intervals - but this works
  colnames(df)[1]="chr"
  df$chr=as.character(df$chr)
  df=df[,c("chr","start","end")]
  ##Now mark these
  if(b_use_cache){
    write.table(df,file=fname,sep="\t",quote=FALSE,row.names = FALSE)
  }
  df
}

get_cna_loh=function(inf){
  do.call("rbind",lapply(c("LOH","CNA"),function(x){
    df=data.frame(chr=character(),start=integer(),end=integer())
    if(length(inf$meta[[x]])>0){
      df=do.call("rbind",lapply(inf$meta[[x]],
                                function(y) data.frame(chr=as.character(y$chr),start=y$start,end=y$end,stringsAsFactors = FALSE)
      )
      )
    }
    df
  }))
}

add_baitset=function(PD){

  alldat.snv=read_tsv(BAITSET_PILEUP_SNV,col_names = TRUE,guess_max=1e6)
  alldat.indel=read_tsv(BAITSET_PILEUP_INDEL,col_names = TRUE,guess_max=1e6)
  alldat=rbind(alldat.snv,alldat.indel)
  mtr=grep("_MTR$",colnames(alldat))
  wtr=grep("_WTR$",colnames(alldat))
  mtr=alldat[,mtr]
  wtr=alldat[,wtr]
  dep=wtr+mtr
  colnames(mtr)=gsub("_MTR","",colnames(mtr))

  colnames(dep)=colnames(mtr)
  alldat$key=with(alldat,sprintf("%s:%s:%s:%s",CHROM,POS,REF, ALT))
  ##We add in all to see background error rate..
  PD$pdx$dat$details$key=with(PD$pdx$dat$details,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref, Alt))
  idx=match(PD$pdx$dat$details$key,alldat$key)
  PD$pdx$dat$mtr=cbind(PD$pdx$dat$mtr,mtr[idx,])
  PD$pdx$dat$dep=cbind(PD$pdx$dat$dep,dep[idx,])
  rownames(PD$pdx$dat$mtr)=NULL
  rownames(PD$pdx$dat$dep)=NULL
  idxname=c(match(PD$pdx$cfg$SHORT_LABEL[which(PD$pdx$cfg$TYPE %in% c("N","T"))],colnames(PD$pdx$dat$mtr)),grep(PD$patient,colnames(PD$pdx$dat$mtr)))

  #PD$pdx$meta$baitsample=colnames(PD$pdx$dat$mtr)[grep(PD$patient,colnames(PD$pdx$dat$mtr))]
  PD$pdx$meta$baitsample=colnames(PD$pdx$dat$mtr)[idxname]
  ##Add in a VAF too great filter
  #PD$pdx$dat$details$p_cna=with(PD$pdx$dat,pbinom(details$mtr_s-1,details$dep_s,prob=0.5,lower.tail = FALSE))
  PD
}

get_chr_info=function(build=SPECIES_BUILD){
  chr_info=read.table(sprintf("../data/genome.%s.fa.fai",build),head=F,stringsAsFactors = FALSE)
  colnames(chr_info)=c("chr","end","cumulative_start","l1","l2")
  chr_info$start=0
  chr_info=chr_info[,c("chr","start","end")]
  chr_info$chr=gsub("^chr","",chr_info$chr)
  ##Restrict to Autosomes + X,Y
  chr_info=chr_info %>% filter(chr %in% c(sprintf("%d",1:22),c("X","Y")))
  chr_info
}

##Adds the node ID to LOH and CNA metadata
add_node_to_cna=function(pdx){
  dtd=get_all_tree_drivers(pdx)
  for(y in c("CNA","LOH")){
    pdx$meta[[y]]=lapply(pdx$meta[[y]],function(x) {
      idx=which(gsub("^LOH_","",x$LABEL)==dtd$label)
      if(length(idx)==1){
        x$node=dtd$node[idx]
      }else{
        x$node=-1
        cat("Missing CNA/LOH on tree\n")
      }
      x
    })
  }
  pdx
}

get_cna_profile=function(pdx){
  if(length(pdx$meta$CNA)==0){
    return(NULL)
  }
  pdx$meta$CNA=pdx$meta$CNA[sapply(pdx$meta$CNA,function(x) length(x$samples)>0)]
  loh=sapply(pdx$meta$CNA,function(x){zeros=rep(0,length(pdx$df$samples))
  zeros[match(x$samples,pdx$df$samples)]=1;
  c(label=x$LABEL,profile=paste(zeros,collapse=""))
  })
  loh=as.data.frame(t(loh),stringsAsFactors=FALSE)
  loh$cv="cna"
  loh$node=pdx$tree_ml$edge[pdx$summary$edge_ml[match(loh$profile,pdx$summary$profile)],2]
  loh
}

get_loh_profile=function(pdx){
  #browser()
  if(length(pdx$meta$LOH)==0){
    return(NULL)
  }
  pdx$meta$LOH=pdx$meta$LOH[sapply(pdx$meta$LOH,function(x) length(x$samples)>0)]
  if(length(pdx$meta$LOH)==0){
    return(NULL)
  }
  loh=sapply(pdx$meta$LOH,function(x){zeros=rep(0,length(pdx$df$samples))
  zeros[match(x$samples,pdx$df$samples)]=1;
  c(label=x$LABEL,profile=paste(zeros,collapse=""))
  })
  loh=as.data.frame(t(loh),stringsAsFactors=FALSE)
  loh=loh %>% filter(grepl("1",profile))
  if(dim(loh)[1]==0){
    return(NULL)
  }
  ###browser()
  loh=loh[grep("1",loh$profile),]
  loh$label=gsub("^LOH_","",loh$label)
  loh$label=gsub("chr9UPD","9pUPD",loh$label)
  loh$cv="loh"
  loh$node=pdx$tree_ml$edge[pdx$summary$edge_ml[match(loh$profile,pdx$summary$profile)],2]
  loh
}


##
add_derived_fields=function(PD){

  PD$nodes$child_count=sapply(PD$nodes$node,function(node) length(get_samples_in_clade(node,PD$pdx$tree_ml)))
  
  PD$nodes=PD$nodes[order(PD$nodes$child_count,decreasing = TRUE),]
  ##Double check tree_ml edge lengths are correctly adjusted coming in.]
  PD$pdx$tree_ml$edge.length=PD$localx.correction2*PD$pdx$tree_ml$el.snv.local.filtered/PD$pdx$tree_ml$sensitivity.snv.local.filtered

  ##Discard previous derived fields..
  PD$pdx$agedf=PD$pdx$agedf[,c("tip.label","age_at_sample_pcy","age_at_sample_exact")]

  PD$pdx$agedf$tip.label=as.character(PD$pdx$agedf$tip.label)
  PD$pdx$agedf$patient=PD$patient
  idx=which(PD$pdx$agedf$tip.label!="zeros")
  colonies=PD$pdx$agedf$tip.label[idx]
  for(field in c("meandepth","nsnv_call","nindel_call","nsnv_adj")){
    PD$pdx$agedf[[field]]=NA
  }
  ##browser()ou
  PD$pdx$agedf$meandepth[idx]=colMeans(PD$pdx$dat$dep[which(PD$pdx$dat$details$is_localx_excluded==0),colonies])
  PD$pdx$agedf$meanvaf=0
  PD$pdx$agedf$rho.bb=0
  PD$pdx$agedf$vafd=0
  tree=PD$pdx$tree_ml
  
  td=t(sapply(colonies,function(clone){
    idxx=with(PD$pdx$dat$details,which((node %in% get_parents(tree$edge,node=match(clone,tree$tip.label))) & TYPE=="SNV" & is_localx_excluded==0))
    MTR=PD$pdx$dat$mtr[idxx,clone]
    DEP=PD$pdx$dat$dep[idxx,clone]
    vaf=sum(MTR)/sum(DEP)
    c(vaf=vaf,findrho2(MTR,DEP))
  }))
  PD$pdx$agedf$meanvaf[idx]=td[,"vaf"]
  PD$pdx$agedf$rho.bb[idx]=td[,"rho.bb"]


  colSums(PD$pdx$dat$mtr[which(PD$pdx$dat$details$is_localx_excluded==0),colonies])/colSums(PD$pdx$dat$dep[which(PD$pdx$dat$details$is_localx_excluded==0),colonies])


  PD$pdx$agedf$vafd[idx]=PD$pdx$agedf$meandepth[idx]*PD$pdx$agedf$meanvaf[idx]

  PD$pdx$agedf$nsnv_call[idx]=with(PD$inf$snv,colSums(1-ofs[which(details$is_localx_excluded==0 & details$filter_bglod==0),colonies]))*PD$localx.correction2
  PD$pdx$agedf$nindel_call[idx]=with(PD$inf$indel,colSums(1-ofs[which(details$is_localx_excluded==0 & details$filter_bglod==0),colonies]))*PD$localx.correction2

  nh_snvadj=nodeHeights(PD$pdx$tree_ml)
  PD$pdx$agedf$nsnv_adj[idx]=nh_snvadj[match(idx,PD$pdx$tree_ml$edge[,2]),2]

  ##Add in the driver status of each colony
  tree=get_colored_markup_tree2(PD$pdx$tree_ml,PD$nodes[PD$nodes$status>=0,])
  #browser()
  nodes=rbind(data.frame(node=-1,driver="WT",status=1,child_count=-1,driver2="WT",driver3="WT",stringsAsFactors = FALSE),PD$nodes %>% dplyr::select(node,driver,status,driver2,driver3,child_count))
  nodes=nodes[which(nodes$status>=0),]
  #nodes$patient=PD$patient
  ##nodes$type=type
  M=dim(nodes)[1]
  tipnode=tree$rate[match(1:dim(PD$pdx$agedf)[1],tree$edge[,2])]
  tipdesc=nodes$driver[tipnode]
  #tipdesc=ifelse(is.na(tipdesc),"WT",tipdesc)
  PD$pdx$agedf$driver=tipdesc
  PD$pdx$agedf$driver3=nodes$driver3[tipnode]
  ##Driver is based on manually curated PD$node
  drivers=get_all_tree_drivers(PD$pdx)
  ##Here we want any that are JAK2 to have priority..
  for(driverkey in c("JAK2","DNMT3A","PPM1D","9pUPD","TET2")){
    nodes=drivers$node[grep(driverkey,drivers$label)]
    drivcol=do.call("c",lapply(nodes,function(node) get_samples_in_clade(node,PD$pdx$tree_ml)))
    PD$pdx$agedf[[driverkey]]=0
    if(length(drivcol)>0){
      PD$pdx$agedf[[driverkey]][match(drivcol,PD$pdx$agedf$tip.label)]=1
    }
  }
  PD$pdx$agedf$driver2=ifelse(PD$pdx$agedf$JAK2>0,"JAK2",ifelse(PD$pdx$agedf$driver=="WT","WT","Other"))
  ##+Add in sensitivities... if necessary
  df=build_sensitivity_PD(PD)
  if(dim(df)[2]!=18){
    #browser()
    #stop("Unexpected result from build_sensivity_PD")
    cat("dim=",dim(df)[2],"\n")
  }
  PD$pdx$agedf=PD$pdx$agedf %>% left_join(df[,-c(3,4,5)],by=c("patient"="patient","tip.label"="sample_short"))
  PD$pdx$agedf$smraw=get_smetric(PD$pdx$tree_ml)
  PD
}

set_driver2=function(df){
  telo=df
  telo$driver2="Other"
  idx=grep("JAK2",telo$driver)
  telo$driver2[idx]="JAK2"
  idx=grep("WT",telo$driver)
  telo$driver2[idx]="WT"
  cf=getcolor_df(telo)
  telo2=telo %>% inner_join(cf)
  if(dim(telo2)[1]!=dim(telo)[1]){
    stop("Lost info")
  }
  telo2
}

findrho2=function(nmuts,depth){
  if(length(nmuts)==0){
    cat("findrho2:warning.. empty nmuts!\n")
    return(c(mu.bb=1e-6,rho.bb=1e-6,L.bb=-Inf,L.bi=-Inf,p.lr=1,mu.bi=1e-6))
  }
  pseudo=1e-6
  loglik=function(par,nmuts,depth){
    idx=which(!is.na(nmuts/depth))
    out=sum(VGAM::dbetabinom(x =nmuts[idx],size = depth[idx],prob = par[2],rho = par[1],log = T))
    if(is.na(out) || is.infinite(out)){
      browser()
    }
    out
  }
  optres=optim(par=c(0.1,0.1),loglik, gr = NULL,method="L-BFGS-B",lower =c(pseudo,pseudo),upper=c(0.8,0.9),control=list(fnscale=-1),nmuts=nmuts,depth=depth)
  #c(rho=optres$par[1],p=optres$par[2])
  L=optres$value
  LNULL=sum(dbinom(nmuts,depth,prob = sum(nmuts)/sum(depth),log = TRUE))
  p.lr=pchisq(2*(L-LNULL),df=1,lower.tail = FALSE)
  c(mu.bb=optres$par[2],rho.bb=optres$par[1],L.bb=L,L.bi=LNULL,p.lr=p.lr,mu.bi=sum(nmuts)/sum(depth))
}



#' Sets sensitivity fields "sensitivity.local.vafd" in PD$pdx$tree_ml and also the adjusted total burden PD$pdx$agedf$nsnv_adj_vafd.
#' @param PD per patient list structure.
add_vaf_depthmodel=function(PD,model){
  PD$pdx$agedf$sensitivity.local.vafd=1
  PD$pdx$agedf$sensitivity.local.vafd[idx]=predict(fm2,newdata=PD$pdx$agedf[idx,])
  tree=PD$pdx$tree_ml
  tree$sensitivity.local.vafd=sapply(tree$edge[,2],function(node){
  clade=get_samples_in_clade(node,tree)
  with(PD$pdx$agedf,1-prod(1-sensitivity.local.vafd[match(clade,tip.label)]))
})
  tree$edge.length=PD$localx.correction2*tree$el.snv.local.filtered/tree$sensitivity.local.vafd
  nh_snvadj_vafd=nodeHeights(tree)
  PD$pdx$agedf$nsnv_adj_vafd=0
  PD$pdx$agedf$nsnv_adj_vafd[idx]=nh_snvadj_vafd[match(idx,PD$pdx$tree_ml$edge[,2]),2]
  PD$pdx$tree_ml$sensitivity.local.vafd=tree$sensitivity.local.vafd
  PD
}

#' Sets sensitivity fields "sensitivity.local.vafd" in PD$pdx$tree_ml and also the adjusted total burden PD$pdx$agedf$nsnv_adj_vafd.
#' @param PD per patient list structure.
add_adjustment_models=function(PD,predictfn=function(vaf,depth,vardepth){(1-0.0604683)/(1+exp(-(-2.3317652+0.5582742*vaf*depth+0*vardepth/depth+0*depth)))}){
  ##First thing is to check for material non-clonality.
  tree=PD$pdx$tree_ml
  tree$edge.length=PD$localx.correction2*tree$el.snv.local.filtered/tree$sensitivity.snv.local.filtered
  treeo=tree
  tips=which(tree$tip.label!="zeros")
  tipc=tree$tip.label[tips]
  ##Make this rerunnable get rid of per.sample...
  idx=grep("per.sample.sensitivity",colnames(PD$pdx$agedf))
  if(length(idx)>0){
    PD$pdx$agedf=PD$pdx$agedf[,-idx]
  }
  idx=match(tipc,PD$pdx$agedf$tip.label)
  meanvaf=PD$pdx$agedf$meanvaf[idx]
  meandepth=PD$pdx$agedf$meandepth[idx]
  vardepth=PD$pdx$agedf$vardepth[idx]
  adjfactor=function(vaf,depth,vardepth){predictfn(ifelse(vaf>0.5,0.5,vaf),depth,vardepth)/predictfn(0.5,depth,vardepth)}
  sdf=data.frame(tips,
                 tip.label=tipc,
                 meanvaf=meanvaf,
                 meandepth=meandepth,
                 vardepth=vardepth,
                 per.sample.vafd=meanvaf*meandepth,
                 vaf05.sample.sensitivity=tree$sensitivity.snv.local.filtered[match(tips,tree$edge[,2])]
  )

  sdf$per.sample.sensitivity.reg=predictfn(
    ifelse(sdf$meanvaf>0.5,0.5,sdf$meanvaf),
    sdf$meandepth,
    sdf$vardepth)
  sdf$per.sample.sensitivity.hybrid=adjfactor(ifelse(sdf$meanvaf>0.5,0.5,sdf$meanvaf),
                                              sdf$meandepth,
                                              sdf$vardepth)*sdf$vaf05.sample.sensitivity

  ##
  tree$per.branch.sensitivity.reg=sapply(tree$edge[,2],function(node){
    clade=get_samples_in_clade(node,tree)
    1-prod(1-sdf$per.sample.sensitivity.reg[match(clade,sdf$tip.label)])
  })

  tree$per.branch.sensitivity.germline=sapply(tree$edge[,2],function(node){
    clade=get_samples_in_clade(node,tree)
    1-prod(1-sdf$vaf05.sample.sensitivity[match(clade,sdf$tip.label)])
  })
  ##
  tree$per.branch.sensitivity.hybrid=sapply(tree$edge[,2],function(node){
    clade=get_samples_in_clade(node,tree)
    1-prod(1-sdf$per.sample.sensitivity.hybrid[match(clade,sdf$tip.label)])
  })
  tree$per.branch.sensitivity.reg[which(is.na(tree$per.branch.sensitivity.reg))]=1
  tree$per.branch.sensitivity.germline[which(is.na(tree$per.branch.sensitivity.germline))]=1
  tree$per.branch.sensitivity.hybrid[which(is.na(tree$per.branch.sensitivity.hybrid))]=1
  tree$per.branch.sensitivity.germline.multi=tree$sensitivity.snv.local.filtered
  ##This measures the extent to which the parameteric lookup shared across samples differs
  ##from that expected by an independent model (perhaps because sites that are difficult call result in correlated failure)
  h2=tree$sensitivity.snv.local.filtered/tree$per.branch.sensitivity.germline
  #h2 will be used to adjust the multi-sample sensitivity accordingly
  modsens=function(factor,sensitivity){ifelse(is.na(factor*sensitivity) | factor*sensitivity>1,1,factor*sensitivity)}
  tree$per.branch.sensitivity.hybrid.multi=modsens(h2,tree$per.branch.sensitivity.hybrid)
  for(type in c("reg","germline.multi","hybrid","hybrid.multi")){
    field=sprintf("per.branch.sensitivity.%s",type)
    field2=sprintf("nsnv_adj.%s",type)
    tree$edge.length=PD$localx.correction2*tree$el.snv.local.filtered/tree[[field]]
    nh=nodeHeights(tree)
    ##agedf is parallel to tree..
    idx=which(PD$pdx$agedf$tip.label!="zeros")
    PD$pdx$agedf[[field2]]=0
    PD$pdx$agedf[[field2]][idx]=nh[match(idx,PD$pdx$tree_ml$edge[,2]),2]


  }
  PD$pdx$agedf=PD$pdx$agedf %>% left_join(sdf[,c("tip.label",sprintf("per.sample.sensitivity.%s",c("reg","hybrid")))],by="tip.label")


  PD$pdx$tree_ml=tree
  PD=apply_adjustment(PD,type="hybrid.multi")
  PD
}

#' Adjusts edge.length to reflect the specified sensitivity adjusment
apply_adjustment=function(PD,type="hybrid.multi"){
  tree=PD$pdx$tree_ml
  field=sprintf("per.branch.sensitivity.%s",type)
  if(is.null(tree[[field]])){
    stop("Invalid type requested for adjustment")
  }
  PD$pdx$tree_ml$edge.length=PD$localx.correction2*tree$el.snv.local.filtered/tree[[field]]
  PD
}


#' For a specified tree and list of nodes (given in nodeinfo$node) an integer valued vector "rate" that classifies the branches by the most recent
#' ancestor in the node list.
#' @param tree
#'  @param nodes  A data.frame contain a "node" column and a "status" column. Only nodes with status>=0 are are given a rate.
get_colored_markup_tree=function(tree,nodes){
  cols=RColorBrewer::brewer.pal(10,"Paired")
  tmp=markup_tree(tree,nodes)
  tree$color=cols[tmp$type+1]
  tree
}

#' For a specified tree and list of nodes (given in nodeinfo$node) an integer valued vector "rate" that classifies the branches by the most recent
#' ancestor in the node list.
#' @param tree
#'  @param nodinfo  A data.frame contain a "node" column and a "status" column. Only nodes with status>=0 are are given a rate.
get_colored_markup_tree2=function(tree,nodeinfo){
  cols=RColorBrewer::brewer.pal(10,"Paired")
  nodes=nodeinfo$node[which(nodeinfo$status>=0)]
  ##The following returns a list
  tmp=markup_tree(tree,nodes)
  excluded=which(nodeinfo[which(nodeinfo$status>=0),"status"]==0)
  tree$color=ifelse(tmp$type %in% c(excluded),"grey",cols[tmp$type+1])
  tree$rate=tmp$type+1
  tree
}

markup_tree=function(tree,
                     switch_nodes##<< child node ids of edges with a switch in lambda.
){
  ###
  parentidx=match(tree$edge[,1],tree$edge[,2])-1
  parentidx=ifelse(is.na(parentidx),-1,parentidx)
  N=length(tree$tip.label)
  #Index of tips
  tipidx=match(1:N,tree$edge[,2])
  edges=tree$edge
  m=tree$edge.length
  ##Parents of every edge including the edge itself [order from root to tip]
  parentlist=lapply(edges[,2],function(node) {
    ##We need to make sure that parent edges are ordered from root to tip
    parents=rev(get_parents(node,edges))##
    if(length(parents)>0){
      #parents
      match(parents,edges[,2])
    }else{
      parents
    }
  }
  )
  ##Index of applicable lambda rate...
  rates=rep(0,length(m))
  if(length(switch_nodes)>0){
    for(i in 1:length(switch_nodes)){
      snode=switch_nodes[i]
      ##We want to include the parent node as well
      nodes=c(snode,get_all_node_children(node = snode,tree = tree))
      ##Check that nodes in correct order.
      if(length(unique(rates[match(nodes,tree$edge[,2])]))>1){
        stop("switch nodes provided in wrong order!")
      }
      ##browser()
      idxanc=match(nodes,tree$edge[,2])
      if(length(idxanc)>0){
          rates[idxanc]=i
      }
    }
  }
  #Rates of parent node.
  ratesp=rates[match(tree$edge[,1],tree$edge[,2])]
  #If node has no parent assume the wild type rate
  ratesp[which(is.na(ratesp))]=0
  #Index of edges where rate switches occur.
  #if(b_pool_rates){
  nlambda=max(rates)+1
  list(type=rates,ntype=nlambda)
}

build_sensitivity_PD=function(PD,type="snv",mtrmin=5,minvaf=0){
  x=PD
  idx=get_germline_sites(x$inf,type,mtrmin = mtrmin,minvaf=minvaf)
  samples=x$inf$meta$clones_short
  depth=apply(x$inf[[type]]$dep[idx,samples],2,median,na.rm=TRUE)
  meandepth=apply(x$inf[[type]]$dep[idx,samples],2,mean,na.rm=TRUE)
  vardepth=apply(x$inf[[type]]$dep[idx,samples],2,var,na.rm=TRUE)
  sensitivity=apply(1-x$inf[[type]]$ofs[idx,samples],2,mean,na.rm=TRUE)
  tmp=x$inf[[type]]$details[idx,]
  tmp=add_trin(tmp)
  idx2=grepl("\\[C>T\\]G",tmp$trin)
  is_ctcpg=idx2
  sensitivity_ctcpg=apply(1-x$inf[[type]]$ofs[idx[idx2],samples],2,mean,na.rm=TRUE)
  #local
  idx2=which(tmp$is_localx_excluded==0)
  sensitivity_local=apply(1-x$inf[[type]]$ofs[idx[idx2],samples],2,mean,na.rm=TRUE)

  nofs=apply(1-x$inf[[type]]$ofs[idx[idx2],samples],2,sum,na.rm=TRUE)
  ofscnt=length(idx[idx2])





  idx2=which(tmp$is_localx_excluded==0 & is_ctcpg)
  sensitivity_local_ctcpg=apply(1-x$inf[[type]]$ofs[idx[idx2],samples],2,mean,na.rm=TRUE)

  #global
  idx2=which(tmp$is_globalx_excluded==0)
  sensitivity_global=apply(1-x$inf[[type]]$ofs[idx[idx2],samples],2,mean,na.rm=TRUE)
  idx2=which(tmp$is_globalx_excluded==0 & is_ctcpg)
  sensitivity_global_ctcpg=apply(1-x$inf[[type]]$ofs[idx[idx2],samples],2,mean,na.rm=TRUE)


  idxpass=which(x$pdx$dat$details$TYPE==toupper(type) & x$pdx$dat$details$Chrom %in% 1:22)
  
  idxpass=setdiff(idxpass,do.call("c",lapply(c(x$inf$meta$LOH,x$inf$meta$CNA),
                                             function(y) with(x$pdx$dat$details,which(TYPE==toupper(type) & Chrom==y$chr & Pos>=y$start & Pos<y$end))
  )
  ))


  npass=apply(1-x$pdx$dat$ofs[idxpass,samples],2,sum,na.rm=TRUE)
  vaf=apply(x$pdx$dat$geno[idxpass,samples]==1,2,mean,na.rm=TRUE)


  idx=get_non_germline_sites(x$inf,type)
  ##Number that are called that are not germline...
  ncall=apply(1-x$inf[[type]]$ofs[idx,samples],2,sum,na.rm=TRUE)
  tmp=x$inf[[type]]$details[idx,]
  tmp=add_trin(tmp)
  idx2=grepl("\\[C>T\\]G",tmp$trin)
  is_ctcpg=idx2
  ncall_ctcpg=apply(1-x$inf[[type]]$ofs[idx[idx2],samples],2,sum,na.rm=TRUE)
  ##Now filter according to local and global
  #local
  idx2=which(x$inf[[type]]$details$is_localx_excluded[idx]==0)
  ncall_local=apply(1-x$inf[[type]]$ofs[idx[idx2],samples],2,sum,na.rm=TRUE)

  idx2=which(x$inf[[type]]$details$is_localx_excluded[idx]==0 & is_ctcpg)
  ncall_local_ctcpg=apply(1-x$inf[[type]]$ofs[idx[idx2],samples],2,sum,na.rm=TRUE)

  ##glabal
  idx2=which(x$inf[[type]]$details$is_globalx_excluded[idx]==0)
  ncall_global=apply(1-x$inf[[type]]$ofs[idx[idx2],samples],2,sum,na.rm=TRUE)
  idx2=which(x$inf[[type]]$details$is_globalx_excluded[idx]==0 & is_ctcpg)
  ncall_global_ctcpg=apply(1-x$inf[[type]]$ofs[idx[idx2],samples],2,sum,na.rm=TRUE)
  data.frame(patient=x$inf$meta$prefix,
             sample_short=samples,
             sample=x$inf$meta$clones,
             lcm=x$inf$cfg$LCM[match(samples,x$inf$cfg$SHORT_LABEL)],
             depth=depth,
             #meandepth=meandepth,
             vardepth=vardepth,
             sensitivity=sensitivity,
             sensitivity_ctcpg=sensitivity_ctcpg,
             sensitivity_local=sensitivity_local,
             sensitivity_global=sensitivity_global,
             sensitivity_local_ctcpg=sensitivity_local_ctcpg,
             sensitivity_global_ctcpg=sensitivity_global_ctcpg,
             npass=npass,
             ncall=ncall,
             ncall_ctcpg,
             ncall_local=ncall_local,
             ncall_global=ncall_global,
             ncall_local_ctcpg=ncall_local_ctcpg,
             ncall_global_ctcpg=ncall_global_ctcpg,
             nofs=nofs,
             ofscnt=rep(ofscnt,length(nofs)),
             meandepth2=meandepth,
             stringsAsFactors=FALSE)
}

get_germline_sites=function(inf,type="snv",mtrmin=5,minvaf=0){
  details=inf[[type]]$details
  filters=setdiff(colnames(details)[grep("^filter",colnames(details))],c("filter_bglod","filter_count","filter_basic_germline","filter_beta_binomial"))
  normal=inf$cfg$SHORT_LABEL[match(inf$meta$normal,inf$cfg$LABEL)]
#  browser()
  if(!is.na(normal)){
    nvaf=inf[[type]]$mtr[,normal]/inf[[type]]$dep[,normal]
    nmtr=inf[[type]]$mtr[,normal]
    min1000g=0
  }else{
    nvaf=rep(minvaf+1,dim(inf[[type]]$mtr)[1])
    nmtr=rep(mtrmin+1,dim(inf[[type]]$mtr)[1])
    min1000g=1
  }

  if(grepl("hg",SPECIES_BUILD)){
    idx=which(rowSums(details[,filters])==0 & details$filter_bglod==1 & nmtr>mtrmin & details$Chrom %in% (get_chr_info() %>% filter(!grepl("X|Y",chr)) %>% pull(chr)) & nvaf>minvaf & details$is_localx_excluded==0 & details$G1000_AC >= min1000g)##Need at least 5 mtr reads in normal
  }else{
    idx=which(rowSums(details[,filters])==0 & details$filter_bglod==1 & nmtr>mtrmin & details$Chrom %in% (get_chr_info() %>% filter(!grepl("X|Y",chr)) %>% pull(chr)) & nvaf>minvaf & details$is_localx_excluded==0 )##Need at least 5 mtr reads in normal
    
  }
  cat("Found",inf$meta$prefix,": Found",length(idx),"germline loci")
  idx
}

get_non_germline_sites=function(inf,type="snv"){
  details=inf[[type]]$details
  filters=intersect(c("filter_bglod","filter_count","filter_basic_germline"),colnames(details))
  normal=inf$cfg$SHORT_LABEL[match(inf$meta$normal,inf$cfg$LABEL)]
  idx=which(rowSums(details[,filters])==0  & details$Chrom %in% 1:22 & details$is_localx_excluded==0)##Need at least 5 mtr reads in normal
  cat("Found",inf$meta$prefix,": Found",length(idx),"autosomal non-germline loci")
  idx
}

get_recapinfo=function(){
  ##test=read.table("../data/Recapture-colonies.txt",head=T,stringsAsFactors = FALSE,sep="\t")
  test=read.table("../data/Recapture-colonies-updated-28sep20.txt",head=T,stringsAsFactors = FALSE,sep="\t")

  patients=read.table("../data/patient_info.txt",head=TRUE,stringsAsFactors = FALSE)
  test$PDID=str_extract(test$Supplier.Name,"PD[0-9]+")
  test=test %>% left_join(patients)
  test$SAMPLE_DATE=as.Date(test$Date.of.sample,"%d/%m/%Y")
  test$Age.of.sample=as.integer(as.Date(test$Date.of.sample,"%d/%m/%Y")-as.Date(test$DOB,"%d/%m/%Y"))/365.25
  test
}

get_patient_summary=function(subject){
  colonyinfo=get_sample_age_table() %>% group_by(patient,age_at_sample_exact,age_at_sample,DOB,DATE_OF_DIAGNOSIS,SAMPLE_DATE) %>% summarise(N=n()) %>% filter(patient==subject)
  recapinfo=get_recapinfo() %>% filter(PDID==subject)
  ##browser()
  c1=colonyinfo[,c("patient", "age_at_sample_exact","SAMPLE_DATE")]
  c1$cell_type="BFU-E-Colony"
  c1$BaitLabel=NA
  c1$phase="Colony"
  c1$ID=sprintf("COLONY%s",round(c1$age_at_sample_exact))
  if(dim(recapinfo)[1]>0){
    r1=recapinfo[,c("PDID","Age.of.sample","SAMPLE_DATE","Cell.type","BaitLabel")]
    r1$phase="Recapture"
    r1$ID=recapinfo$Supplier.Name
    colnames(r1)=colnames(c1)
  }else{
    r1=NULL
  }
  ##browser()
  out=rbind(as.data.frame(c1),as.data.frame(r1))
  out[order(out$age_at_sample_exact),c("patient","ID","age_at_sample_exact","cell_type","phase","SAMPLE_DATE","BaitLabel")]
}

expand_short_branches=function(pdx,prop=0.05){
  tree=plot_tree(pdx$tree_ml,b_do_not_plot = TRUE)
  minmuts=round(prop*max(tree$coords$a1))
  pdx$tree_ml$color=ifelse(pdx$tree_ml$edge.length<minmuts,"grey","black")
  pdx$tree_ml$edge.length.unexpanded=pdx$tree_ml$edge.length
  pdx$tree_ml$edge.length=ifelse(pdx$tree_ml$edge.length<minmuts,minmuts,pdx$tree_ml$edge.length)
  pdx$minmuts=minmuts
  pdx
}

source("project.R")
