
get_germline_idx=function(pdx,b_fall_over_if_error=TRUE){
  germline_profile=paste(ifelse(pdx$df$samples=="zeros",0,1),collapse="")
  idx=which(pdx$tree_ml$label==germline_profile)
  if(length(idx)!=1){
    if(b_fall_over_if_error){
      stop("Can't find germline profile!")
    }else{
      return(-1)
    }
  }
  idx
}

set_germline_branch_to_zero=function(pdx){
  idx=get_germline_idx(pdx)
  cat("WARNING: Setting germline of length",length(idx),"to zero!")
  pdx$tree_ml$edge.length[idx]=0
  pdx
}


convert_nodes_to_el=function(tree,nodelist,idx.germline){
  el=rep(0,dim(tree$edge)[1])
  counts=as.data.frame(table(nodelist))
  #browser()
  el[match(counts$nodelist,tree$edge[,2])]=counts$Freq
  el[idx.germline]=0
  el
}




#' Creates germline adjustment factors..
#' Note need to rebuild "../cache/mutcountbins.Rds" (MUTCOUNTBIN) whenever additional patients are added to the cohort.
#' @param PD per patient list structure.
apply_adjustment_germline=function(PD){
  #browser()
  pdx=PD$pdx
  inf=PD$inf
  type="snv"
  colonies=inf$cfg %>% filter(TYPE=="C") %>% pull(SHORT_LABEL)
  ##browser()
  inf[[type]]$details$totalvaf=rowSums(inf[[type]]$mtr[,colonies])/rowSums(inf[[type]]$dep[,colonies])
  #pdx$summary$type=pdx$dat$details$TYPE
  if(PD$pdx$meta$b.remove.germline){
    idx.germline=get_germline_idx(pdx)
  }else{
    idx.germline=c()
  }
  ##get raw SNV and filtered SNV
  counts.snv=with(pdx$dat$details,as.data.frame(table(node[which(TYPE=="SNV")])))
  el.snv=convert_nodes_to_el(pdx$tree_ml,
                             with(pdx$dat$details,node[which(TYPE=="SNV")]),
                             idx.germline
  )
  el.snv.global.filtered=convert_nodes_to_el(pdx$tree_ml,
                                             with(pdx$dat$details,node[which(TYPE=="SNV" & is_globalx_excluded==0 )]),
                                             idx.germline
  )
  el.snv.local.filtered=convert_nodes_to_el(pdx$tree_ml,
                                            with(pdx$dat$details,node[which(TYPE=="SNV" & is_localx_excluded==0 )]),
                                            idx.germline
  )
  
  sensitivity.snv=get_edge_sensitivities(pdx,inf,NULL)
  sensitivity.snv.global.filtered=get_edge_sensitivities(pdx,inf,"is_globalx_excluded")
  sensitivity.snv.local.filtered =get_edge_sensitivities(pdx,inf,"is_localx_excluded")
  
  pdxs=pdx
  pdxs$tree_ml$edge.length=el.snv
  pdxsa=pdxs
  pdxsa$tree_ml$el.snv=el.snv
  pdxsa$tree_ml$el.snv.global.filtered=el.snv.global.filtered
  pdxsa$tree_ml$el.snv.local.filtered=el.snv.local.filtered
  
  
  pdxsa$tree_ml$sensitivity.snv=sensitivity.snv
  pdxsa$tree_ml$sensitivity.snv.global.filtered=sensitivity.snv.global.filtered
  pdxsa$tree_ml$sensitivity.snv.local.filtered=sensitivity.snv.local.filtered
  
  pdxsa$tree_ml$edge.length=el.snv/sensitivity.snv
  pdxsa$tree_ml$elo=pdx$tree_ml$edge.length
  
  list(pdx_orig=pdx,pdx_snv=pdxs,pdx_snv_adj=pdxsa)
}

get_edge_sensitivities=function(pdx,inf,exclude_field){
  sapply(1:length(pdx$tree_ml$edge[,2]),function(i){
    colonies=setdiff(get_samples_in_clade(pdx$tree_ml$edge[i,2],pdx$tree_ml),"zeros")
    if(length(colonies)==0){
      ##zeros
      return(1)
    }
    get_germline_sensitivity_multi_sample(inf,colonies,exclude_field=exclude_field)
  })
}

adjust_edge_length=function(pdx,inf,exclude_field){
  pdx$cfg$depth=apply(pdx$dat$dep[,pdx$cfg$SHORT_LABEL],2,mean,na.rm=T)
  el=sapply(1:length(pdx$tree_ml$edge[,2]),function(i){
    #browser()
    colonies=setdiff(get_samples_in_clade(pdx$tree_ml$edge[i,2],pdx$tree_ml),"zeros")
    if(length(colonies)==0){
      return(pdx$tree_ml$edge.length[i])
    }
    sensitivity=get_germline_sensitivity_multi_sample(inf,colonies,exclude_field = exclude_field)
    #cat(sensitivity,":",colonies,"\n")
    #adj=1/sensitivity
    ##Use germline....
    
    
    if(is.na(sensitivity)){
      browser()
    }
    if(sensitivity<0.2){
      cat("Capping edge length adjustment !\n")
      sensitivity=0.2
    }
    
    pdx$tree_ml$edge.length[i]/sensitivity
    
  })
  pdx$tree_ml$elo=pdx$tree_ml$edge.length
  ##pdx$tree_ml$
  pdx$tree_ml$edge.length=el
  pdx=set_germline_branch_to_zero(pdx)
  pdx
}

get_germline_sensitivity_multi_sample=function(inf,samples,type="snv",exclude_field){
  ###browser()
  details=inf[[type]]$details
  filters=setdiff(colnames(details)[grep("^filter",colnames(details))],c("filter_bglod","filter_count","filter_basic_germline","filter_beta_binomial" ))
  normal=inf$cfg$SHORT_LABEL[match(inf$meta$normal,inf$cfg$LABEL)]
  vaf=inf[[type]]$details$totalvaf
  ##Don't want homozygous variants
  if(is.na(inf$meta$normal)){
    idx=which(rowSums(details[,filters])==0 & details$filter_bglod==1 & vaf<0.6)##
  }else{
    if(!is.null(exclude_field)){
      idx=which(rowSums(details[,filters])==0 & details$filter_bglod==1 & inf[[type]]$mtr[,normal]>5 & details[[exclude_field]]==0 & vaf<0.6)##Need at least 5 mtr reads in normal
    }else{
      idx=which(rowSums(details[,filters])==0 & details$filter_bglod==1 & inf[[type]]$mtr[,normal]>5 & vaf<0.6)##Need at least 5 mtr reads in normal
    }
  }
  if(length(idx)<100){
    #browser()
	  warning("Not enough SNPs to estimated sensitivity.. Setting to 1")
	  return(1)
  }
  ##cat("Estimating sensitivity over",length(idx),"sites\n")
  if(length(samples)>1){
    mean(ifelse(rowSums(inf[[type]]$ofs[idx,samples]==0)>0,1,0))
  }else{
    ##browser()
    mean(inf[[type]]$ofs[idx,samples]==0)
  }
}

