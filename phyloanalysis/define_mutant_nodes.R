##This reflects the manual curation of driver nodes.
require("digest")
require("treemut")
get_hash=function(PD){
  digest(list(PD$patient,PD$pdx$tree_ml$edge,PD$pdx$tree_ml$tip.label))
}
check_hash=function(PD,hash){
  if(get_hash(PD)!=hash){
    stop(sprintf("HASH Mismatch for %s",PD$patient))
  }
}

##Special function for 4781
add_duplicate=function(pdx,idx,edges){
  N=dim(pdx$summary)[1]
  #browser()
  idxedge=pdx$summary$edge_ml[idx]
  pdx$summary$edge_ml[idx]=edges[1]
  pdx$summary=rbind(pdx$summary,pdx$summary[idx,])
  pdx$summary$edge_ml[N+1]=edges[2]
  pdx$geno=rbind(pdx$geno,pdx$geno[idx,])
  for(field in names(pdx$dat)){
    pdx$dat[[field]]=rbind(pdx$dat[[field]],pdx$dat[[field]][idx,])
  }
  ##Adjust branch length
  pdx$tree_ml$edge.length[idxedge]=ifelse(pdx$tree_ml$edge.length[idxedge]>1,pdx$tree_ml$edge.length[idxedge]-1,0)
  pdx$tree_ml$edge.length[edges]=pdx$tree_ml$edge.length[edges]+1
  ##reset node..
  pdx$dat$details$node=pdx$tree_ml$edge[pdx$summary$edge_ml,2]
  ##Don't bother with df
  pdx
}

add_driver_node_info=function(PD){
  dff=get_all_tree_drivers(PD$pdx,genes=GENES,cv = CV)
  ##browser()
  #browser()
  if(dim(dff)[1]>0){
    dff$status=-1
    nodes=do.call("rbind",lapply(dff$node,function(node) {
      idx=which(dff$node==node)
      driver=paste(dff$label[idx],collapse=",")
      data.frame(node=node,driver=driver,status=1,stringsAsFactors = FALSE)
    }))
    if(any(is.na(nodes$node))){
      print(dff)
      browser()
      stop("fix up annotations in yaml they are not consistent with tree...")
    }
    nodes$child_count=sapply(nodes$node,function(node) length(get_samples_in_clade(node=node,tree=PD$pdx$tree_ml)))
  }else{
    nodes=data.frame(node=integer(),driver=character(),status=integer(),child_count=integer(),stringsAsFactors = FALSE)
  }
 
  PD$nodes=nodes
  PD
}


add_driver_node_info_BESPOKE_EXAMPLE=function(PD){
  patient=PD$patient
  
  if(patient=="PD5147"){
    HASH="a97d897b58f84a560685d763681ecf73"
    check_hash(PD,HASH)
    PD$nodes=data.frame(
      node=c(70),
      #driver=c("TET2:p.S657fs*42+PPM1D:p.T483fs*3"),
      driver=c("PPM1D,TET2"),
      status=c(1),
      stringsAsFactors = FALSE
    )
    return(PD)
  }
}


