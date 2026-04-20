require("MutationalPatterns")
require("hdp")
require("dplyr")
require("tibble")
require("stringr")
if(!exists("SPECIES_BUILD")){
  stop("Please create SPECIES (Hsapiens, Mmusculus ) SPECIES_BUILD and set to e.g. hg19,mm10")
}
if(!exists("ref_genome") ){
  ref_genome= sprintf("BSgenome.%s.UCSC.%s",SPECIES,SPECIES_BUILD)
  library(ref_genome, character.only = TRUE)
}


get_pcawg60=function(csvfile="https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.2_SBS_GRCh37.txt"){
  res=read.table(csvfile,head=TRUE,stringsAsFactors = FALSE,sep="\t")  %>% column_to_rownames("Type")
  res[MutationalPatterns:::TRIPLETS_96,]
  
}


##Get mutation count matrix from VCFs
get_mut_matrix_from_vcfs=function(vcf_files,sample_names){
  vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
  mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
}

extract_details_from_hdpres=function(res,priorsigs,mutcount,b.is.denovo=TRUE,threshold=0.9){
  N=dim(mutcount)[2]
  ##check integrity of mutcount
  dps=dp(final_hdpState(chains(res)[[1]]))
  chkcount=tail(sapply(dps, function(x) x@numdata),N)
  idx=which(colSums(mutcount)!=chkcount)
  if(length(idx)>0){
    browser()
    stop("matrix mismatch!")
  }
  comp_distn= comp_categ_distn(res)
  dp_distn=comp_dp_distn(res)
  sigs=comp_distn$mean
  cdc=comp_distn$cred.int
  ##browser()
  sigs_lb95=t(sapply(1:length(cdc),function(x) cdc[[x]][1,]))
  sigs_ub95=t(sapply(1:length(cdc),function(x) cdc[[x]][2,]))
  contr=tail(dp_distn$mean,N)
  cdp=dp_distn$cred.int
  contr_lb95=tail(t(sapply(1:length(cdp),function(x) cdp[[x]][1,])),N)
  contr_ub95=tail(t(sapply(1:length(cdp),function(x) cdp[[x]][2,])),N)
  #Only keep signatures that have non-zero lower bound in at least one sample...
  idx=which(apply(contr_lb95,2,max)>0)
  if(length(idx)==0){
    stop("No significant sigs found!")
  }
  if(length(idx)==1){
    cat("Warning. Only 1 signature found!")
  }
  sigs_lb95=sigs_lb95[idx,,drop=FALSE]
  sigs_ub95=sigs_ub95[idx,,drop=FALSE]
  sigs=sigs[idx,,drop=FALSE]
  contr_lb95=contr_lb95[,idx,drop=FALSE]
  contr_ub95=contr_ub95[,idx,drop=FALSE]
  contr=contr[,idx,drop=FALSE]
  rownames(contr)=colnames(mutcount)
  if(b.is.denovo){
    ##Try ma
    csm=cos_sim_matrix(t(sigs),priorsigs)
    labels=apply(csm,1,function(x){
      idx=which.max(x)
      if(x[idx]>threshold){
        colnames(csm)[idx]
      }else{
        NA
      }})
    idx=which(!is.na(labels))
    if(length(idx)>0){
      rownames(sigs)[idx]=labels[idx]
    }
    idx2=which(is.na(labels))
    if(length(idx2)>0){
      rownames(sigs)[idx2]=sprintf("N%s",1:length(idx2))
    }
    csmmax=apply(csm,1,function(x){max(x)})
    names(csmmax)=rownames(sigs)
  }else{
    csmmax=NULL
  }
  colnames(sigs)=rownames(mutcount)
  colnames(contr)=rownames(sigs)
  rownames(contr_lb95)=rownames(contr)
  rownames(contr_ub95)=rownames(contr)
  colnames(contr_lb95)=colnames(contr)
  colnames(contr_ub95)=colnames(contr)
  cnt=colSums(mutcount)
  sigcount=t(sigs) %*% diag(as.vector(cnt %*% contr),nrow=ncol(contr),ncol=ncol(contr))
  activities=diag(cnt) %*% contr
  rownames(activities)=rownames(contr)
  reconstruction=contr %*% sigs
  sample_csm=sapply(1:dim(reconstruction)[1],function(i) cos_sim(reconstruction[i,],mutcount[,i]))
  list(signatures=t(sigs),
       sigs_lb95=t(sigs_lb95),
       sigs_ub95=t(sigs_ub95),
       contr=contr,
       contr_lb95=contr_lb95,
       contr_ub95=contr_ub95,
       reconstruction=reconstruction,
       csm_max=csmmax,
       sample_csm=sample_csm,
       mutcount=mutcount,
       sigcount=sigcount,
       activities=activities
  )
}

##### HDP functions...

##Denove
runHDP_denovo=function(mut_count,burnin=10000,spacing=500,n_per_chain=250,seed=-1,mtype=96){
  mut_count=t(mut_count)
  n_samp=dim(mut_count)[1]
  n_cat=dim(mut_count)[2]
  hdp_mut <- hdp_init(ppindex=c(0,1,rep(2,n_samp)),
                      cpindex=c(1,2,rep(3,n_samp)),
                      hh=rep(1,mtype),
                      alphaa=rep(1,3),
                      alphab=rep(1,3))
  hdp_mut = hdp_setdata(hdp_mut,
                        dpindex=3:(n_samp+2),
                        mut_count)
  if(seed<0){
    seed_offset=round(as.numeric(Sys.time()))
  }else{
    seed_offset=seed
  }
  ##Now run the sampling chains
  chlist=mclapply(1:4,function(i){
    hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=seed_offset+i*2000)
    hdp_posterior(hdp_activated,
                  burnin=burnin,
                  n=n_per_chain,
                  space=spacing,
                  cpiter=3,
                  seed=seed_offset+i*1e5)
  },mc.cores = 4)
  mut_multi=hdp_multi_chain(chlist)
  hdp_extract_components(mut_multi)
}

get_signatures_color_scheme=function(signames){
  cols=RColorBrewer::brewer.pal(9,"Set1")
  cols[9]="lightgrey"
  scheme=data.frame(
    signames=c("SBS1","SBS5","SBS32","SBS23","SBS2","SBS8","SBS9","SBS40","SBS19","SBS_5_40","SBS_19_23","SBS_5_9_40"),
    col=c(cols,cols[2],cols[9],cols[2]),stringsAsFactors = FALSE
  )
  idx=match(signames,scheme$signames)
  #browser()
  if(length(which(is.na(idx))>0)){
    cat("Can't match all signatures to preferred colour scheme. Using General one instead...\n")
    return(hue_pal(h = c(0, 360) + 15, c = 100, l = 65, h.start = 0,direction = 1)(length(signames)))
  }else{
    return(scheme$col[idx])
  }
}



expanded_plot=function(contribution,title,extracol=NULL,css=NULL,lb=NULL,ub=NULL,cex.label=1,renorm=TRUE,b.add.verticals=TRUE){
  k=dim(contribution)[1]
  N=dim(contribution)[2]
  if(renorm){
    qm=sweep(contribution,2,colSums(contribution),"/")
    if(!is.null(lb)){
      lb=sweep(lb,2,colSums(contribution),"/")
      ub=sweep(ub,2,colSums(contribution),"/")
    }
  }else{
    qm=contribution
  }
  cols=get_signatures_color_scheme(rownames(contribution))
  
  
  par(mar=c(5, 4, 4, 2) + 0.1+c(2,4,0,0))
  plot(NULL,xlim=c(0,N),ylim=c(0,k),xaxt="n",yaxt="n",
       main=title,xlab="",ylab="",cex.main=2)
  for(i in 1:k){
    for(j in 1:N){
      rect(xleft=j-0.5,xright=j+0.5,ybottom=i-1,ytop=i+qm[i,j]-1,col=cols[i],border=NA)
      if(!is.null(lb)){
        rect(col="black",xleft=j-0.02,xright=j+0.02,ybottom = i+lb[i,j]-1,ytop=i+ub[i,j]-1,border=NA)
      }
    }
  }
  axis(2,at=seq(0.5,k),labels=rownames(qm),las=2)
  #idx=which(nchar(rownames(qm))>0)
  #axis(2,at=0.5+idx-1,labels=qm[idx],las=2)
  
  if(!is.null(extracol)){
    axis(1,at=seq(1,N),labels=FALSE,las=2)
    mtext(colnames(qm),side=1,line=1,at=seq(1,N),col=extracol,las=2,cex.axis=cex.label)
  }else{
    #axis(1,at=seq(1,N),labels=colnames(qm),las=2,cex.axis=cex.label)
    idx=which(nchar(colnames(qm))>0)
    axis(1,at=idx,labels=colnames(qm)[idx],las=2,cex.axis=cex.label)
  }
  if(b.add.verticals){
    arrows(y0=par("usr")[3],y1=par("usr")[4],x0=0.5:(N+1),length=0,lwd=0.5)
  }
  arrows(y0=0:k,y1=0:k,x1=par("usr")[1],x0=N+0.5,length=0,lwd=0.5)
}

addSpacedAxisVert=function(at,labels,text.size=1,min.space.factor=0.4,xfac=0.05,x0=NULL){
  #browser()
  xl=par("usr")[1:2]
  yl=par("usr")[3:4]
  par(xpd=NA)
  if(is.null(x0)){
    x0=xl[1]
  }
  xm=x0-(xl[2]-x0)*xfac
  
  n=length(at)
  
  
  #d.to=seq(xl[1]+((xl[2]-xl[1])/(2*(n))),xl[2],(xl[2]-xl[1])/(n))
  tick.at=at
  
  d.to=applyMinimumSpacing(tick.at,min.space.factor*(yl[2]-yl[1])/length(at))
  arrows(y0=tick.at,y1=d.to,x1=xm,x0=xl[1],length=0)
  #TODO!
  
  text(y=d.to,x=xm,labels=labels,srt=0,cex=text.size,pos=2)
}
applyMinimumSpacing=function(test,min.space,test.min=min(test)-min.space,test.max=max(test)+2*min.space){
  prev=0
  n=length(test)
  for(i in 1:n){
    if((test[i]-prev)>min.space){
      prev=test[i]
    }else{
      test[i]=prev+min.space
      prev=test[i]
    }
  }
  if(max(test)>test.max){
    ##TODO deal with overextension on the RHS of plot
  }
  test
}

plot_signature_annotated_tree=function(tree,pdx,sig.col.scheme,maxlen=NULL,b.include.csm.text=TRUE){
  if(is.null(pdx$sigbynode)){
    stop("Need to build sigbynode first: see add_sig_to_pdx")
  }
  ##need df of SIG,col
  control=list(col.scheme=sig.col.scheme)
  if(!is.null(maxlen)){
    control$maxlen=maxlen
    control$b.include.csm.text=b.include.csm.text
  }
  res=add_annotation(pdx,
                     tree,
                     add_signature_decomp,control=control)
  leg=legend("topleft",sig.col.scheme$SIG,col = sig.col.scheme$COL,pch=15,pt.cex=2)$rect
  tree
}
add_signature_decomp=function(pdx,##<< PDX object or list including details matrix
                              tree,##<< enhanced phylo returned from plot_tree
                              node,##<< Tree node - maps to pdx$details$node
                              control,##<< Control parameters.. Requires bfield
                              ...)
{
  ##tree,node,res,ptarget){
  info=get_edge_info(pdx,tree,node)
  i=which(tree$edge[,2]==node)
  
  decomp=pdx$sigbynode[[i]]$contr
  csm=pdx$sigbynode[[i]]$csm
  control=add_defaults(control,defaults=list( maxlen=tree$top*0.2,b.include.csm.text=FALSE),
                       mandatory_fields=c("col.scheme"))
  
  
  y0=info$yb
  y1=info$yt
  
  x=info$x
  mh=y1-y0
  maxlen=control$maxlen
  #minlen=mh;##100
  maxlen=min(maxlen,mh)
  
  
  #if(mh<minlen){
  #  return(NULL)
  #}
  #frac=minlen/mh
  ##dd=0.5*(1-frac) ##0.01
  start=y1-0.01*mh
  start=y1-(mh-maxlen)*0.5
  start=y0+maxlen##ifelse(maxlen<(mh-0.1*maxlen),1.1*maxlen,maxlen)
  
  unit=0.98*mh##frac
  unit=maxlen
  #browser()
  if(is.na(decomp[1])){
    rect(xleft=x-0.25,xright=x+0.25,ytop = start,ybottom=start-unit*1,col="grey",border=NA)
    return(NULL)
  }
  sigs=rownames(decomp)
  col.scheme=control$col.scheme
  cols=col.scheme$COL[match(sigs,col.scheme$SIG)]
  for(i in 1:length(decomp)){
    if(decomp[i]>0){
      rect(xleft=x-0.25,xright=x+0.25,ytop = start,ybottom=start-unit*decomp[i],col=cols[i],border=NA)
      start=start-unit*decomp[i]
    }
  }
  if(control$b.include.csm.text){
    text(x-0.25,y0+maxlen-0.5*unit,sprintf("%3.2f",csm))
  }
  NULL
}

add_sig_to_pdx=function(PD,prior,mutprob=NULL,minmut=50,method="MP"){
  sigbynode=list()
  sigs=colnames(prior)
  empty=matrix(NA,ncol=1,nrow=length(sigs))
  rownames(empty)=sigs
  ##data.frame(SIG=sigs,count=NA)
  if(is.null(mutprob)){
    ##fit using mutational patterns...
    for(i in 1:length(PD$pdx$tree_ml$edge.length)){
      cat("#")
      node=PD$pdx$tree_ml$edge[i,2]
      mutc=get_mut_matrix_for_node(PD$pdx,node)
      if(!is.null(mutc) & sum(mutc)>minmut){
        qf=quickfit(as.vector(mutc),prior,method)
        sigbynode[[i]]=list(node=node,contr=qf$contribution,csm=qf$csm,count=sum(mutc))
      }else{
        sigbynode[[i]]=list(node=node,contr=empty,csm=NA,count=sum(mutc))
      }
    }
    cat("\n")
    PD$pdx$dat$sigbynode=sigbynode
    PD$pdx$dat$siglist=sigs
    return(PD)
  }else{
    stop("Please check this is applicable to current setup!")
    sigs=colnames(mutprob)[-c(1,2)]
    prior=prior[,sigs]
    ##Expect mutprob to have form sample <PATIENT>.<node index> (wildtype =0, others are presented in same order as PD$nodes$node)
    descdf=data.frame(id=sprintf("%s.%s",PD$patient,0:length(PD$nodes$driver3[PD$nodes$status>=0])),desc=c("WT",PD$nodes$driver3[PD$nodes$status>=0]))
    mutprob=mutprob %>% filter(sample %in% descdf$id)
    ##Save in the form 
    tmp=markup_tree(PD$pdx$tree_ml,PD$nodes$node[PD$nodes$status>=0])
    tmp$sample=sprintf("%s.%s",PD$patient,tmp$type)
    sigbynode=list()
    for(i in 1:length(tmp$sample)){
      cat("#")
      idx=which(mutprob$sample==tmp$sample[i])
      node=PD$pdx$tree_ml$edge[i,2]
      mutc=get_mut_matrix_for_node(PD$pdx,node)
      if(!is.null(mutc)){
        contr=t(t(mutc) %*% as.matrix(mutprob[idx,sigs]))
        contr=contr/sum(contr)
        csm=cos_sim(as.vector(mutc),as.vector(prior %*% contr))
        sigbynode[[i]]=list(node=node,contr=contr,csm=csm,count=sum(mutc))
      }else{
        sigbynode[[i]]=list(node=node,contr=empty,csm=NA,count=sum(mutc))
      }
      
    }
    cat("\n")
    PD$pdx$dat$sigbynode=sigbynode
    PD$pdx$dat$siglist=sigs
    return(PD)
  }
  
}

get_cohort_count_matrix=function(PDD,...){
  do.call("cbind",lapply(PDD,get_pd_count_matrix,...))
}

get_pd_count_matrix=function(PD,b.pool.early.branches=FALSE,early.threshold=50,b.pool.short.branches=TRUE){
  cat("Processing",PD$patient,"\n")
  nodes=PD$pdx$tree_ml$edge[,2]
  tree=PD$pdx$tree_ml
  tree$edge.length=tree$el.snv.local.filtered ## Only want actual SNVs not adjusted counts.
  if(b.pool.early.branches){
    nh=nodeHeights(tree)
    early.nodes=tree$edge[nh[,2]<=early.threshold & tree$edge.length>0.001,2]
    if(length(early.nodes)==0){
      EARLY=matrix(rep(0,96),ncol=1)
    }else{
      EARLY=matrix(rowSums(do.call("cbind",lapply(early.nodes,function(x) get_mut_matrix_for_node(PD$pdx,x)))),ncol=1)
    }
    colnames(EARLY)=sprintf("%s.EARLY",PD$pdx$meta$prefix)
    cbind(EARLY,do.call("cbind",lapply(setdiff(nodes,early.nodes),function(x) get_mut_matrix_for_node(PD$pdx,x))))
  }else{
    if(b.pool.short.branches){
      ###. RECALL:
      #' @param tree  An ape::phylo object (ultrametric or with meaningful edge lengths).
      #' @param s     Maximum edge length for a branch to be considered "short".
      #' @param t     Minimum total branch-length sum for a group to be reported.
      #' @param z     Maximum time-span (depth window) within which group members
      #'              must fall. Measured as max(end) - min(start) across edges.
      #' @param e     Minimum start depth for candidate edges (the "late" threshold).
      #' @param f     Maximum allowed depth gap between the group's latest start
      #'              and the MRCA of all descendant tips. Enforces topological
      #'              proximity: a small f requires grouped edges to descend from a
      #'              recent common ancestor rather than being scattered across the tree.
      #'
      late_clusters=cluster_late_short_branches(tree,
                                  40,# Try and merge branches shorter than this
                                  30,# Accept this as total branch sum. Sometimes branches will be merged with themselves.
                                  100,# max width of grouping.
                                  early.threshold, # Early threshold
                                  100 ## distance from MRCA - note this should be greater than 40.
      )
      nh=nodeHeights(tree)
      
      early.nodes=tree$edge[nh[,2]< early.threshold & tree$edge.length>0.001,2]
      
      ####browser()
      if(length(early.nodes)>0){
        EARLY=list(EARLY=data.frame(edge_idx=match(early.nodes, tree$edge[,2]),child=early.nodes))
        if(length(late_clusters)>0){
          names(late_clusters)=sapply(late_clusters,function(cl) paste(as.character(cl$child),collapse="_"))
          GROUPS=c(EARLY,late_clusters)
        }else{
          GROUPS=EARLY
        }
        
      }else if(length(late_clusters)>0){
        names(late_clusters)=sapply(late_clusters,function(cl) paste(as.character(cl$child),collapse="_"))
        GROUPS=late_clusters
      }else{
        return(do.call("cbind",lapply(nodes,function(x) get_mut_matrix_for_node(PD$pdx,x))))
      }
      ## check uniqueness
      exclude.nodes=do.call("c",lapply(GROUPS,function(x) x$child))
      if(length(exclude.nodes) != length(unique(exclude.nodes))){
        stop("Error forming disjoint groups in get_pdd_counts..!")
      }
      GROUP_MATS=lapply(GROUPS,function(grp) 
        matrix(rowSums(do.call("cbind",lapply(grp$child,function(x) get_mut_matrix_for_node(PD$pdx,x)))),ncol=1)
      )
      names(GROUP_MATS)=names(GROUPS)
      the_rest=do.call("cbind",lapply(setdiff(nodes,exclude.nodes),
                                      function(x) get_mut_matrix_for_node(PD$pdx,x))
                       )
      out=do.call("cbind",c(GROUP_MATS,list(the_rest=the_rest)))
      colnames(out)[1:length(GROUP_MATS)]=sprintf("%s.%s",PD$patient,names(GROUP_MATS))
      out
    }else{
      do.call("cbind",lapply(nodes,function(x) get_mut_matrix_for_node(PD$pdx,x)))
    }
  }
}



get_mut_matrix_for_node=function(pdx,node){
  idx=which(pdx$dat$details$node==node & pdx$dat$details$TYPE=="SNV")
  if(length(idx)==0){
    return(NULL)
  }
  ##out=with(pdx$dat$details[idx,] %>% filter(TYPE=="SNV"),get_mut_matrix(Chrom,Pos,Ref,Alt))
  out=with(pdx$dat$details[idx,],get_mut_matrix(Chrom,Pos,Ref,Alt))
  colnames(out)=sprintf("%s.%s",pdx$meta$prefix,node)
  out
}

quickfit=function(sig,priorsigs,method="MP"){
  if(method=="MP"){
    fitfunc=MutationalPatterns:::fit_to_signatures
  }else if(method=="MPS"){
    fitfunc=function(mut_matrix,signatures){MutationalPatterns:::fit_to_signatures_strict(mut_matrix,signatures,max_delta = 0.002)$fit_res}
  }else{
    fitfunc=fixed_weight_em
  }
  decomp=fitfunc(matrix(sig,ncol=1),priorsigs)
  ##browser()
  decomp$contribution[,1]=decomp$contribution[,1]/sum(decomp$contribution[,1])
  decomp$csm=cos_sim(decomp$reconstructed[,1],sig)
  decomp
}
#' Cluster late-occurring short branches on a phylogenetic tree (Developed with the help of claude-ai)
#'
#' Identifies groups of short branches that co-occur late in a tree within a
#' narrow time window and share a recent common ancestor. 
#'
#' Algorithm:
#'   1. Compute node depths (root-to-node distance) for the tree.
#'   2. Select candidate edges: length < s AND start depth > e.
#'   3. Sort candidates by start depth (ascending).
#'   4. Greedy sweep: accumulate consecutive candidates into a group while
#'      the group's time span (max end - min start) stays < z and the MRCA
#'      of all descendant tips is within f of the group's latest start.
#'      When a candidate violates either constraint, finalise the current
#'      group (if total branch length > t) and begin a new one.
#'
#' @param tree  An ape::phylo object (ultrametric or with meaningful edge lengths).
#' @param s     Maximum edge length for a branch to be considered "short".
#' @param t     Minimum total branch-length sum for a group to be reported.
#' @param z     Maximum time-span (depth window) within which group members
#'              must fall. Measured as max(end) - min(start) across edges.
#' @param e     Minimum start depth for candidate edges (the "late" threshold).
#' @param f     Maximum allowed depth gap between the group's latest start
#'              and the MRCA of all descendant tips. Enforces topological
#'              proximity: a small f requires grouped edges to descend from a
#'              recent common ancestor rather than being scattered across the tree.
#'
#' @return A list of data.frames, one per qualifying group. Each data.frame
#'   contains columns: edge_idx, parent, child, start, end, len, desc_tips,
#'   and mrca_depth (the root-to-MRCA distance, constant within a group).
#'   Returns an empty list if no groups are found.
#'
#' @details
#' The greedy sweep is order-dependent: candidates are processed by start
#' depth, and each candidate is assigned to at most one group. This is a
#' heuristic — it may not find all valid groupings, but is fast and
#' deterministic.
#'
#' @examples
#' \dontrun{
#' library(ape)
#' groups <- cluster_late_short_branches(
#'   tree = my_tree,
#'   s = 50,    # branches shorter than 50 mutations
#'   t = 100,   # groups must sum to > 100 mutations
#'   z = 5,     # within a 5-year depth window
#'   e = 40,    # starting after age 40
#'   f = 10     # MRCA within 10 years of latest branch start
#' )
#' }
cluster_late_short_branches <- function(tree, s, t, z, e, f) {
  
  n_tips <- length(tree$tip.label)
  
  # --- 1. Node depths (root-to-node distance) ---
  depths <- node.depth.edgelength(tree)
  
  edges   <- tree$edge
  lengths <- tree$edge.length
  
  start <- depths[edges[, 1]]  # depth of parent node
  end   <- depths[edges[, 2]]  # depth of child node
  
  # --- 2. Filter candidate edges: short AND late-starting ---
  candidates <- which(lengths < s & start > e)
  if (length(candidates) == 0) return(list())
  
  df <- data.frame(
    edge_idx = candidates,
    parent   = edges[candidates, 1],
    child    = edges[candidates, 2],
    start    = start[candidates],
    end      = end[candidates],
    len      = lengths[candidates]
  )
  
  # --- 3. Descendant tips per candidate edge ---
  desc_list <- prop.part(tree)
  
  get_desc_tips <- function(node) {
    if (node <= n_tips) {
      return(node)
    } else {
      return(desc_list[[node - n_tips]])
    }
  }
  
  df$desc_tips <- lapply(df$child, get_desc_tips)
  
  # --- 4. Sort by start depth and greedy-sweep into groups ---
  df <- df[order(df$start), ]
  
  groups <- list()
  
  current_idx       <- integer(0)
  current_start     <- NULL
  current_max_start <- NULL
  current_end       <- NULL
  current_sum       <- 0
  current_tips      <- integer(0)
  
  # Helper: check that the MRCA of a tip set is within f of a reference depth
  # (the group's latest start). Returns the mrca_depth on success, or NULL
  # if the constraint fails. For a single tip, the MRCA is the tip itself,
  # which is always at or below max_start, so the constraint is trivially
  # satisfied; we return the tip's depth.
  mrca_check <- function(tips, max_start) {
    if (length(tips) < 2) {
      return(depths[tips[1]])
    }
    m <- getMRCA(tree, tips)
    if (is.null(m)) return(NULL)
    mrca_depth <- depths[m]
    if ((max_start - mrca_depth) < f) {
      return(mrca_depth)
    }
    return(NULL)
  }
  
  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    
    if (length(current_idx) == 0) {
      current_idx       <- i
      current_start     <- row$start
      current_max_start <- row$start
      current_end       <- row$end
      current_sum       <- row$len
      current_tips      <- row$desc_tips[[1]]
      next
    }
    
    # Proposed group after absorbing this edge
    # (start is non-decreasing due to sort, so min is always current_start,
    #  and max_start is always row$start)
    new_max_start <- row$start
    new_end       <- max(current_end, row$end)
    new_width     <- new_end - current_start
    new_sum       <- current_sum + row$len
    new_tips      <- unique(c(current_tips, row$desc_tips[[1]]))
    
    mrca_depth <- mrca_check(new_tips, new_max_start)
    
    if (new_width < z && !is.null(mrca_depth)) {
      # Absorb into current group
      current_idx       <- c(current_idx, i)
      current_max_start <- new_max_start
      current_end       <- new_end
      current_sum       <- new_sum
      current_tips      <- new_tips
    } else {
      # Finalise current group if it passes the total-length threshold
      mrca_depth_final <- mrca_check(current_tips, current_max_start)
      if (current_sum > t && !is.null(mrca_depth_final)) {
        grp <- df[current_idx, ]
        grp$mrca_depth <- mrca_depth_final
        groups[[length(groups) + 1]] <- grp
      }
      
      # Start new group from this edge
      current_idx       <- i
      current_start     <- row$start
      current_max_start <- row$start
      current_end       <- row$end
      current_sum       <- row$len
      current_tips      <- row$desc_tips[[1]]
    }
  }
  
  # Finalise last group
  mrca_depth_final <- mrca_check(current_tips, current_max_start)
  if (current_sum > t && !is.null(mrca_depth_final)) {
    grp <- df[current_idx, ]
    grp$mrca_depth <- mrca_depth_final
    groups[[length(groups) + 1]] <- grp
  }
  
  return(groups)
}


color_tree=function(tree,s, t, z, e, f){
  nh=nodeHeights(tree)
  short_branches=which(tree$edge.length < s & nh[,2]>e)
  tree$color=rep("black",length(tree$edge.length))
  tree$color[short_branches]="red"
  if(length(short_branches)==0){
    return(tree)
  }else{
    clusters=cluster_late_short_branches(tree, s, t, z, e, f)
  }
  if(is.null(clusters) || length(clusters)==0){
    return(tree)
  }
  colors=RColorBrewer::brewer.pal(9,"Set1")
  k=2
  for(cluster in clusters){
    tree$color[cluster$edge_idx]=colors[k]
    k=k+1
  }
  tree
}

report_variants_per_colony=function(PD){
  colonies=setdiff(PD$pdx$tree_ml$tip.label,"zeros")##PD$pdx$meta$clones_short
  colonies_long=PD$pdx$cfg$LABEL[match(colonies,PD$pdx$cfg$SHORT_LABEL)]
  tree=PD$pdx$tree_ml
  details=PD$pdx$dat$details
  nodes=match(colonies,tree$tip.label)
  do.call("rbind",lapply(1:length(nodes),function(i){
    node=nodes[i]
    df=data.frame(node=get_parents(node,tree$edge))
    df$nshared=sapply(df$node,function(node) length(get_samples_in_clade(node,tree=tree)))
    #vars=details %>% filter(node %in% anc.nodes) %>% dplyr::select(Chrom,Pos,Ref,Alt,TYPE,is_localx_excluded)
    vars=details %>% dplyr::select(Chrom,Pos,Ref,Alt,TYPE,GENE,CCDS,HGVS_PROTEIN,VC,is_localx_excluded,node) %>% inner_join(df,by="node")
    vars$colony=colonies_long[i]
    vars$internal_id=PD$INTERNAL_ID
    vars
  }))
}

get_extended_pcawg=function(){
  ## Retrieve HSPC signature reported by Machado et al. 
  lym=read.table("../data/lympho_sbs_sigs.txt",head=T,stringsAsFactors = FALSE)
  pcawg=get_pcawg60()
  pcawg=cbind(pcawg,SBSblood=lym[,"SBSblood"]) %>% as.matrix()
  pcawg
}