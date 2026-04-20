#Invoke this to build the background rate model based on 
init_bins=function(increment=1e5){
  inf=get_chr_info()
  inf=inf[which(inf$chr %in% 1:22),]
  bins=do.call("rbind",
               lapply(1:dim(inf)[1],function(i) {
                 res=data.frame(chr=inf$chr[i],start=seq(0,inf$end[i],by=increment))
                 res$end=res$start+increment
                 res}
               )
  )
  bins$count=0
  bins
}

get_binned_mutcount=function(PD,bins,b.plot=FALSE){
  ##The following identifies which mutations do not belong to CNA/LOH samples and then will count each mutation just once.
  cndat=rbind(get_loh_profile(PD$pdx),get_cna_profile(PD$pdx))
  if(!is.null(cndat)){
    topnodes=unique(cndat$node)
  }else{
    topnodes=c()
  }
  excludenodes=unique(c(unlist(lapply(topnodes,get_all_node_children,tree=PD$pdx$tree_ml)),topnodes))
  cat("ignoring",excludenodes,"\n")
  tree=PD$pdx$tree_ml
  nodes=setdiff(tree$edge[which(tree$edge[,1]!=(length(tree$tip.label)+1)),2],excludenodes)
  if(b.plot){
    #Useful for checking the correct branches have been selected..
    tree=plot_tree(tree,b_do_not_plot = TRUE)
    tree$coords$color="black"
    tree$coords$color[which(tree$edge[,2] %in% excludenodes)]="red"
    tree$coords$color[which(tree$edge[,2] %in% nodes)]="green"
    plot_tree(tree,lwd=2,cex.label = 0)
    title(PD$patient)
  }
  ncolony=length(nodes[which(nodes<=length(PD$pdx$tree_ml$tip.label))])
  #scount=scount+ncolony
  cat(PD$patient,"ncolony=",ncolony,"\n")
  if(length(nodes)>0){
    details=PD$pdx$dat$details[which(PD$pdx$dat$details$node %in% nodes),]
    ##Add to counts..
    for(i in 1:22){
      idx2=which(bins$chr==i)
      if(length(idx2)>0){
        bins$count[idx2]=bins$count[idx2]+sapply(idx2,function(j) length(which(details$Chrom==i & details$Pos>=bins$start[j] & details$Pos<bins$end[j])))
      }
    }
  }
  bins
}

get_all_mutcount_bins=function(PDD){
  bins=init_bins()
  for(x in PDD){
    cat("Processing:",x$patient,"\n")
    bins=get_binned_mutcount(x,bins);
    cat("found",sum(bins$count),"muts\n")
  }
  bins$prob=bins$count/sum(bins$count)
  bins
}

plot_bins=function(bins){
  require(karyoploteR)
  gr=GRanges(sprintf("chr%s",bins$chr), IRanges(start=bins$start,end=bins$end))
  mcols(gr)=data.frame(y=bins$prob)
  ymax=max(bins$prob)
  ymin=0
  kp <- plotKaryotype(genome=SPECIES_BUILD,chromosomes=sprintf("chr%s",1:22),plot.type=1)
  kpAbline(kp, h=0, ymin=ymin, ymax=ymax, lty=2, col="#666666")
  kpAxis(kp, ymin = ymin, ymax=ymax)
  kpLines(kp, data=gr, ymin=ymin, ymax=ymax)
}


if(FALSE){
#source("load_and_annotate_tree.R")
PDD=lapply(PATIENTS,get_pd,b.raw=TRUE)
##Build the model in 100kb bins
bins=init_bins(increment=1e5)
for(x in PDD){
  cat("Processing:",x$patient,"\n")
  bins=get_binned_mutcount(x,bins,b.plot = TRUE);
  cat("found",sum(bins$count),"muts\n")
}
bins$prob=bins$count/sum(bins$count)
PHENO="CML"
saveRDS(bins,sprintf("../cache/mutcountbin.%s.RDS",SPECIES_BUILD))#format(Sys.time(), "%Y%m%d")))

}

