require("RColorBrewer")
require("scales")
require("tidyverse")

set_cedge=function(parent,tree){
  for(child in get_node_children(parent,tree)){
    #cat("\nsetting parent - child ",parent,child,"\n")
    child.idx=which(tree$edge[,2]==child)
    parent.idx=which(tree$edge[,2]==parent)
    if(length(parent.idx)==0){
      pedge=0
    }else{
      pedge=tree$cedge[parent.idx]
    }
    #cat("before:",length(which(tree$cedge>0)),"\n")
    tree$cedge[child.idx]=pedge+tree$edge.length[child.idx];
    #cat("after:",length(which(tree$cedge>0)),"\n")
    #cat(parent,child,":",tree$cedge)
    tree=set_cedge(child,tree)
  }
  tree
}

##Not very efficient -- recursively calculates height as average of children's height.
get_height=function(tree,node){
  #if(is.null(tree$height)){
  #  tree$height=rep(NA,1:length(tree$edge.length))
  #}
  N=length(tree$tip.label)
  if(node<=N){
    return(node)#tree$height[which(tree$edge[,2]==node)]=node
  }else{
    children=get_node_children(node,tree)
    return(mean(sapply(children,function(child) get_height(tree,child))))
  }
}

set_height=function(tree){
  tree$height_end=sapply(tree$edge[,2],function(i) get_height(tree,i))
  tree$height_start=tree$height[match(tree$edge[,1],tree$edge[,2])]
  N=length(tree$tip.label)
  root=N+1
  idx=which(tree$edge[,1]==root)
  tree$height_start[idx]=mean(tree$height_end[idx])
  tree
}

##horizontal elbow
elbow=function(x0,x1,y0,y1,...){
  arrows(x0=x0,y0=y0,y1=y1,length=0,...)
  arrows(x0=x0,x1=x1,y0=y1,length=0,...)
}
##vertical elbow
elbowv=function(x0,x1,y0,y1,...){
  #browser()
  arrows(x0=x0,x1=x1,y0=y0,length=0,...)
  arrows(x0=x1,y0=y0,y1=y1,length=0,...)

}

set_tree_coords=function(atree){
  ##get the cumulative distance from the root.
  tree=atree
  tree$cedge=rep(0,length(tree$edge.length))
  N=length(tree$tip.label)
  root=N+1
  tree=set_cedge(root,tree)
  tt=set_height(tree)
  atree$coords=data.frame(a0=tt$cedge-tt$edge.length,a1=tt$cedge,
                          b0=tt$height_start,b1=tt$height_end,stringsAsFactors = FALSE)
  if(!is.null(atree$color)){
    atree$coords$color=atree$color
  }
  atree
}
plot_tree=function(tree,direction="down",cex.label=1,offset=0,b_do_not_plot=FALSE,lwd=1,bars=NULL,default_edge_color="darkgrey",ymax=NULL,cex.terminal.dots=0,mar=NULL,b.add.scale=TRUE,left.margin.prop=0.1,vspace.reserve=0){
  if(is.null(mar)){
    par(mar=c(1, 1, 1, 3) + 0.1)
  }else{
    par(mar=mar)
  }
  #browser()
  if(!(direction %in% c("down","across"))){
    stop("Unsupported direction provided")
  }
  N=length(tree$tip.label)
  if(is.null(tree$coords)){
    tree=set_tree_coords(tree)
  }
  coords=tree$coords
  
  if(direction=="across"){
    xmax=max(coords$a1)*1.05
    ymax=max(coords$b1)+1
    offset=offset*xmax
  }else{
    if(is.null(ymax)){
      ytop=max(coords$a1)
      ymax=max(coords$a1)*1.05
    }else{  #add 24jun
      #coords[,1:2]=coords[,1:2]+(ymax-max(coords$a1))
      #tree$coords=coords
      ytop=ymax
      ymax=1.05*ymax
    }  #end add
    xmax=max(coords$b1)+1
    offset=offset*ymax
  }
  if(b_do_not_plot){
    return(tree)
  }
  
  if(is.null(bars)){
    YMIN=0-ymax*0.05-vspace.reserve*ymax
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*left.margin.prop),xmax),ylim=c(YMIN,ymax),xlab="",ylab="")
    #plot(NULL,axes=FALSE,xlim=c(0-(xmax*left.margin.prop),xmax),ylim=c(0-(ymax*0.05),ymax),xlab="",ylab="")
  }else{
    ##ymax not supported..
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*left.margin.prop),xmax),ylim=c(0-(ymax*0.15),ymax),xlab="",ylab="")
  }
  idx.tip=match(1:N,tree$edge[,2])
  if(direction=="across"){
    apply(coords,1,function(x) elbow(x[1],x[2],x[3],x[4]))
    text(tree$tip.label,x =coords$a1[idx.tip]+offset ,y=coords$b1[idx.tip],cex = cex.label,pos = 4)
  }else{
    top=ytop
    #top=max(coords$a1)
    m=dim(coords)[1]
    if(is.null(coords$color)){
      col=rep(default_edge_color,m)
    }else{
      col=coords$color
    }
    sapply(1:m,function(i) {x=as.numeric(coords[i,1:4]);elbowv(x[3],x[4],top-x[1],top-x[2],col=col[i],lwd=lwd)})
    if(is.null(tree$tip.color)){
      tipcol="black"
    }else{
      tipcol=tree$tip.color
    }
    if(cex.label>0){
      text(tree$tip.label,y =top-(coords$a1[idx.tip]+offset) ,x=coords$b1[idx.tip],cex = cex.label,pos = 1,col=tipcol)
    }
    if(cex.terminal.dots>0){
      points(y =top-(coords$a1[idx.tip]) ,x=coords$b1[idx.tip],col=tipcol,cex=cex.terminal.dots,pch=19)
    }
  }
  tree$direction=direction
  tree$top=top
  scales=c(0,1,10,50,100,200,500,1000,2000,5000)
  #scale=scales[max(which(ymax/2>scales))]
  scale=scales[max(which(ymax/5>=scales))]
  #browser()
  #cat("scale=",scale,"\n")
  if(b.add.scale){
    axis(side = 4,at=seq(top,-scale,-scale),label=seq(0,top+scale,scale),las=2)
  }
  #arrows(x0=length(tree$tip.label)+0.5,y0=0,y1=scale,length=0.1,code=3,angle=90)
  #text(sprintf("%s Muts",scale),x=length(tree$tip.label)-0.5,y=0.5*scale,pos=4,cex=cex.label,offset=0.1)
  if(!is.null(bars)){
    maxbar=max(bars)
    idx=match(names(bars),tree$tip.label)
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = -ymax*0.15,ytop=-ymax*0.15+ymax*0.1*bars/maxbar,col = "grey")
    
  }
  tree$ymax=ymax
  tree$vspace.reserve=vspace.reserve
  tree
}

get_all_node_children=function(node,tree){
  children=tree$edge[which(tree$edge[,1]==node),2]
  offspring=children
  for(child in children){
    offspring=c(offspring,get_all_node_children(child,tree))
  }
  offspring
}
get_node_children=function(node,tree){
  tree$edge[which(tree$edge[,1]==node),2]
}

get_samples_in_clade=function(node,tree){
  if(node<=length(tree$tip.label)){
    return(tree$tip.label[node])
  }
  tree$tip.label[intersect(get_all_node_children(node,tree),1:length(tree$tip.label))]
}

get_y_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$top-tree$coords[idx,c("a1","a0")])
}

get_x_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$coords[idx,c("b1","b0")])
}

get_parents=function(node,mut_table,exclude_root=TRUE){
  idx=which(mut_table[,2]==node)
  parents=node##Include the node
  while(length(idx)>0){
    if(length(idx)>1){
      stop("multiple parents!")
    }
    parent=mut_table[idx,1]
    parents=c(parents,parent)
    idx=which(mut_table[,2]==parent)
  }
  if(exclude_root){
    parents[-length(parents)]
  }else{
    parents
  }
}
##The following ties the PDX datastructure to the

get_idx_for_node=function(pdx,node){
  which(pdx$details$node==node)
}
###


get_edge_info=function(pdx,tree,node){
  y=get_y_range(tree,node)
  x=get_x_range(tree,node)
  idx=get_idx_for_node(pdx,node)
  samples=get_samples_in_clade(node,tree)
  list(yb=y[1],yt=y[2],x=x[1],xm=x[2],idx.in.details=idx,samples=samples)
}

add_annotation=function(pdx,tree,annot_function,control=NULL){
  N=dim(tree$edge)[1]
 # browser()
  lapply(1:N,function(i) annot_function(pdx,tree,tree$edge[i,2],control=control))
}



add_defaults=function(control,defaults,mandatory_fields=c()){
  if(is.null(control)){
    control=list()
  }
  for(field in mandatory_fields){
    if(is.null(control[[field]])){
      stop(sprintf("Required parameter %s not supplied in control list",field))
    }
  }

  for(field in names(defaults)){
    if(is.null(control[[field]])){
      control[[field]]=defaults[[field]]
    }
  }
  control
}


node_labels=function(tree,col="blue",cex=1,b_include_tips=FALSE){
  idx=which(tree$edge[,2]>length(tree$tip.label))
  if(b_include_tips){
    idx=1:dim(tree$edge)[1]
  }else{
    idx=which(tree$edge[,2]>length(tree$tip.label))
  }
  text(tree$edge[idx,2],x=tree$coords$b1[idx],y=tree$top-tree$coords$a1[idx],col=col,cex=cex)
}


add_binary_proportion=function(pdx,##<< PDX object or list including details matrix
                               tree,##<< enhanced phylo returned from plot_tree
                               node,##<< Tree node - maps to pdx$details$node
                               control,##<< Control parameters.. Requires bfield
                               ...
                               ){
  control=add_defaults(control,defaults=list( b.add.line=TRUE,
                                              b.add.text=FALSE),
                       mandatory_fields="bfield")

  bfield=control$bfield
  #browser()
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(pdx,tree,node)
  bdat=pdx$details[[bfield]][info$idx]
  if(is.null(bdat) || class(bdat)!="logical"){
    stop("Error in provided bfield (does it exist and is it boolean?)")
  }
  pass=sum(bdat,na.rm=TRUE)
  fail=sum(!bdat,na.rm=TRUE)
  tot=pass+fail
  ycross=info$yb+(fail/tot)*(info$yt-info$yb)
  ##Could add in a third category NA
  #missing=sum(is.na(bdat))
  if(control$b.add.line){
    arrows(y0=info$yb,y1=ycross,x0=info$x,length = 0,col="red",lend=1,...)
    arrows(y0=ycross,y1=info$yt,x0=info$x,length = 0,col="blue",lend=1,...)
  }
  if(control$b.add.text){
    text(y=ycross,x=info$x,label=pass,pos = 4,offset = 0,...)
  }
}





##Chromosome plot
test_chr=function(col_ctr,ptarget){
  tot=sum(col_ctr)
  if(tot>20){
    chires=chisq.test(col_ctr,p=ptarget,simulate.p.value = TRUE,B=500)
    pval=chires$p.value
    if(is.na(pval)){
      browser()
    }
    idx=which.max(abs(chires$stdres))
    stdres=chires$stdres[idx]
  }else{
    pval=1
    idx=-1
    stdres=0
  }
  list(idx.max.deviation=idx,pval=pval,stdres=stdres)
}

get_chr_cols=function(){
  cc=hue_pal(h = c(0, 360) + 15, c = 100, l = 65, h.start = 0,direction = 1)(24)
  cc[as.vector(rbind(1:12,13:24))]##Alternate
}

add_chromosome_decomp=function(pdx,##<< PDX object or list including details matrix
                               tree,##<< enhanced phylo returned from plot_tree
                               node,##<< Tree node - maps to pdx$details$node
                               control,##<< Control parameters.. Requires bfield
                               ...)
{
                               ##tree,node,res,ptarget){
  info=get_edge_info(pdx,tree,node)
  control=add_defaults(control,defaults=list( ),
                       mandatory_fields=c("targetprop"))
  ptarget=control$targetprop

  legend=FALSE
  if(is.null(tree$direction)){
    stop("Need to add to existing plotted tree")
  }
  chrs=sprintf("%s",c(1:22,"X","Y"))
  chr=match(pdx$details$Chrom[info$idx.in.detail],chrs)
  col_ctr=tabulate(chr,nbins=24)
  cols=get_chr_cols()
  tot=sum(col_ctr)
  y0=info$yb
  y1=info$yt
  x=info$x
  mh=y1-y0
  start=y1-0.05*mh
  #browser()
  unit=0.9*mh/tot
  at=c()

  chrtest=test_chr(col_ctr,ptarget)

  for(i in 1:length(col_ctr)){
    if(chrtest$idx.max.deviation==i){
      y0=start
      y1=start-unit*col_ctr[i]
      mh=y1-y0
    }
    if(col_ctr[i]>0){
      rect(xleft=x-0.25,xright=x+0.25,ytop = start,ybottom=start-unit*col_ctr[i],col=cols[i],border="black")
      start=start-unit*col_ctr[i]

    }

  }

  if(chrtest$pval<0.01){
    ##idx=which.max(abs(chires$stdres))
    ccol="blue"
    voff=(runif(1)-0.5)*0.5
    ##browser()
    yy=y0+mh*0.5+voff*mh
    ##text(x=x-0.25,y=tree$top-((y0+y1)/2),labels = sprintf("chr%s",col_ctr[idx]/sum(col_ctr)*ptarget[idx],chrs[idx],pval),cex = 1,pos = 2,offset = 0,col=ccol)
    #text(x=x-0.25,y=tree$top-((y0+y1)/2),labels = sprintf("%schr%s",ifelse(chires$stdres>0,"+","-"),chrs[idx]),cex = 1,pos = 2,offset = 0,col=ccol)
    txtbox(x=x-0.25,y=yy,txt = sprintf("%schr%s\nP=%3.2g",ifelse(chrtest$stdres>0,"+","-"),chrs[chrtest$idx],chrtest$pval),ccex = 1,pos = NULL,offset = 0,col=ccol)
    points(x,yy,cex=1,pch=4,col=ccol)
  }
  col_ctr
}







get_consequence_scheme=function(allowed.value=c("missense","nonsense","ess_splice","frameshift","inframe","cna","loh","start_lost","nc_ess_splice")){
  mut.scheme=read.table("consequence_col_scheme.txt",head=T,stringsAsFactors = FALSE)
  colnames(mut.scheme)=c("value","col","pch")
  if(length(allowed.value)==0){
    mut.scheme
  }else{
    mut.scheme %>% filter(value %in% allowed.value)
  }
}


add_simple_labels=function(
                    pdx,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
                    tree,##<< enhanced phylo returned from plot_tree
                    node,##<< Node (see details)
                    control,
                    query.allowed.df=get_consequence_scheme(),
                    ... ##<< paremeters for points (not color)
){
  control=add_defaults(control,defaults=list( query.field="VC",
                                              label.field="GENE",
                                              query.allowed.df=query.allowed.df,
                                              cex.label=1,
                                              b.add.label=TRUE,
                                              b.add.marker=TRUE)
                       )
  info=get_edge_info(pdx,tree,node)
  ##browser()
  details=pdx$details
  query.field=control$query.field
  query.allowed.df=control$query.allowed.df
  label.field=control$label.field
  idx=info$idx[which(details[[query.field]][info$idx] %in% query.allowed.df$value)]
  if(length(idx)>50){
    stop("too many variants to annotate")
  }
  query.value=details[[query.field]][idx]
  idx.match=match(query.value,query.allowed.df$value)
  cols=query.allowed.df$col[idx.match]
  pch=query.allowed.df$pch[idx.match]

  vlabels=details[[label.field]][idx]
  ## spread out
  N=length(idx)
  ##Vertical offset so that labels sit slightly above the markers.
  voffset=0.0075*(par("usr")[4]-par("usr")[3])
  if(N>0){
    yd=info$yt-info$yb
    if(N==1){
      y=0.5*(info$yb+info$yt)
    }else{
      y=seq(info$yb+(1/(N+1))*yd,info$yt-(1/(N+1))*yd,length.out = N)
    }
    if(control$b.add.marker){
      points(rep(info$x,N),y,col=cols,pch=pch,...)
    }
    if(control$b.add.label){
      text(rep(info$x,N),y+voffset,labels = vlabels,pos = 2,offset = 0,cex=control$cex.label)
    }
  }
  list(node=node,value=query.value)
}


add_specified_labels=function(
  pdx,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
  tree,##<< enhanced phylo returned from plot_tree
  node,##<< Node (see details)
  control,
  ... ##<< paremeters for points (not color)
){
  info=get_edge_info(pdx,tree,node)
  thisnode=node
  if(is.null(control$vars)){
    return(NULL)
  }
  vars=control$vars %>% filter(node==thisnode)
  ###browser()
  if(dim(vars)[1]==0){
    return(NULL)
  }
  if(control$b.annot.just.shared && thisnode<=length(tree$tip.label)){
    return(NULL)
  }
  cols=vars$col
  pch=vars$pch
  vlabels=vars$label
  #vlabels=details[[label.field]][idx]
  ## spread out
  N=length(vlabels)
  ##Vertical offset so that labels sit slightly above the markers.
  voffset=0.0075*(par("usr")[4]-par("usr")[3])
  y=c()
  if(N>0){
    yd=info$yt-info$yb
    ##browser()
    if(N==1){
      y=info$yb+rbeta(1,shape1 = 10,shape2 = 10)*yd
    }else{
      y=seq(info$yb+(1/(N+1))*yd,info$yt-(1/(N+1))*yd,length.out = N)
    }
    if(control$b.add.marker){
      points(rep(info$x,N),y,col=cols,pch=pch,bg=cols,...)
    }
    if(control$b.add.label){
      text(rep(info$x,N),y+voffset,labels = vlabels,pos = 2,offset = 0,cex=control$cex.label,col=control$col.label)
    }
  }
  list(node=node,value=N,y=y,x=info$x)
}

add_length=function(
  pdx,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
  tree,##<< enhanced phylo returned from plot_tree
  node,##<< Node (see details)
  control,
  ... ##<< paremeters for points (not color)
){
  info=get_edge_info(pdx,tree,node)
  thisnode=node
  #browser()
  ll=tree[[control$el_param]]
  if(is.null(ll)){
    stop(sprintf("tree$%s does not exist",control$el_param))
  }
  len=ll[which(tree$edge[,2]==node)]
  y=info$yb+(0.2+0.6*runif(1))*(info$yt-info$yb)
  if(len<control$cutoff && length(info$samples)>=control$min.shared){
  text(info$x,y,labels = sprintf("%3.1f",len),
       pos = 2,offset = 0,
       cex=control$cex.annot.label)
  }
  list(node=node,value=len)
}


add_vaf=function(pdx,##<< PDX object or list including details matrix
                 tree,
                 node,
                 control,
                 ...
){
  
  control=add_defaults(control,defaults=list( samples=c(),
                                              b.plot.bars=TRUE,
                                              lwd.rect=1,
                                              min.depth=1,
                                              filter.on=NULL,
                                              cex.label=1,
                                              min.mean.vaf=0.35,
                                              mode="recap"))
  ##browser()
  if(control$mode=="recap"){
    barcol="red"
    txtcol="blue"
  }else{
    barcol="darkgrey"
    txtcol="black"
  }
  
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(pdx,tree,node)
  
  if(length(info$idx)==0){
    return(NULL)
  }
  
  if(!is.null(control$filter.on)){
    info$idx=info$idx[which(pdx$details[[control$filter.on]][info$idx])]
    if(length(info$idx)==0){
      return(NULL)
    }
  }
  
  
  if(control$b.plot.bars){
    plotF=plotBars
  }else{
    plotF=plotPie
  }
  samples=control$samples
  if(is.null(samples) || length(samples)==0){
    samples=info$samples
  }
  
  if(length(samples)>1){
    if(length(info$idx)>1){
      df=data.frame(mtr=rowSums(pdx$mtr[info$idx,samples],na.rm = TRUE),
                    dep=rowSums(pdx$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }else{
      df=data.frame(mtr=sum(pdx$mtr[info$idx,samples],na.rm = TRUE),
                    dep=sum(pdx$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }
    adep=rowSums(pdx$dep[,samples])
  }else{
    
    df=data.frame(mtr=pdx$mtr[info$idx,samples],
                  dep=pdx$dep[info$idx,samples],stringsAsFactors = FALSE)
    adep=pdx$dep[,samples]
  }
  
  
  mdep=mean(df$dep)
  depth_zscore=sqrt(length(df$dep))*(mdep-mean(adep))/sd(adep)
  
  df=cbind(df,pdx$details[info$idx,])
  df=df[which(df$dep>=control$min.depth),]
  
  df$vaf=df$mtr/df$dep
  df=df[which(!is.na(df$vaf)),]
  N=dim(df)[1]
  if(N==0){
    return(df)
  }
  
  
  df=df[order(df$vaf),]
  yd=info$yt-info$yb
  
  
  if(N==1){
    y=0.5*(info$yb+info$yt)
    width=yd
  }else{
    y1=seq(info$yb,info$yt,length.out = N+2)
    #Don't use the ends..
    y=y1[2:(N+1)]
    width=y[2]-y[1]
  }
  
  if(!control$b.plot.bars){
    r=0.4  ##
  }else{
    r=0.4
  }
  for(i in 1:N){
    vaf=min(df$vaf[i],0.999)
    #if(!is.na(df$label[i]) && nchar(df$label[i])>0){
    #  text(x=info$x-r,y=y[i],labels = df$label[i],cex=control$cex.label,pos=2,offset=0)
    #}
    if(is.na(vaf)){
      plotF(x=info$x,y = y[i],radius=r,col=c("lightgray","lightgray"),prop=c(0,1),border="lightgray",width=width)
    }else{
      plotF(x=info$x,y = y[i],radius = r,col=c(barcol,"white"),prop = c(vaf,1-vaf),width = width)
    }
  }
  for(i in 1:N){
    if(!is.null(df$label) && !is.na(df$label[i]) && nchar(df$label[i])>0){
      text(x=info$x-r,y=y[i],labels = df$label[i],cex=control$cex.label,pos=2,offset=0)
      points(x=info$x,y=y[i],pch=df$pch[i],col=df$col[i],bg=df$col[i])
      #cat(df$label[i],":",df$vaf[i],": rank=",length(which(df$vaf[i]>df$vaf))/length(df$vaf),":",df$vaf,"\n")
    }
  }
  
  
  if( !control$b.plot.bars){
    return(df)
  }
  ##Now check to see if we need to highlight
  ##Test if VAF is significantly > 0.05 or significantly < 0.45
  ##Can also do a binomial test...
  MTR=sum(df$mtr)
  DEP=sum(df$dep)
  if(control$mode=="recap"){
    rect(xleft=info$x-r,xright=info$x+r,ybottom = info$yb, ytop=info$yt,border="darkgrey",lwd=control$lwd.rect)
    if(MTR/DEP>0.01){
      txt=gsub("^0\\.",".",sprintf("%3.2f",MTR/DEP))
      text(txt,x=info$x,y=info$yb+0.3*(info$yt-info$yb),col=txtcol,cex=0.6)
    }
    
  }else{
    min.mean.vaf=control$min.mean.vaf
    z=binom.test(MTR,DEP,alternative = "less",p=min.mean.vaf)
    z2=binom.test(MTR,DEP,alternative = "greater",p=0.05)
    z$p.value=max(z$p.value,z2$p.value)
    txt=gsub("^0\\.",".",sprintf("%3.2f",MTR/DEP))
    if(z$p.value<0.05){
      if(z$p.value<0.05/dim(tree$edge)[1]){
        border.color="red"
      }else{
        border.color="blue"
      }
    }else{
      border.color="darkgrey"
    }
    
    rect(xleft=info$x-r,xright=info$x+r,ybottom=y[1]-width/2,ytop=y[N]+width/2,border=border.color,lwd=control$lwd.rect)
    if(border.color!="darkgrey"){
      text(txt,x=info$x,y=y[1]+0.3*(y[N]-y[1]),col="black",cex=0.6)
    }
    
    if(!is.na(depth_zscore)){
      p=pnorm(depth_zscore,lower.tail=T)
      if(p<0.01){
        text(sprintf("%3.1f",mdep),x=info$x,y=y[1]+0.6*(y[N]-y[1]),col="red",cex=0.6)
        #text(sprintf("%2.1f",z),x=coords$b1,y=yy[1]+0.8*(yy[A]-yy[1]),col="red",cex=0.6)
      }
    }
  }
  arrows(x0=info$x,y0=info$yb,y1=info$yt,lwd=0.5,col="black",length=0,lend=2)
  
  df
  
}



add_vaf_old=function(pdx,##<< PDX object or list including details matrix
                 tree,
                 node,
                 control,
                 ...
){

  control=add_defaults(control,defaults=list( samples=c(),
                                              b.plot.bars=TRUE,
                                              lwd.rect=1,
                                              min.depth=1,
                                              filter.on=NULL,
                                              cex.label=1))
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(pdx,tree,node)
  
  if(length(info$idx)==0){
    return(NULL)
  }

  if(!is.null(control$filter.on)){
    info$idx=info$idx[which(pdx$details[[control$filter.on]][info$idx])]
    if(length(info$idx)==0){
      return(NULL)
    }
  }


  if(control$b.plot.bars){
    plotF=plotBars
  }else{
    plotF=plotPie
  }
  samples=control$samples
  if(is.null(samples) || length(samples)==0){
    samples=info$samples
  }
  if(length(samples)>1){
    if(length(info$idx)>1){
      df=data.frame(mtr=rowSums(pdx$mtr[info$idx,samples],na.rm = TRUE),
                    dep=rowSums(pdx$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }else{
      df=data.frame(mtr=sum(pdx$mtr[info$idx,samples],na.rm = TRUE),
                    dep=sum(pdx$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }
  }else{
    df=data.frame(mtr=pdx$mtr[info$idx,samples],
                  dep=pdx$dep[info$idx,samples],stringsAsFactors = FALSE)
  }
  df=cbind(df,pdx$details[info$idx,])
  df=df[which(df$dep>=control$min.depth),]




  df$vaf=df$mtr/df$dep
  df=df[which(!is.na(df$vaf)),]
  N=dim(df)[1]
  cat(N,"\n")
  if(N==0){
    return(df)
  }


  df=df[order(df$vaf),]
  yd=info$yt-info$yb


  if(N==1){
    y=0.5*(info$yb+info$yt)
    width=yd
  }else{
    y1=seq(info$yb,info$yt,length.out = N+2)
    #Don't use the ends..
    y=y1[2:(N+1)]
    width=y[2]-y[1]
  }

  if(!control$b.plot.bars){
    r=0.8  ##r>0.5 will cause likely overlap problems
  }else{
    r=0.4
  }
  for(i in 1:N){
    vaf=min(df$vaf[i],0.999)
    #if(!is.na(df$label[i]) && nchar(df$label[i])>0){
    #  text(x=info$x-r,y=y[i],labels = df$label[i],cex=control$cex.label,pos=2,offset=0)
    #}
    if(is.na(vaf)){
      plotF(x=info$x,y = y[i],radius=r,col=c("lightgray","lightgray"),prop=c(0,1),border="lightgray",width=width)
    }else{
      plotF(x=info$x,y = y[i],radius = r,col=c("black","white"),prop = c(vaf,1-vaf),width = width)
    }
  }
  for(i in 1:N){
    if(!is.null(df$label) && !is.na(df$label[i]) && nchar(df$label[i])>0){
      text(x=info$x-r,y=y[i],labels = df$label[i],cex=control$cex.label,pos=2,offset=0)
      points(x=info$x,y=y[i],pch=df$pch[i],col=df$col[i],bg=df$col[i])
      #cat(df$label[i],":",df$vaf[i],": rank=",length(which(df$vaf[i]>df$vaf))/length(df$vaf),":",df$vaf,"\n")
    }
  }


  if( !control$b.plot.bars){
    return(df)
  }
  ##Now check to see if we need to highlight
  ##Test if VAF is significantly > 0.05 or significantly < 0.45
  ##Can also do a binomial test...
  MTR=sum(df$mtr)
  DEP=sum(df$dep)
  min.mean.vaf=0.45
  z=binom.test(MTR,DEP,alternative = "less",p=min.mean.vaf)
  z2=binom.test(MTR,DEP,alternative = "greater",p=0.05)
  z$p.value=max(z$p.value,z2$p.value)
  txt=gsub("^0\\.",".",sprintf("%3.2f",MTR/DEP))
  if(z$p.value<0.05){
    if(z$p.value<0.05/dim(tree$edge)[1]){
      border.color="red"
    }else{
      border.color="blue"
    }
  }else{
    border.color="darkgrey"
  }

  rect(xleft=info$x-r,xright=info$x+r,ybottom=y[1]-width/2,ytop=y[N]+width/2,border=border.color,lwd=control$lwd.rect)
  if(border.color!="darkgrey"){
    text(txt,x=info$x,y=y[1]+0.3*(y[N]-y[1]),col="black",cex=0.6)
  }
  arrows(x0=info$x,y0=info$yb,y1=info$yt,lwd=0.5,col="black",length=0,lend=2)

  df

}

add_colored_vaf=function(pdx,##<< PDX object or list including details matrix
                         tree,
                         node,
                         control,
                         ...
){
  control=add_defaults(control,defaults=list( samples=c(),
                                              lwd.rect=1,
                                              min.depth=1,
                                              filter.on=NULL,
                                              mean.type="agg",
                                              cex.label=1,
                                              cex.vaf.label=0.6,
                                              mode="recap",
                                              b.add.vaf.label=TRUE))
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(pdx,tree,node)
  ##browser()
  
  if(length(info$idx)==0){
    return(NULL)
  }
  
  if(!is.null(control$filter.on)){
    info$idx=info$idx[which(pdx$details[[control$filter.on]][info$idx])]
    if(length(info$idx)==0){
      return(NULL)
    }
  }
  
  samples=control$samples
  if(is.null(samples) || length(samples)==0){
    samples=info$samples
  }
  
  ##browser()
  if(length(samples)>1){
    if(length(info$idx)>1){
      df=data.frame(mtr=rowSums(pdx$mtr[info$idx,samples],na.rm = TRUE),
                    dep=rowSums(pdx$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }else{
      df=data.frame(mtr=sum(pdx$mtr[info$idx,samples],na.rm = TRUE),
                    dep=sum(pdx$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }
  }else{
    df=data.frame(mtr=pdx$mtr[info$idx,samples],
                  dep=pdx$dep[info$idx,samples],stringsAsFactors = FALSE)
  }
  df=cbind(df,pdx$details[info$idx,])
  if(control$mean.type!="agg"){
    df=df[which(df$dep>=control$min.depth),]
  }
  df$vaf=df$mtr/df$dep
  df=df[which(!is.na(df$vaf)),]
  vaf=mean(df$vaf)
  
  
  N=dim(df)[1]
  if(N==0){
    return(df)
  }
  
  # cat("Fix me in add_color_vaf\n")
  r=0.3 ##0.25*length(tree$tip.label)/20 #0.3
  #rect(xleft=info$x-r,xright=info$x+r,ybottom = info$yb, ytop=info$yt,
  #     col=thiscol)
  
  ##Test if VAF is significantly > 0.05 or significantly < 0.45
  ##Can also do a binomial test...
  MTR=sum(df$mtr)
  DEP=sum(df$dep)
  
  if(control$mean.type=="agg"){
    vaf=MTR/DEP
  }
  thiscol=getVafColor(vaf)$col
  if(control$mean.type=="agg" & DEP<control$min.agg.depth){
    thiscol="grey"
  }
  #thiscol2=getVafColor(MTR/DEP)$col
  
  rect(xleft=info$x-r,xright=info$x+r,ybottom = info$yb, ytop=info$yt,
       col=thiscol)
  if(control$mode=="recap"){
    if(MTR/DEP>0.01){
      txt=gsub("^0\\.",".",sprintf("%3.2f",MTR/DEP))
      if(control$b.add.vaf.label){
        text(txt,x=info$x,y=info$yb+0.3*(info$yt-info$yb),col="black",cex=control$cex.vaf.label)
      }
    }
    rect(xleft=info$x-r,xright=info$x+r,ybottom = info$yb, ytop=info$yt,border="darkgrey",lwd=control$lwd.rect)
  }else{
    min.mean.vaf=0.45
    z=binom.test(MTR,DEP,alternative = "less",p=min.mean.vaf)
    z2=binom.test(MTR,DEP,alternative = "greater",p=0.05)
    z$p.value=max(z$p.value,z2$p.value)
    txt=gsub("^0\\.",".",sprintf("%3.2f",MTR/DEP))
    if(z$p.value<0.05){
      if(z$p.value<0.05/dim(tree$edge)[1]){
        border.color="red"
      }else{
        border.color="blue"
      }
    }else{
      border.color="darkgrey"
    }
    
    rect(xleft=info$x-r,xright=info$x+r,ybottom = info$yb, ytop=info$yt,border=border.color,lwd=control$lwd.rect)
    if(border.color!="darkgrey"){
      text(txt,x=info$x,y=info$yb+0.3*(info$yt-info$yb),col="black",cex=0.6)
    }
  }
  
  
  
  
  #arrows(x0=info$x,y0=info$yb,y1=info$yt,lwd=0.5,col="black",length=0,lend=2)
  df
}


plotBars=function(x,y,radius,col,prop,border="black",width=1){
  #cat(prop,"\n")
  if(prop[1]<0.001){
    return(1)
  }
  if(width<diff(par("usr")[c(2,4)])/200){
    arrows(x0 = x-radius,y0=y,x1=x-radius+2*radius*prop[1],col=col,lend=2,length=0)
    arrows(x0 = x-radius+2*radius*prop[1],y0=y,x1=x-radius+2*radius,col=rgb(0.98,0.98,0.98),lend=2,length=0)
  }else{
    rect(xleft = x-radius,xright =x-radius+2*radius*prop[1],ybottom = y-width/2,ytop=y+width/2,border = NA,col=col)
    rect(xleft =  x-radius+2*radius*prop[1],xright =x+radius,ybottom = y-width/2,ytop=y+width/2,border = NA,col=rgb(0.98,0.98,0.98))
  }
  1
}

plotBars2=function(x,y,radius,col,prop,border="black",width=1){

    rect(xleft = x-radius,xright =x-radius+2*radius*1,ybottom = y-width/2,ytop=y+width/2,border = NA,col=rgb(red = prop,green=(1-prop),blue = 0.01))

  1
}


plotPie=function(x,y,radius,col,prop,llwd=0.5,border="black",width=1){
  lims=par("usr")
  as=dev.size()
  asr=as[1]/as[2]
  yscale=asr*(lims[4]-lims[3])/(lims[2]-lims[1])
  prop=prop/sum(prop)
  cutpoint=c(0,cumsum(prop)*2*pi)
  max.y.radius=width/2
  if(radius*yscale>max.y.radius){
    radius=max.y.radius/yscale
  }

  N=2*pi/0.05
  n=ceiling(N*diff(cutpoint)/(2*pi))
  d=diff(cutpoint)/n
  if(length(prop)>1){
    for(i in 2:length(cutpoint)){
      polygon(x+c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
              y+yscale*c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
              border=border,col=col[i-1],lwd=llwd)
    }
  }else{
    i=2
    polygon(x+c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
            y+yscale*c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
            border=border,col=col[1],lwd=llwd)
  }
  yscale
}


##Gets unique colour pch combos and returns in dataframe with columns "col" and "pch"
get_color_pch_df=function(n){
  pch.list=c(18,17,16,15,0:6)
  if(n>length(pch.list)*8){
    stop("Too many colours requested")
  }
  cols=rep(RColorBrewer::brewer.pal(8,"Set1"),times=length(pch.list))
  pch=rep(pch.list,each=8)
  data.frame(col=cols,pch=pch,stringsAsFactors = FALSE)[1:n,]

}

get_qdf=function(values){
  if(length(values)>length(unique(values))){
    stop("get_qdf: please provide values without duplication")
  }
  cbind(data.frame(value=values,stringsAsFactors = FALSE),
        get_color_pch_df(length(values)))
}

plot_tree_labels_consequence=function(tree,details,consequences,
                                      query.allowed.df=get_qdf(consequences),
                          query.field="VC",
                          label.field="GENE",
                          cex.label=1){
  plot_tree_labels(tree,details,
                   query.allowed.df = query.allowed.df,
                   query.field=query.field,
                   label.field=label.field,
                   cex.label=cex.label)
}


plot_tree_labels=function(tree,pdx,
                          query.field="VC",
                          query.allowed.df=get_consequence_scheme(),
                          label.field="GENE",
                          cex.label=1){

  ##annot_function(pdx,tree,tree$edge[i,2],control=control)
  control=add_defaults(list(),defaults=list( query.field=query.field,
                                              label.field=label.field,
                                              query.allowed.df=query.allowed.df,
                                              cex.label=cex.label,
                                              b.add.label=TRUE,
                                              b.add.marker=TRUE)
  )
  res=add_annotation(pdx,
                     tree,
                     add_simple_labels,control=control)

  with(control$query.allowed.df,legend("topleft",legend=value,col=col,pch=pch))
}

add_vertical_labels=function(
  pdx,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
  tree,##<< enhanced phylo returned from plot_tree
  node,##<< Node (see details)
  control,
  ... ##<< paremeters for points (not color)
){
  info=get_edge_info(pdx,tree,node)
  thisnode=node
  if(is.null(control$vars)){
    return(NULL)
  }
  vars=control$vars %>% filter(node==thisnode)
  #browser()
  if(dim(vars)[1]==0){
    return(NULL)
  }
  scheme=get_driver_scheme2()
  vars$group=gsub(":.+","",vars$label2)
  #col=scheme$colour[match(vars$label2,scheme$driver)]
  col=scheme$colour[match(vars$group,scheme$group)]
  pch=vars$pch
  
  vlabels=vars$label
  if(control$b.show.just.group){
    vlabels=unique(vars$group)
  }
  N=length(vlabels)
  
  col1=col
  col=unique(col)
  N1=length(col)
  if(N1>4){
    stop("too many distinct driver types <= 4. Please fix code..")
  }
  ##Vertical offset so that labels sit slightly above the markers.
  voffset=0.0075*(par("usr")[4]-par("usr")[3])
  ww=min(0.25,0.25*length(tree$tip.label)/20)   ##half Rect width   ## Or less if very small tree
  offsetodd=c(0,-0.5,0.5)*ww/0.25
  offseteven=c(0.25,-0.25,0.75,-0.75)*ww/0.25
  #idx2=info$idx.in.details
  
  if(N1 %% 2){
    offset=offsetodd[1:N1]
    
  }else{
    offset=offseteven[1:N1]
    ##oidx=(offset)
  }
  offset=sort(offset)
  #browser()
  
  if(N>0){
    yd=info$yt-info$yb
    if(N==1){
      #y=0.5*(info$yb+info$yt)
      y=info$yb+rbeta(1,shape1 = 5,shape2 = 5)*yd
    }else{
      y=seq(info$yb+(1/(N+1))*yd,info$yt-(1/(N+1))*yd,length.out = N)
    }
    
    if(control$b.add.label){
      text(rep(info$x+min(offset)-ww,N),y+voffset,labels = vlabels,pos = 2,offset = 0,cex=control$cex.label,col=col1)
    }
    for(i in 1:N1){
      rect(xleft=info$x+offset[i]-ww,xright=info$x+offset[i]+ww,
           ytop =info$yt,ybottom = info$yb,col=col[i],border=NA)
    }
  }
  list(node=node,value=N)
}
plot_labelled_tree=function(tree,
                            pdx,
                            cv=c("missense","nonsense","ess_splice","frameshift","inframe","loh","cna"),
                            cex.annot.label=0.5,legpos="topleft",style="classic",genes=NULL,seed=-1,b.annot.just.shared=FALSE,col.label="black",b.show.just.group=FALSE){
  if(seed>0){
    set.seed(seed)
  }else{
    seed=sample(100000,1)+100000
    set.seed(seed)
  }
  muts=get_all_tree_drivers(pdx,cv = cv,genes=genes)
  print(muts)
  
  scheme=get_consequence_scheme(allowed.value=cv)
  colnames(scheme)[1]="cv"
  vars=muts %>% left_join(scheme,by="cv")
  ##browser()
  control=add_defaults(list(),defaults=list(vars=vars,
                                            cex.label=cex.annot.label,
                                            b.add.label=cex.annot.label>0,
                                            b.add.marker=TRUE,
                                            b.annot.just.shared=b.annot.just.shared,
                                            col.label=col.label,
                                            b.show.just.group=b.show.just.group))
  if(style=="classic"){
    annotfun=add_specified_labels
  }else{
    annotfun=add_vertical_labels
  }
  res=add_annotation(pdx,
                     tree,
                     annotfun,control=control)
  
  if(!is.null(legpos) && !is.null(vars) && dim(vars)[1]>0){
    ##leg=legend(legpos,mut.scheme$cv,col = mut.scheme$col,pch=mut.scheme$pch,pt.bg=mut.scheme$col,pt.cex = cex.marker)$rect
    ##browser()
    leg=with(unique(control$vars[,c("cv","col","pch")]),legend(legpos,legend=cv,col=col,pch=pch,pt.bg=col)$rect)
  }
  tree$seed=seed
  tree
}

plot_chromosome_annotated_tree=function(tree,pdx){

  chrs=sprintf("%s",c(1:22,"X","Y"))
  chr=match(pdx$details$Chrom,chrs)
  col_ctr=tabulate(chr,nbins=24)

  ptarget=col_ctr/sum(col_ctr)
  ptarget=(ptarget+0.0001)/sum(ptarget+0.0001)
  control=list(targetprop=ptarget)
  res=add_annotation(pdx,
                     tree,
                     add_chromosome_decomp,control=control)

  leg=legend("topleft",sprintf("chr%s",c(1:22,"X","Y")),col = get_chr_cols(),pch=15,
             pt.cex=2)$rect
}


#add_vaf_bar
#add_vof_pie
#add_label

##Need to hack at get_all_tree_drivers to make a satisfactory interface for specifying
## what drivers to annotate the tree with!
get_all_tree_drivers=function(pdx,drivers=get_drivers("drivers_red.bed"),
                              cv=c("missense","nonsense","ess_splice","frameshift","inframe","loh","cna"),genes=NULL){
  drivers=get_tree_drivers(pdx,drivers=drivers,cv = cv,genes=genes)

  drivers$label=with(drivers,sprintf("%s:%s",GENE,HGVS_PROTEIN))
  ##map copy number
  cna=rbind(get_cna_profile(pdx),get_loh_profile(pdx))[,c("label","profile","cv")]
  drivers=rbind(drivers[,c("label","profile","cv")],cna)
  #drivers$node=pdx$tree_ml$edge[pdx$summary$edge_ml[match(drivers$profile,pdx$summary$profile)],2]
  drivers$node=get_node_for_profile(pdx,drivers$profile)
  
  drivers$label2=gsub("_[A-Z]$","",drivers$label)
  drivers$label2=gsub("_[a-z]$","",drivers$label2)
  drivers
}

get_node_for_profile=function(pdx,profile){
  ##df and the tree are aligned.
  df=pdx$df$df
  if(any(pdx$df$samples!=pdx$tree_ml$tip.label)){
    stop("profile df is inconsistent with tree!")
  }
  df$node=pdx$tree_ml$edge[,2]
  df$node[match(profile,df$profile)]
}


get_tree_drivers=function(pdx,b_just_drivers=FALSE,genes=c("EIF6","SBDS","TP53","RPL5"),
                          cv=c("missense","nonsense","ess_splice","frameshift",
                               "inframe","start_lost","nc_ess_splice"),
                          drivers=drivers){
  res=pdx
  res$dat$details$ID=with(res$dat$details,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))

  defdrivers=drivers
  if(is.null(res$dat$details$color)){
    res$dat$details$color="black"
  }
  ##Restrict to sites matching drivers
  res$dat$details$profile=res$df$df$profile[res$summary$edge_ml]
  allprofiles=unique(res$dat$details$profile)
  ##browser()
  if(b_just_drivers){
    idx.keep=which(res$dat$details$ID %in% defdrivers$ID)
  }else{
    if(length(genes)>0){
      idx.keep=which(res$dat$details$GENE %in% genes  & res$dat$details$VC %in% cv)
    }else{
      idx.keep=which(res$dat$details$VC %in% cv)
    }
  }
  if(length(idx.keep)==0){
    #cat("No drivers to add to tree!\n")
    return(data.frame(CHROM=character(),
                      POS=integer(),
                      REF=character(),
                      ALT=character(),
                      cv=character(),
                      GENE=character(),
                      HGVS_PROTEIN=character(),
                      VC=character(),
                      profile=character()))
  }
  #browser()
  drivers=res$dat$details[idx.keep, c("Chrom","Pos","Ref","Alt","color","GENE","VC","HGVS_PROTEIN","profile")]#
  colnames(drivers)[c(1:4,7)]=c("CHROM","POS","REF","ALT","cv")#
  drivers$label=with(drivers,sprintf("%s:%s",GENE,HGVS_PROTEIN))
  #drivers$node=res$tree_ml$edge[res$summary$edge_ml[match(drivers$profile,res$summary$profile)],2]
  drivers$node=get_node_for_profile(res,drivers$profile)
  drivers$label2=gsub("_[A-Z]$","",drivers$label)
  drivers$label2=gsub("_[a-z]$","",drivers$label2)

  drivers
}


get_drivers=function(driver_file){
  defdrivers=read.table(driver_file,head=FALSE,stringsAsFactors = FALSE)
  colnames(defdrivers)=c("CHROM","POS","REF","ALT")
  defdrivers$ID=sprintf("%s:%s:%s:%s",defdrivers$CHROM,defdrivers$POS,defdrivers$REF,defdrivers$ALT)
  ##defdrivers=defdrivers[which(defdrivers$CHROM=="ZZZ"),] ## HACK to ignore driver file for now

  defdrivers
}

##Invokes high level pdx....
plot_basic_tree=function(pdx,label,style="classic",cex.annot.label=1,legpos="topleft",genes=NULL,cv=c(),seed=-1,b.annot.just.shared=FALSE,col.label="black",b.show.just.group=FALSE,...){
  tree=plot_tree(pdx$tree_ml,...);title(label)
  tree=plot_labelled_tree(tree,pdx,style=style,cex.annot.label=cex.annot.label,legpos=legpos,genes=genes,cv=cv,seed=seed,b.annot.just.shared=b.annot.just.shared,col.label=col.label,
                          b.show.just.group = b.show.just.group)
  tree
}


txtbox=function(xt, yt, txt,ccex=1,bcol="white",...){
  sw   <- strwidth(txt,cex=ccex)
  sh   <- strheight(txt,cex=ccex)
  frsz <- 0.3
  rect(
    xt - 2*sw*(0.5+ frsz),
    yt - sh*(0.5+ frsz),
    xt ,
    yt + sh*(0.5+ frsz),col=bcol
  )
  text(x=xt -sw*(0.5+ frsz),y=yt,label=txt,cex=ccex,...)
}

plot_chromosome_annotated_tree_pdx=function(pdx,label){
tree=plot_tree(pdx$tree_ml);title(label)
plot_chromosome_annotated_tree(tree,pdx$dat)
}
plot_vaf_tree=function(pdx,label,samples=NULL,cex.label=1,b.add.drivers=TRUE,filter.on=NULL,min.depth=1,mode="qc",b.plot.bars=TRUE,lwd=1,min.mean.vaf=0.4,...){
  tree=plot_tree(pdx$tree_ml,cex.label = cex.label,mar=c(1,1,2,3)+0.1,...)
  title(label)
  ##mtext(label,side=3,xpd=TRUE)
  idx=c()
  if(b.add.drivers){
    
    drivers=get_tree_drivers(pdx,drivers=get_drivers("drivers_red.bed"))
    colnames(drivers)[1:4]=c("Chrom", "Pos", "Ref", "Alt")
    scheme=get_consequence_scheme()
    colnames(scheme)[1]="cv"
    drivers=drivers %>% left_join(scheme,by="cv")
    drivers$ID=with(drivers,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
    #browser()
    if(dim(drivers)[1]>0){
      pdx$dat$details=pdx$dat$details %>% left_join(drivers[,c("Chrom", "Pos", "Ref", "Alt","node","label", "label2","col","pch")],
                                                    by=c("Chrom", "Pos", "Ref", "Alt","node"))
      
      
    }
    #control=list(cex.label=0.5,samples=samples,filter.on=filter.on,min.depth=min.depth)
    #}else{
    # control=list(filter.on=filter.on,min.depth=min.depth)
  }else{
    idx=which(pdx$dat$details$keep_embryonic==1)
    pdx$dat$details$label=NA
    if(length(idx)>0){
      #pdx$dat$details$label[idx]=sprintf("emb%s",1:length(idx))
    }
  }
  control=list(cex.label=cex.label,samples=samples,filter.on=filter.on,min.depth=min.depth,mode=mode,b.plot.bars=b.plot.bars,lwd.rect=lwd,min.mean.vaf=min.mean.vaf)
  #browser()
  res=add_annotation(pdx$dat,tree,add_vaf,control=control)
  ##tree=plot_labelled_tree(tree,pdx,style="classic",cex.annot.label=1,legpos="topleft")
  tree$embryonic_idx=idx
  return(tree)
}

plot_vaf_treeO=function(pdx,label,samples=NULL,cex.label=1,b.add.drivers=TRUE,filter.on=NULL,min.depth=1){
  tree=plot_tree(pdx$tree_ml,cex.label = cex.label,mar=c(1,1,2,3)+0.1)
  title(label)
  if(b.add.drivers){

    drivers=get_tree_drivers(pdx,drivers=get_drivers("drivers_red.bed"))
    colnames(drivers)[1:4]=c("Chrom", "Pos", "Ref", "Alt")
    scheme=get_consequence_scheme()
    colnames(scheme)[1]="cv"
    drivers=drivers %>% left_join(scheme,by="cv")
    drivers$ID=with(drivers,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
    #browser()
    if(dim(drivers)[1]>0){
    pdx$dat$details=pdx$dat$details %>% left_join(drivers[,c("Chrom", "Pos", "Ref", "Alt","node","label", "label2","col","pch")],
                                                  by=c("Chrom", "Pos", "Ref", "Alt","node"))


    }
    #control=list(cex.label=0.5,samples=samples,filter.on=filter.on,min.depth=min.depth)
  #}else{
   # control=list(filter.on=filter.on,min.depth=min.depth)
  }
  control=list(cex.label=cex.label,samples=samples,filter.on=filter.on,min.depth=min.depth)
  #browser()
  res=add_annotation(pdx$dat,tree,add_vaf,control=control)
}

plot_color_vaf_tree=function(pdx,label,samples=NULL,cex.label=1,b.add.drivers=TRUE,
                             filter.on=NULL,min.depth=1,min.agg.depth=50,mean.type="agg",mode="recap",b.add.vaf.label=TRUE,cex.vaf.label=0.6,legpos="topleft",ymax=NULL){
  tree=plot_tree(pdx$tree_ml,cex.label = cex.label,mar=c(1,1,2,3)+0.1,ymax=ymax)
  title(label)
  control=list(cex.label=0.5,
               samples=samples,
               filter.on=filter.on,
               min.depth=min.depth,
               mean.type=mean.type,
               min.agg.depth=min.agg.depth,
               mode=mode,
               b.add.vaf.label=b.add.vaf.label,
               cex.vaf.label=cex.vaf.label
  )
  res=add_annotation(pdx$dat,tree,add_colored_vaf,control=control)
  if(b.add.drivers){
    tree=plot_labelled_tree(tree,pdx,style="classic",cex.annot.label=0.6,legpos="topleft")
  }
  add_vaf_legend(legpos)
}
add_vaf_legend=function(legpos){
  if(!is.null(legpos)){
    vafc=getVafColor()
    legend(legpos,vafc$strVaf,col=vafc$colR,pch=15,pt.cex =3,title="Mean VAF")
  }
}

plot_len=function(pdx,label,edge.len.param="edge.length",cex.label=1,cex.annot.label=0.5){
  tree=plot_tree(pdx$tree_ml,cex.label = cex.label,mar=c(1,1,2,3)+0.1)
  title(label)
  control=list(el_param=edge.len.param,cex.annot.label=cex.annot.label)
  #browser()
  res=add_annotation(pdx$dat,tree,add_length,control=control)
  tree
}



getVafColor=function(vaf=-1,brewer.palette="YlOrRd"){
  colRange=brewer.pal(8,brewer.palette)
  #colorRampPalette(c("yellow", "red"),bias=10)(5)
  vafRange=c(0,0.02,0.05,seq(0.1,0.5,length.out = 5))
  strVaf=sprintf("%3.2f-%3.2f",vafRange[-length(vafRange)],vafRange[-1])
  if(!is.na(vaf) && vaf>=0){
    thiscol=colRange[findInterval(vaf,vafRange,rightmost.closed = TRUE,all.inside = TRUE)]
  }else{
    thiscol="grey"
  }
  list(colR=colRange,vafR=vafRange,strVaf=strVaf,col=thiscol)
}

plot_pooled_vaf_tree=function(pdx,label,cex.label=1,b.add.drivers=TRUE,lwd=1,min.mean.vaf=0.45){
  plot_vaf_tree(pdx,label,samples=NULL,cex.label=cex.label,b.add.drivers=b.add.drivers,lwd=lwd,min.mean.vaf=min.mean.vaf)
}

plot_all_vaf=function(pdx,label,vtype="SNV",samples=colnames(pdx$dat$mtr),min.mean.vaf=0.45){
  pdx$dat$details$is_vtype=pdx$dat$details$TYPE %in% vtype
  ##Reorder so that samples are in tip order...
  samples=intersect(c(setdiff(samples,pdx$tree_ml$tip.label),setdiff(pdx$tree_ml$tip.label,"zeros")),samples)

  for(x in samples){
    plot_vaf_tree(pdx,
                  sprintf("%s: Annotated with VAF from %s\n Mean Depth=%3.2f [MinVaf=%3.2f]",label,x,mean(pdx$dat$dep[,x],na.rm=T),min.mean.vaf),
                  samples=x,filter.on="is_vtype",min.depth=8,min.mean.vaf=min.mean.vaf)
  }
}

add_heatmap=function(tree,heatmap,heatvals=NULL,border="black",cex.label=2,pos=4,force.count=-1,txtcols=rep("black",dim(heatmap)[1]),b.add.lines=FALSE){
  ymax=tree$ymax
  idx=match(colnames(heatmap),tree$tip.label)
  top=-0.05*ymax
  if(force.count<0){
    force.count=dim(heatmap)[1]
  }
  gap=ymax*tree$vspace.reserve/force.count
  labels=rownames(heatmap)
  for(i in 1:dim(heatmap)[1]){
    bot=top-gap
    if(!is.na(border) && border=="none"){
      ##collapse 
      cc=unique(heatmap[i,])
      for(ccc in cc){
        idxx=which(heatmap[i,]==ccc)
        start=idxx[1]
        end=start
        for(j in idxx){
          if(j>(end+1)){
            rect(xleft=start-0.5,xright=end+0.5,ybottom = bot,ytop=top,col = heatmap[i,end],border=NA)
            start=j
            end=j
          }else{
            end=j
          }
        }
        rect(xleft=start-0.5,xright=end+0.5,ybottom = bot,ytop=top,col = heatmap[i,end],border=NA)
      }
    }else{
      rect(xleft=idx-0.5,xright=idx+0.5,ybottom = bot,ytop=top,col = heatmap[i,],border=border,lwd=0,ljoin=2)
      if(TRUE){
        offsets=0.05
      ##draw lines for boundaries...
      ch=heatmap[i,1]
      k=1
      LWD=1
      LTY="dotted"
      if(b.add.lines && ch!="white"){
        segments(x0 = idx[k]-0.5,y0=bot,y1=tree$top,col=ch,lwd=LWD,lty=LTY)
      }
      for(k in 2:dim(heatmap)[2]){
        
        if(heatmap[i,k]!=ch){
          
          if(b.add.lines && ch!="white"){
            segments(x0 = idx[k]-0.5-offsets,y0=bot,y1=tree$top,col=ch,lwd=LWD,lty=LTY)
          }
          ch=heatmap[i,k]
          if (b.add.lines && ch!="white"){
            segments(x0 = idx[k]-0.5+offsets,y0=bot,y1=tree$top,col=ch,lwd=LWD,lty=LTY)
          }
        }
      }
      if(b.add.lines && ch!="white"){
        segments(x0 = idx[k]+0.5,y0=bot,y1=tree$top,col=ch,lwd=LWD,lty=LTY)
      }
      }
    }
    
    if(!is.null(heatvals)){
      text(xx=idx,y=0.5*(top+bot),labels = sprintf("%3.2f",heatvals[i,]))
    }
    if(!is.null(labels)){
      text(labels[i],x=par("usr")[1]+1,y=0.5*(top+bot),pos = pos,cex = cex.label,col=txtcols[i])
    }
    top=bot
  }
  tree
}


magnify_top_of_tree=function(tree,a=0.75,scale.factor=5){
  nh=nodeHeights(tree)
  ht=max(nh)
  f2=1;##(ht-scale.factor*a)/(ht-a)# rescale factor for rest of tree (exact for ultrametric)
  factor=ifelse(nh[,1]<=a & nh[,2]-nh[,1]>0,
                ifelse(nh[,2]>a,(scale.factor*(a-nh[,1])+f2*(nh[,2]-a))/(nh[,2]-nh[,1]),scale.factor),
                f2)
  tree$edge.length=factor*tree$edge.length
  tree
}

magnify_top_of_tree_pdx=function(pdx,a=0.75,scale.factor=5){
  pdx$tree_ml=magnify_top_of_tree(pdx$tree_ml,a=a,scale.factor = scale.factor)
  pdx
}

add_scale=function(tree,a,scale.factor,breaks1,breaks2,labels1=breaks1,labels2=breaks2,col1="grey",col2="black"){
  xa=length(tree$tip.label)*1.05
  xd=0.01*length(tree$tip.label)
  ##Add the first axis line
  segments(y0=tree$top-a*scale.factor,y1=tree$top,x0 = xa,col=col1,lwd=2,xpd=NA)
  ##Add the divider between scales (can be made more fancy)
  segments(y0=tree$top-a*scale.factor,x0 = xa-xd,x1=xa,lwd=3,col="black",xpd=NA)
  ##Add division on first axis
  segments(y0=tree$top-breaks1*scale.factor,x0=xa-xd,x1=xa,col=col1,xpd=NA)
  ##Add labels first axis
  text(y=tree$top-breaks1*scale.factor,x=xa,pos=4,labels=labels1,xpd=NA)
  ##Add second axis
  segments(y1=tree$top-a*scale.factor,y0=0,x0 = xa,col=col2,lwd=2,xpd=NA)
  ##Add divider
  segments(y0=tree$top-a*(scale.factor-1)-breaks2,x0=xa-xd,x1=xa,col=col2,xpd=NA)
  cat(a*scale.factor+breaks2)
  ##Add labels
  text(y=tree$top-a*(scale.factor-1)-breaks2,x=xa,pos=4,labels=labels2,xpd=NA)
  
}
