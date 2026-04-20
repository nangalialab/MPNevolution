##Utilities
require("MutationalPatterns")
require("metafor")
require("lubridate")
require("data.table")
require("dndscv")
if(!exists("SPECIES_BUILD")){
  stop("Please create SPECIES_BUILD and set to e.g. hg19,mm10 (or load project.R first)")
}

get_ref_genome_name=function(){
  if(grepl("Hglaber",SPECIES)){
    ref_genome= sprintf("BSgenome.%s.NCBI.%s",SPECIES,"Nakedmoleratpaternal")
  }else{
    ref_genome= sprintf("BSgenome.%s.UCSC.%s",SPECIES,SPECIES_BUILD)
  }
}

chk_ref_genome_name=get_ref_genome_name()
if(!exists("ref_genome") || ref_genome != chk_ref_genome_name){
    ref_genome= get_ref_genome_name()
    library(ref_genome, character.only = TRUE)
    if(grepl("Hglaber",SPECIES)){
      gene_rds="/nfs/casm/team273jn/resources/heterocephalus_glaber_male_core_111_1_GENES.RDS"
      if(file.exists(gene_rds)){
        genes_txdb <- readRDS(gene_rds)
      }else{
        genes_txdb <- GenomicFeatures::genes(GenomicFeatures::makeTxDbFromEnsembl(organism="Heterocephalus glaber male",release=112))
      }
      seqlevelsStyle(genes_txdb) = "UCSC"
    } else{
      require(sprintf("TxDb.%s.UCSC.%s.knownGene",SPECIES,SPECIES_BUILD),character.only = TRUE)
      genes_txdb <- genes(get(sprintf("TxDb.%s.UCSC.%s.knownGene",SPECIES,SPECIES_BUILD)))
    }
}

if(system("uname",intern = TRUE)=="Darwin"){
  MC_CORES=4
  LUSTRE="/Users/nw14/volumes/nw14_network/LUSTRE"
  NFS="/Users/nw14/volumes/nw14_network/NFS"
  TMPDIR="~/tmp"

}else{
  MC_CORES=4
  LUSTRE="/lustre"
  NFS="/nfs"
  TMPDIR=sprintf("%s/scratch119/casm/team273jn/nw14/tmp",LUSTRE)
}


get_loh=function(sample,minlen=2e6,maxgap=1e6,project="1805"){
  cn=get_cn(cn_summary_file = sprintf("/lustre/scratch124/casm/references/nst_links/live/%s/%s/%s.ascat_ngs.summary.csv",project,sample,sample))
  #idx=which(cn$minor==0 & cn$end-cn$start>minlen & cn$chr %in% 1:22)
  if(is.null(cn)){
    return(NULL)
  }
  idx=which(cn$minor==0 & cn$end-cn$start>minlen )
  cn$sample=sample
  cn=cn[idx,]
  merge_cn(cn,maxgap)
}

get_cna=function(sample,minlen=2e6,maxgap=1e6,project="1805"){
  cn=get_cn(cn_summary_file = sprintf("/lustre/scratch124/casm/references/nst_links/live/%s/%s/%s.ascat_ngs.summary.csv",project,sample,sample))
  if(is.null(cn)){
    cat("No copy number file for ",sample,"\n")
    return(NULL)
  }
  ##browser()
  #idx=which(cn$minor!=cn$major & cn$major!=1 & cn$end-cn$start>minlen & cn$chr %in% 1:22)
  idx=which(cn$major>1 & cn$end-cn$start>minlen )
  cn$sample=sample
  cn=cn[idx,]
  merge_cn(cn,maxgap)
}

merge_cn=function(cn,maxgap=1e6){
  if(dim(cn)[1]<2){
    return(cn)
  }
  uchr=sort(unique(cn$chr))
  if(length(uchr)>1){
    do.call("rbind",lapply(uchr,function(chr){cn_tmp=cn[which(cn$chr==chr),];merge_cn(cn_tmp,maxgap)}))
  }else{
    ##Must be a single chr
    cn=cn[order(cn$start),]
    N=dim(cn)[1]
    ##print(cn)
    ###gap=cn$start[2:N]-cn$end[1:(N-1)]
    j=1
    for(i in 1:(N-1)){
      if((cn$start[j+1]-cn$end[j])<maxgap){
        cn$end[j]=cn$end[j+1]
        cn=cn[-(j+1),]
      }else{
        j=j+1
      }
    }
    cn
  }

}

merge_mcna=function(cn,maxgap=1e6){
  if(is.null(cn) || dim(cn)[1]<2){
    return(cn)
  }
  uchr=sort(unique(cn$chr))
  if(length(uchr)>1){
    do.call("rbind",lapply(uchr,function(chr){cn_tmp=cn[which(cn$chr==chr),];merge_mcna(cn_tmp,maxgap)}))
  }else{
    ##
    cn=cn[order(cn$start),]
    cs=cn$start[1]
    
    for(i in 1:nrow(cn)){
      if(cn$start[i]<cs+maxgap){
        cn$start[i]=cs
      }else{
        cs=cn$start[i]
      }
    }
    cn=cn[order(-cn$end),]
    ce=cn$end[1]
    for(i in 1:nrow(cn)){
      if(cn$end[i]>ce-maxgap){
        cn$end[i]=ce
      }else{
        ce=cn$end[i]
      }
    }
    cn[order(cn$start,-cn$end),]
  }
}

report_cna_bedfiles_from_ascat=function(cfg,## e.g.read.table("../config/PD63423.cfg",head=T,sep="\t") %>% filter(TYPE=="C") %>% filter(LABEL!="PD63423a_lo0002")
maxgap=5e5,minlength=2e6){
  ##browser()
  CNA=do.call("rbind",lapply(1:nrow(cfg),
                             function(i){
                               tmp=get_cna(cfg$LABEL[i],project=cfg$PROJECT[i],minlen=minlength)
                               if(!is.null(tmp) && nrow(tmp)==0){
                                 tmp %>% mutate(sample=character(0),SHORT_LABEL=character(0))
                               }else{
                                  if(is.null(tmp)){
                                    NULL
                                  }else{
                                    tmp %>% mutate(sample=cfg$LABEL[i],SHORT_LABEL=cfg$SHORT_LABEL[i])
                                  }
                               }
                               }
                             )
  )
  LOH=do.call("rbind",lapply(1:nrow(cfg),
                             function(i){
                               tmp=get_loh(cfg$LABEL[i],project=cfg$PROJECT[i],minlen=minlength)
                               if(!is.null(tmp) && nrow(tmp)==0){
                                 tmp %>% mutate(sample=character(0),SHORT_LABEL=character(0))
                               }else{
                                 if(is.null(tmp)){
                                   NULL
                                 }else{
                                   tmp %>% mutate(sample=cfg$LABEL[i],SHORT_LABEL=cfg$SHORT_LABEL[i])
                                 }
                               }
                             }
  )
  )
  ##browser()
  list(CNA=merge_mcna(CNA,maxgap = maxgap) %>% group_by(chr,start,end) %>%
         summarise(major_rng=paste(range(major),collapse="-"),minor_rng=paste(range(minor),collapse="-"),
                   major=max(major),minor=max(minor),
                   long_sample=paste(sample,collapse=","),
                   sample=paste(SHORT_LABEL,collapse=",")) %>% mutate(label=sprintf("cna:chr%s:%dMb-%dMb",chr,round(start/1e6),round(end/1e6))
                   ),
       LOH=merge_mcna(LOH,maxgap = maxgap) %>% group_by(chr,start,end) %>%
         summarise(major_rng=paste(range(major),collapse="-"),minor_rng=paste(range(minor),collapse="-"),
                   major=max(major),minor=0,
                   long_sample=paste(sample,collapse=","),
                   sample=paste(SHORT_LABEL,collapse=",")) %>% mutate(label=sprintf("loh:chr%s:%dMb-%dMb",chr,round(start/1e6),round(end/1e6))
                   ),
       LOH_ALL=LOH,
       CNA_ALL=CNA
  )
}

plot_all_cn=function(patient,ignore_samples=c(),b.plot.to.file=TRUE){
  xx=read.table(sprintf("../config/%s.cfg",patient),head=T,stringsAsFactors = FALSE)
  
  xx=xx[which(xx$TYPE=="C"),]
  xx=xx[order(xx$LABEL),]
  samples=xx$LABEL
  project=xx$PROJECT
  idxkeep=which(!(samples %in% ignore_samples))
  samples=samples[idxkeep]
  project=project[idxkeep]
  projectlabel=paste(sort(unique(project)),collapse="_")
  M=ceiling(length(samples)/2)
  if(b.plot.to.file){
    pngfile=sprintf("../post_qc/%s_project%s_cn.png",patient,projectlabel)
    png(pngfile,h=M*180,w=360*3*2)
  #
    par(mfcol=c(M,2))
  }else{
    par(mfcol=c(1,2))
  }
  for(i in 1:length(samples)){
    sample=samples[i]
    summaryf=sprintf("%s/cancer_ref01/nst_links/live/%s/%s/%s.ascat_ngs.summary.csv",NFS,project[i],sample,sample)
    if(!file.exists(summaryf)){
      cat("can't find cn for",sample,"\n")
    }else{
    plot_cn(get_cn(cn_summary_file = summaryf));
    title(sample,cex=3)
    }
  }
  if(b.plot.to.file){
  dev.off()
  }
}


report_cn=function(patient,chr,ignore_samples=c(),do.baf=TRUE){
  xx=read.table(sprintf("../config/%s.cfg",patient),head=T,stringsAsFactors = FALSE)
  xx=xx[which(xx$TYPE=="C"),]
  xx=xx[order(xx$LABEL),]
  samples=xx$LABEL
  project=xx$PROJECT
  idxkeep=which(!(samples %in% ignore_samples))
  samples=samples[idxkeep]
  project=project[idxkeep]
  projectlabel=paste(sort(unique(project)),collapse="_")
  #


  cn=do.call("rbind",lapply(1:length(samples),function(i) get_loh(samples[i],maxgap=2e6,project = project[i])))
  cna=do.call("rbind",lapply(1:length(samples),function(i) get_cna(samples[i],maxgap=2e6,project = project[i])))

  if(do.baf){
    bafs=getbafdat(samples,chr,mc.core=8,project=project)
    sn=sapply(1:length(bafs),function(i) bafs[[i]]$sample)
    #s1=sort(unique(cn$sample[which(cn$chr==chr)]))
    ##sn=sn[which(sn %in% s1)]
    ##idx=match(s1,sn)
    idx=1:length(sn)
    M=ceiling(length(idx)/2)
    ##system("mkdir ../post_qc/%s")

    pngfile=sprintf("../post_qc/%s_chr%s_project%s_bafs.png",patient,chr,projectlabel)
    cat("plotting too",pngfile,"\n")
    png(pngfile,h=M*120,w=360*2)
    #
    par(mfcol=c(M,2))
    for(i in idx){
      if(!is.null(bafs[[i]]$df)){
        baf=bafs[[i]]$df;
        sample=bafs[[i]]$sample;
        plot(y=baf$BAF,x=baf$Position/1e6,cex=0.5,xlab="Position (Mb)",ylab="BAF",pch=19)
        #
        title(sample)
      }
    }
    dev.off()


    list(cn=cn,bafs=bafs)
  }else{
    list(cn=cn,cna=cna)
  }
}

get_cnfile=function(PD,colony,type="ascat"){
  if(!type %in% c("ascat","baf")){
    stop("get_cnfile:unsupported type")
  }
  idx=which(PD$inf$cfg$SHORT_LABEL==colony)
  sample=PD$inf$cfg$LABEL[idx]
  project=PD$inf$cfg$PROJECT[idx]
  if(type=="ascat"){
    suffix="ascat_ngs.summary.csv"
  }else{
    ##baf
    suffix="ascat_ngs.cn.tsv.gz"
  }
  target=sprintf("%s/%s_%s_%s",CACHE,project,sample,suffix)
  if(!file.exists(target)){
    src=sprintf("%s/cancer_ref01/nst_links/live/%s/%s/%s.%s",NFS,project,sample,sample,suffix)
    status=file.copy(src,target)
    if(!status){
      cat(sprintf("local copy: %s -> %s",src,target))
      stop("unable to copy file..")
    }
  }
  target
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


liftover2=function(snp,chr,position,hgflag="hg18ToHg19"){
  out=data.frame(chr=sprintf("chr%s",chr),
                 pos1=sprintf("%d",position),pos2=sprintf("%d",position+1) ,
                 snp=snp,stringsAsFactors=FALSE)
  tmp1=tempfile(tmpdir="/lustre/scratch119/casm/team273jn/nw14/tmp")
  tmp2=tempfile(tmpdir="/lustre/scratch119/casm/team273jn/nw14/tmp")
  write.table(out,tmp1,row.names=FALSE,col.names=FALSE,quote=FALSE)
  cmd=sprintf("/nfs/casm/team273jn/nw14/software/liftOver %s /nfs/casm/team273jn/nw14/software/liftover.chains/%s.over.chain.gz %s %s.error",
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

plot_cn=function(dat,label="",xlab="",yylim=c(-1,4)){
  if(is.null(dat)){
    warning("Empty input to plot_cn")
    return(NULL)
  }
  inf=get_chr_info2()
  p.start=inf$start_pos[match(dat$chr,inf$chr)]+dat$start
  p.end=inf$start_pos[match(dat$chr,inf$chr)]+dat$end
  ###browser()
  plot(NULL,
       pch=19,cex=0.5,xaxt="n",
       ylim=yylim,ylab="",
       yaxt="n",xlab=xlab,xlim=c(0,p.end[length(p.end)]),main=label)
  abline(h=0:10,col="grey")
  axis(side = 2,at = 0:4,las=1)
  rect(xleft = p.start,xright=p.end,
       ybottom = dat$minor-0.1,ytop=dat$minor-0.01,col="green",border="green")

  rect(xleft = p.start,xright=p.end,
       ybottom = dat$major+0.01,ytop=dat$major+0.1,col="red",border = "red")
  abline(v=inf$start_pos)
  mtext(at=inf$start_pos+0.5*inf$end,side = 1,cex=0.8,line = 1,text = inf$chr)

}

##Would be much better if the following just looked at the genome.fa.fai file rather than hard-coded chromosome lengths.
get_chr_info2=function(){
  chr_info=data.frame(chr=c(sprintf("%d",1:22),c("X","Y")) ,end=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,
                                                                  146364022,141213431,135534747,134996516,133851895,115169878,107349540,
                                                                  102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,
                                                                  155270560,59373566))
  chr_info$start_pos=c(0,cumsum(chr_info$end)[-length(chr_info$end)])
  chr_info
}

get_cn=function(cn_summary_file){
  ##cat("opening ",cn_summary_file,"\n")
  if(!file.exists(cn_summary_file)){
    warning(sprintf("%s: does not exist",cn_summary_file))
    return(NULL)
  }
  cn=read.csv(cn_summary_file,header = FALSE)
  cn$start=cn$V3
  cn$end=cn$V4
  cn$chr=gsub("chr","",cn$V2)
  ###V7=total copy number
  ## V8=minor allele copy number
  cn$major=cn$V7-cn$V8
  cn$minor=cn$V8
  cn[,-grep("^V",colnames(cn))]
}
require("readr")
getbaf=function(x,chr,project="1805"){

  baf=read_tsv(sprintf("%s/cancer_ref01/nst_links/live/%s/%s/%s.ascat_ngs.cn.tsv.gz",NFS,project,x,x));
  baf[which(baf$Chromosome==chr),]
}



getbafdat=function(samples,chr,project,mc.cores=4){
  mclapply(1:length(samples),function(i){
    key=sprintf("%s/%s_%s_chr%s.RDS",TMPDIR,samples[i],project[i],chr)
    if(file.exists(key)){##cat("Using cached",key,"\n")
      return(readRDS(key))}
    baf=try(getbaf(samples[i],chr,project=project[i]),silent=TRUE)
    if(class(baf)=="try-error"){
      baf=list(sample=samples[i],df=NULL)}
    else{
      baf=list(sample=samples[i],df=baf[,c("Position","BAF")])}
    ;
    saveRDS(baf,key)
    baf
  },mc.cores = mc.cores)

}


get_hethom=function(PD,refcolony,b_do_plot=FALSE){
  loh=PD$pdx$meta$LOH[which(sapply(PD$pdx$meta$LOH,function(x) x$ploidy)==2)]
  if(length(loh)==0){
    return(NULL)
  }
  baff=get_cnfile(PD,refcolony,type="baf")
  bafref=read_tsv(baff)
  bafrefhet=bafref[!is.na(bafref$BAF) & bafref$BAF>0.2 & bafref$BAF<0.8,]
  ##Here we are only interested in copy number neutral LOH
  bins=readRDS(MUTCOUNTBIN)

  hethomres=lapply(loh,function(x){

    label=sprintf("%s (%s:%s-%s)",x$LABEL,x$chr,x$start,x$end)
    ##representative sample
    colony=x$samples[1]

    #kable(hethom$inf, "html") %>% kable_styling("striped")
    ascatfile=get_cnfile(PD,colony,type="ascat")
    baff=get_cnfile(PD,colony,type="baf")
    baf=read_tsv(baff);
    baf=baf[which(!is.na(baf$BAF) & baf$Chromosome==x$chr),]
    bafhet=bafrefhet[which(bafrefhet$Chromosome==x$chr),]
    baf=baf[baf$Position %in% bafhet$Position,]
    mbaf=frollmean(2*abs(baf$BAF-0.5),20,align = "right")
    mbaf2=frollmean(2*abs(baf$BAF-0.5),20,align = "left")
    idx=which(mbaf>0.8)
    Pos1=baf$Position[min(which(mbaf2>0.8 & baf$Position>=x$start))]
    Pos2=baf$Position[max(which(mbaf>0.8 & baf$Position<=x$end))]
    cat("Pos1=",Pos1,",Pos2=",Pos2,"\n")
    x$start=Pos1
    x$end=Pos2
    ###browser()
    hethom=get_loh_tree(PD,x$node,x)
    hethom$label=x$LABEL
    hethom$node=x$node
    hethom$het.sensitivity=PD$pdx$tree_ml$sensitivity.snv[match(x$node,PD$pdx$tree_ml$edge[,2])]
    hethom$chr=x$chr
    hethom$start=x$start
    hethom$end=x$end

    if(b_do_plot){
      kable(hethom$inf, "html") %>% kable_styling("striped")
      plot_cn(get_cn(cn_summary_file = ascatfile),label = sprintf("%s: Representative of %d colonies : %s",colony,length(x$samples),label))

      plot(baf$Position/1e6,baf$BAF,pch=19,cex=0.5,xlab="Position(Mb)",ylab="BAF");title(sprintf("%s: BAF: Representative of %d colonies : %s",colony,length(x$samples),label))
      lines(baf$Position/1e6,mbaf,col="magenta")
      abline(v=hethom$inf$Pos[hethom$inf$Chrom==x$chr & hethom$inf$het]/1e6,col="green")
      abline(v=hethom$inf$Pos[hethom$inf$Chrom==x$chr & hethom$inf$hom]/1e6,col="blue")
      abline(v=c(x$start/1e6,x$end/1e6),col="red")
      abline(v=c(Pos1/1e6,Pos2/1e6),col="red",lty="dashed")
      legend("bottomright",c("HET","HOM","LOH Boundary","Rolling Homozygosity (20 SNPs)"),col=c("green","blue","red","magenta"),lwd=2,bg = "white",cex = 0.8)
    }




    idx=which(bins$chr==x$chr & bins$start>=(1e5*floor(x$start/1e5)) & bins$end<=(1e5*floor(x$end/1e5)))
    hethom$kb=length(idx)*1e5
    hethom$count_in_bin=sum(bins$count[idx])
    hethom$count_se=sqrt(sum(bins$count[idx]))
    hethom$pmut=sum(bins$prob[idx])
    ##Estimate of uncertainty in local density
    hethom$pmut_se=hethom$pmut/sqrt(hethom$count_in_bin)
    ##Use the null model estimates
    ut=PD$fit$poisson_tree$nullmodel$ultratree
    T=ut$edge.length[which(ut$edge[,2]==x$node)]
    ##Use null model lambda
    lambda=20
    lambda=PD$fit$poisson_tree$nullmodel$lambda$mean
    dat=list(N=hethom$nhom,M=hethom$nhet,phom=0.5,phet=hethom$het.sensitivity,L=lambda*T*hethom$pmut,Lstd=lambda*T*hethom$pmut_se)
    ##Fit stan model...
    #browser()
    stanr=stan("stanmodels/cntime.stan", data=dat,chains = 3,iter=20000,cores=3)
    sumc=rstan::summary(stanr,probs = c(0.025,  0.50, 0.975))$summary
    for( z in colnames(sumc)){
      hethom[[sprintf("x%s",z)]]=sumc["x",z]
    }
    for( z in colnames(sumc)){
      hethom[[sprintf("l%s",z)]]=sumc["l",z]
    }
    hethom
  })
  othernames=c("inf")
  keepnames=setdiff(names(hethomres[[1]]),othernames)
  list(hethom=do.call("rbind",lapply(hethomres,function(x) as.data.frame(x[keepnames]))),
       other=do.call("rbind",lapply(hethomres,function(x) as.data.frame(x[othernames]))))
}

##A bit of copy and paste from get_hethom
get_cnainf=function(PD,refcolony,b_do_plot=FALSE){
  if(length(PD$pdx$meta$CNA)==0){
    return(NULL)
  }
  bins=readRDS(MUTCOUNTBIN)
  baff=get_cnfile(PD,refcolony,type="baf")
  bafref=read_tsv(baff)
  bafrefhet=bafref[!is.na(bafref$BAF) & bafref$BAF>0.2 & bafref$BAF<0.8,]

  ##Here we are only interested in copy number neutral LOH
  idxx=grep("\\+",sapply(PD$pdx$meta$CNA,function(x) x$LABEL))
  if(length(idxx)==0){
    return(NULL)
  }
  hethomres=lapply(PD$pdx$meta$CNA[idxx],function(x){
    label=sprintf("%s (%s:%s-%s)",x$LABEL,x$chr,x$start,x$end)
    ##representative sample
    colony=x$samples[1]

    #kable(hethom$inf, "html") %>% kable_styling("striped")
    ascatfile=get_cnfile(PD,colony,type="ascat")
    baff=get_cnfile(PD,colony,type="baf")
    baf=read_tsv(baff);
    baf=baf[which(!is.na(baf$BAF) & baf$Chromosome==x$chr),]
    bafhet=bafrefhet[which(bafrefhet$Chromosome==x$chr),]
    baf=baf[baf$Position %in% bafhet$Position,]
    mbaf=frollmean(2*abs(baf$BAF-0.5),20,align = "right")
    mbaf2=frollmean(2*abs(baf$BAF-0.5),20,align = "left")
    idx=which(mbaf>0.8)

    hethom=get_duplication_dat(PD,x$node,x)
    if(is.null(hethom)){
      return(NULL)
    }
    hethom$label=x$LABEL
    hethom$node=x$node
    hethom$het.sensitivity=PD$pdx$tree_ml$sensitivity.snv[match(x$node,PD$pdx$tree_ml$edge[,2])]
    hethom$chr=x$chr
    hethom$start=x$start
    hethom$end=x$end
    if(b_do_plot){
      kable(hethom$inf, "html") %>% kable_styling("striped")
      plot_cn(get_cn(cn_summary_file = ascatfile),label = sprintf("%s: Representative of %d colonies : %s",colony,length(x$samples),label))
      xx=cbind(baf$Position/1e6,baf$BAF)
      smoothScatter(xx,nrpoints = 0,xlab="Position(Mb)",ylab="BAF")
      abline(v=hethom$inf$Pos[hethom$inf$Chrom==x$chr & hethom$inf$is_dup]/1e6,col="green")
      abline(v=hethom$inf$Pos[hethom$inf$Chrom==x$chr & !hethom$inf$is_dup]/1e6,col="blue")
      abline(v=c(x$start/1e6,x$end/1e6),col="red")
      #abline(v=c(Pos1/1e6,Pos2/1e6),col="red",lty="dashed")
      legend("bottomright",c("VAF=2/3","VAF=1/3","CNA Boundary"),col=c("green","blue","red"),lwd=2,bg = "white",cex = 0.8)
    }




    idx=which(bins$chr==x$chr & bins$start>=(1e5*floor(x$start/1e5)) & bins$end<=(1e5*floor(x$end/1e5)))
    hethom$kb=length(idx)*1e5
    hethom$count_in_bin=sum(bins$count[idx])
    hethom$count_se=sqrt(sum(bins$count[idx]))
    hethom$pmut=sum(bins$prob[idx])
    ##Estimate of uncertainty in local density
    hethom$pmut_se=hethom$pmut/sqrt(hethom$count_in_bin)
    ##Use a global estimate of 20..
    ##Use the null model estimates
    ut=PD$fit$poisson_tree$nullmodel$ultratree
    T=ut$edge.length[which(ut$edge[,2]==x$node)]
    ##Use null model lambda
    lambda=PD$fit$poisson_tree$nullmodel$lambda$mean
    #lambda=20
    dat=list(N=hethom$dupcount,M=hethom$ndupcount,pdup=hethom$het.sensitivity,pndup=hethom$het.sensitivity,L=lambda*T*hethom$pmut,Lstd=lambda*T*hethom$pmut_se)
    ##Fit stan model...
    stanr=stan("stanmodels/cntimedup.stan", data=dat,chains = 3,iter=20000,cores=3)
    sumc=rstan::summary(stanr,probs=c(0.025,0.5,0.975))$summary
    for( z in colnames(sumc)){
      hethom[[sprintf("x%s",z)]]=sumc["x",z]
    }
    for( z in colnames(sumc)){
      hethom[[sprintf("l%s",z)]]=sumc["l",z]
    }
    hethom
  })
  othernames=c("inf")
  keepnames=setdiff(names(hethomres[[1]]),othernames)
  list(hethom=do.call("rbind",lapply(hethomres,function(x) as.data.frame(x[keepnames]))),
       other=do.call("rbind",lapply(hethomres,function(x) x[othernames])))
}

get_duplication_dat=function(PD,node,loh,exclude=c(),hom.threshold=0.8){
  inf=PD$pdx$dat
  mtr=inf$mtr
  dep=inf$dep
  cs=PD$pdx$meta$clones_short
  kids=c(loh$node,get_all_node_children(loh$node,PD$pdx$tree_ml))
  tips=kids[kids<=length(PD$pdx$tree_ml$tip.label)]
  mt=PD$pdx$tree_ml$tip.label[tips]
  wt=setdiff(cs,mt)
  ##check
  wt=setdiff(wt,exclude)
  mt=setdiff(mt,exclude)
  if(length(union(mt,intersect(loh$samples,cs)))!=length(mt)){
    cat("mt:",mt,"\n")
    cat("loh$samples",setdiff(loh$samples,cs),"\n")
    browser()
    stop("Inconsistency between LOH samples and associated node tips")
  }
  idx=with(inf$details,which(Chrom==loh$chr & Pos>loh$start & Pos<loh$end & TYPE=="SNV" & BGLOD==0))
  idx2=which(inf$details$node[idx]==node)
  cat("Number of muts overlapping loh shared in tree =",length(idx2),"\n")
  #cat("WT counts:\n")
  #print(inf$mtr[idx[idx2],wt])
  tmp=inf$details[idx[idx2],c("Chrom","Pos")]
  #print(inf$details[idx[idx2],1:10])
  ##browser()
  cat("WT:",wt,"\n")
  cat("MT:",mt,"\n")

  if(length(mt)==1){
    MT=inf$mtr[idx[idx2],mt]
    DP=inf$dep[idx[idx2],mt]
    OFS=inf$ofs[idx[idx2],mt]==0
  }else{
    if(length(idx[idx2])<2){
      return(NULL)
    }
    MT=rowSums(inf$mtr[idx[idx2],mt])
    DP=rowSums(inf$dep[idx[idx2],mt])
    OFS=rowSums(inf$ofs[idx[idx2],mt]==0)>0
  }
  if(length(wt)==1){
    MTW=inf$mtr[idx[idx2],wt]
    DPW=inf$dep[idx[idx2],wt]
  }else{
    MTW=rowSums(inf$mtr[idx[idx2],wt])
    DPW=rowSums(inf$dep[idx[idx2],wt])
  }
  tmp$MT=MT
  tmp$DP=DP
  tmp$VAF=MT/DP
  tmp$MTW=MTW
  tmp$DPW=DPW
  tmp$VAFW=MTW/DPW
  tmp$OFS=OFS
  #ldup=dbinom(tmp$MT,tmp$DP,prob=2/3)
  #lndup=dbinom(tmp$MT,tmp$DP,prob=1/3)
  #tmp$p_dup=ldup/(ldup+lndup)

  status=classify_2state(tmp$MT,tmp$DP,probs = c(1/3,2/3))
  tmp$is_dup=ifelse(status==2,1,ifelse(status==1,0,NA))
  list(dupcount=sum(status==2,na.rm=TRUE),ndupcount=sum(status==1,na.rm=TRUE),inf=tmp)
}







get_loh_tree=function(PD,node,loh,exclude=c(),hom.threshold=0.8){
  inf=PD$pdx$dat
  mtr=inf$mtr
  dep=inf$dep
  cs=PD$pdx$meta$clones_short
  kids=c(node,get_all_node_children(node,PD$pdx$tree_ml))
  tips=kids[kids<=length(PD$pdx$tree_ml$tip.label)]
  mt=PD$pdx$tree_ml$tip.label[tips]
  wt=setdiff(cs,mt)
  ##check
  wt=setdiff(wt,exclude)
  mt=setdiff(mt,exclude)
  if(length(union(mt,intersect(loh$samples,cs)))!=length(mt)){
    cat("mt:",mt,"\n")
    cat("loh$samples",setdiff(loh$samples,cs),"\n")
    browser()
    stop("Inconsistency between LOH samples and associated node tips")
  }
  idx=with(inf$details,which(Chrom==loh$chr & Pos>loh$start & Pos<loh$end & TYPE=="SNV" & BGLOD==0))
  idx2=which(inf$details$node[idx]==node)
  cat("Number of muts overlapping loh shared in tree =",length(idx2),"\n")
  #cat("WT counts:\n")
  #print(inf$mtr[idx[idx2],wt])
  tmp=inf$details[idx[idx2],c("Chrom","Pos")]
  #print(inf$details[idx[idx2],1:10])
  ##browser()


  if(length(mt)==1){
    MT=inf$mtr[idx[idx2],mt]
    DP=inf$dep[idx[idx2],mt]
  }else{
    #browser()
    MT=rowSums(inf$mtr[idx[idx2],mt])
    DP=rowSums(inf$dep[idx[idx2],mt])
  }
  if(length(wt)==1){
    MTW=inf$mtr[idx[idx2],wt]
    DPW=inf$dep[idx[idx2],wt]
  }else{
    MTW=rowSums(inf$mtr[idx[idx2],wt])
    DPW=rowSums(inf$dep[idx[idx2],wt])
  }
  tmp$MT=MT
  tmp$DP=DP
  tmp$VAF=MT/DP
  tmp$MTW=MTW
  tmp$DPW=DPW
  tmp$VAFW=MTW/DPW

  tmp$het_p=pbinom(MT-1,DP,prob=0.5,lower.tail = FALSE)
  tmp$het=(MT/DP)<hom.threshold & tmp$het_p >0.05
  tmp$hom=(MT/DP)>hom.threshold & tmp$het_p <0.05
  status=classify_2state(MT,DP,probs = c(0.5,1))
  nheto=sum(tmp$het)
  nhomo=sum(tmp$hom)
  nhet=sum(status==1,na.rm = TRUE)
  nhom=sum(status==2,na.rm = TRUE)
  list(nhet=nhet,nhom=nhom,inf=tmp)
  #vaf_tree=with(inf,rowSums(mtr[idx[idx2],mt])/rowSums(dep[idx[idx2],mt]))
  #cat("WT counts")
  #print(inf$mtr[idx[idx2],mt])
  #print(inf$mtr[idx[idx2],wt])
  #print(inf$details[idx[idx2],1:10])
  #cat("MT VAFs\n")
  #idx.tree=idx[idx2]
  #print(inf$mtr[idx[idx2],mt]/inf$dep[idx[idx2],mt])
  #idx.tree
}

classify_2state=function(mtr,depth,probs,epsilon=0.01,threshold=0.95){
  if(length(probs)!=2){
    stop("Need to supply to probs")
  }
  #Maths for bernoulli trial.
  #p(mut)=p(mut|mutant read)p(mutant read)+p(mut|ref read)p(ref read)
  #p(mut)=(1-epsilon)*mu+(epsilon/3)(1-mu)
  probs=probs*(1-epsilon)+(epsilon/3)*(1-probs)
  l1=dbinom(mtr,depth,prob=probs[1] )
  l2=dbinom(mtr,depth,prob=probs[2])
  l1n=l1/(l1+l2)
  ifelse(l1n>threshold,1,ifelse(1-l1n>threshold,2,NA))
}




rawviolinplot=function(x,y,ylim,N=50){
  br=seq(ylim[1],ylim[2],length.out = N)
  vals=cut(y,breaks=br,label=FALSE)
  counts=sapply(1:(length(br)-1),function(x) length(which(vals==x)))
  mc=max(counts)
  for(i in 1:(length(br)-1)){
    dx=0.4*counts[i]/mc
    rect(xleft = x-dx,xright = x+dx,ytop = br[i+1],ybottom = br[i],col = "grey",border = "black")
  }
}










ggplot_color=function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

getcolor_df=function(mut_count){
  col.df=data.frame(patient=unique(as.character(mut_count$patient)),stringsAsFactors = FALSE)
  ##col.df$median_depth=with(mut_count,sapply(col.df$patient,function(x){median(depth[which(patient==x)])}))
  ##return(col.df)

  col.df$col=ggplot_color(dim(col.df)[1])#brewer.pal(12,"Paired")[c(1:10,12,11)][1:dim(col.df)[1]]
  col.df
}




get_grange=function(chr,pos,ref,alt){
  idx=which(nchar(ref)+nchar(alt)==2)
  df=data.frame(chr=chr[idx],pos=pos[idx],ref=ref[idx],alt=alt[idx],stringsAsFactors = FALSE)
  gr=GRanges(sprintf("chr%s",df$chr), IRanges(start=as.numeric(df$pos), width=1))
  #values(gr)=cbind(gr,REF=df$ref,ALT=df$alt)  ##Stopped working with GenomicRanges_1.38.0 or GenomicRanges_1.37.0
  mcols(gr)=data.frame(REF=ref,ALT=alt,stringsAsFactors = FALSE)
  gr
}

get_mut_matrix=function(chr,pos,ref,alt){
  gr=GRanges(sprintf("chr%s",chr), IRanges(start=as.numeric(pos), width=1))
  GenomeInfoDb::genome(gr) = SPECIES_BUILD
  #values(gr)=cbind(gr,REF=ref,ALT=alt)
  mcols(gr)=data.frame(REF=ref,ALT=alt,stringsAsFactors = FALSE)
  mut_matrix(list(gr=gr),ref_genome)
}
add_trin=function(dat){
  dato=dat
  bases=c("A","C","G","T")
  idx.keep=with(dat,which(Ref %in% bases & Alt %in% bases)) ##Keeps SNV
  if(length(idx.keep)==0){
    stop("Unexpected have no SNVs to get context for...")
  }
  dat=dat[idx.keep,]
  ##browser()
  gg=with(dat,get_grange(chr=Chrom,pos = Pos,ref=Ref,alt = Alt))
  GenomeInfoDb::genome(gg) = SPECIES_BUILD
  #ctx=mut_context(gg,ref_genome =ref_genome)
  tc=type_context(gg,ref_genome =ref_genome)
  dat$trin=sprintf("%s[%s]%s",substr(tc$context,1,1),tc$types,substr(tc$context,3,3))
  trans=as.character(mut_strand(gg,genes_txdb))
  trans=c("T","U","-")[match(trans,c("transcribed","untranscribed","-"))]
  dato$trin="INDEL"
  dato$transcription=NA
  dato$trin[idx.keep]=dat$trin
  dato$transcription[idx.keep]=trans
  dato
}


get_context_counts_stranded=function(details){

  gr=GRanges(sprintf("chr%s",details$Chrom), IRanges(start=as.numeric(details$Pos), width=1))
  GenomeInfoDb::genome(gr) = SPECIES_BUILD
  #values(gr)=cbind(gr,REF=details$Ref,ALT=details$Alt)
  mcols(gr)=data.frame(REF=ref,ALT=alt,stringsAsFactors = FALSE)
  ##strand = mut_strand(vcfs[[1]], genes_hg19)
  mut_matrix_stranded(list(gr=gr),ref_genome, genes_txdb)
  ##mut_matrix(list(gr=gr),ref_genome)
}

get_mut_matrix_stranded=function(chr,pos,ref,alt){
  gr=GRanges(sprintf("chr%s",chr), IRanges(start=as.numeric(pos), width=1))
  GenomeInfoDb::genome(gr) = SPECIES_BUILD
  #values(gr)=cbind(gr,REF=ref,ALT=alt)
  mcols(gr)=data.frame(REF=ref,ALT=alt,stringsAsFactors = FALSE)
  ##strand = mut_strand(vcfs[[1]], genes_hg19)
  mut_matrix_stranded(list(gr=gr),ref_genome, genes_txdb)
  ##mut_matrix(list(gr=gr),ref_genome)
}

get_mut_matrix_stranded_all=function(chr,pos,ref,alt){
  tmp=data.frame(Chrom=chr,Pos=pos,Ref=ref,Alt=alt,stringsAsFactors=FALSE)
  if(dim(tmp)[1]==0){
    sbs288=c()
  }else{
    tmp=add_trin(tmp)
    sbs288=sprintf("%s:%s",tmp$trin,tmp$transcription)
  }
  alpha288=rep(MutationalPatterns:::TRIPLETS_96,each=3)
  alpha288=sprintf("%s:%s",alpha288,c("T","U","-"))
  #counts=rep(0,288);
  #counts=tabulate[match(sbs288,alph2)]
  counts=sapply(alpha288,function(x) length(which(sbs288==x)))
  out=matrix(counts,ncol=1)
  rownames(out)=alpha288
  out
}



##sharedness metric..
get_smetric=function(tree){
  n=sapply(tree$edge[,2],function(x) length(get_samples_in_clade(x,tree)))/length(tree$tip.label)
  t2=tree
  t2$edge.length=t2$edge.length*n
  idx=match(1:length(t2$tip.label),t2$edge[,2])
  s=nodeHeights(t2)[idx,2]/nodeHeights(tree)[idx,2]
  #browser()
  s
}

plot_all_adjusted_trees=function(PD){
  par(mfcol=c(2,2))
  for(field in c("germline.multi","reg","hybrid","hybrid.multi")){
    PD=apply_adjustment(PD,field)
    PD$pdx=set_color_by_age(PD$pdx)
    tree=PD$pdx$tree_ml
    nh=nodeHeights(tree)
    ht=nh[match(1:length(tree$tip.label),tree$edge[,2]),2]
    label=sprintf("%s:%s:sd/mean=%5.4f",PD$patient,field,sd(ht[ht>1])/mean(ht[ht>1]))
    tree=plot_basic_tree(PD$pdx,label,cex.label=0,cex.terminal=0.5,left.margin.prop =0.15)
    legend("left",legend = floor(PD$pdx$age.df$age),col=PD$pdx$age.df$color,pch=19,title = "Age at Sample")
  }
}

get_driver_scheme_MPN=function(){
  driver.scheme=read.table("driver_scheme.txt",head=T,stringsAsFactors = FALSE,sep="\t")
  driver.group=read.table("driver_groups.txt",head=T,stringsAsFactors = FALSE,sep="\t",comment.char = "")
  driver.scheme=driver.scheme[,c("driver","number")] %>% inner_join(driver.group,by="number")
  driver.scheme
}

get_driver_scheme=function(){
  driver.scheme=read.table("driver_scheme_simple.txt",head=T,stringsAsFactors = FALSE,sep="\t")
  n=max(driver.scheme$number)
  pallete=c(RColorBrewer::brewer.pal(9,"Set1")[-6],RColorBrewer::brewer.pal(8,"Dark2"))
  driver.scheme$colour=pallete[driver.scheme$number]
  driver.scheme
}



add_child_count=function(PD){
  nodes=PD$nodes
  nodes$child_count=sapply(nodes$node,function(node) length(get_samples_in_clade(node=node,tree=PD$pdx$tree_ml)))
  PD$nodes=nodes
  PD
}




get_cached_result=function(context,hashme){
  hash=digest(hashme)
  path=sprintf("../cache/%s_%s.Rds",context,hash)
  if(file.exists(path)){
    #cat("\nretrieving cached result",path,"\n")
    return(readRDS(path))
  }else{
    return(NULL)
  }
}

cache_result=function(context,hashme,res){
  cat("Caching results..")
  hash=digest(hashme)
  path=sprintf("../cache/%s_%s.Rds",context,hash)
  cat("to",path,"\n")
  saveRDS(res,path)
}

wraptreefit=function(PD,niter=20000,b.fit.null=FALSE,method=c("poisson_tree"),cores=4,b_use_cache=TRUE,stan_control=list(adapt_delta=0.95),xcross=-1){
  res=refit(PD,niter=niter,b.fit.null=b.fit.null,method=method,cores=cores,b_use_cache=b_use_cache,stan_control=stan_control,xcross=xcross)
  if(is.null(PD$fit)){
    PD$fit=list()
    PD$fit[[method]]=list()
  }
  if(b.fit.null){
    PD$fit[[method]]$nullmodel=res
  }else{
    PD$fit[[method]]$altmodel=res
  }
  PD
}

refit=function(PD,niter=20000,b.fit.null=FALSE,method=c("poisson_tree"),cores=4,b_use_cache=TRUE,stan_control=list(adapt_delta=0.95),xcross=-1){
  tree=get_tree_for_fit(PD)
  if(b.fit.null){
    switch_nodes=c()
    split_nodes=PD$nodes$node
    xcross=c()#rep(xcross,length(switch_nodes)
  }else{
    switch_nodes=PD$nodes$node[which(PD$nodes$status>=0)]
    split_nodes=PD$nodes$node
    xcross=rep(xcross,length(switch_nodes))
  }
  hashme=list(xcross,tree$sensitivity,tree$edge.length,tree$edge,tree$tip.label,niter,cores,b.fit.null,method,PD$nodes,PD$pdx$agedf[,c("tip.label","age_at_sample_exact")])
  context=paste0("fit_tree_",PD$patient)
  res=get_cached_result(context,hashme)
  if(!is.null(res) && b_use_cache){
    cat("using cached result\n")
    return(res)
  }
  ptree=fit_tree(tree=tree,switch_nodes = switch_nodes,xcross = xcross,niter = niter,cores=cores,model  = method,split_nodes = split_nodes,stan_control = stan_control)
  fields=c("ultratree","lambda","split_nodes","upper_node_lims","lower_node_lims","fullres")
  res=ptree[fields]
  res$summary=with(res,
                   cbind(PD$nodes,chknode=split_nodes,
                         as.data.frame(lower_node_lims) %>% (function(x){colnames(x)=sprintf("lower_%s",colnames(x));x}),
                         as.data.frame(upper_node_lims) %>% (function(x){colnames(x)=sprintf("upper_%s",colnames(x));x})
                   ))
  res$param=list(niter=niter,switch_nodes=switch_nodes,split_nodes=split_nodes,tree=tree)
  cache_result(context,hashme,res)
  res
}

cache_PD=function(PD){
  cat("caching result..\n")
  saveRDS(PD,sprintf("%s/%s.RDS",CACHE,PD$patient))

}

get_tree_for_fit=function(PD,type="local",suffix="",atype="hybrid.multi"){
  thistree=PD$pdx$tree_ml
  thistree$agedf=PD$pdx$agedf
  if(type %in% c("global","local")){
    thistree$edge.length=thistree[[sprintf("el.snv.%s.filtered%s",type,suffix)]]
    thistree$sensitivity=thistree[[sprintf("per.branch.sensitivity.%s",atype)]]
  }else{
    stop("unsupported")
    thistree$edge.length=thistree[[sprintf("el.snv%s",suffix)]]
    thistree$sensitivity=thistree[[sprintf("sensitivity.snv%s",suffix)]]
  }
  thistree$sensitivity=ifelse(thistree$sensitivity>0.99,0.99,thistree$sensitivity)
  thistree$agedf$age=thistree$agedf$age_at_sample_pcy
  thistree
}



get_stan_result_df=function(PD,model="poisson_tree",bIsNull=FALSE){
  if(!bIsNull){
    fitres=PD$fit[[model]]$altmodel
    tree=get_colored_markup_tree2(PD$fit[[model]]$altmodel$ultratree,PD$nodes[PD$nodes$status>=0,])
    nodes=rbind(data.frame(node=-1,driver="WT",status=1,child_count=-1,driver2="WT",driver3="WT",stringsAsFactors = FALSE),
                PD$nodes %>% dplyr::select( node, driver ,status, child_count, driver2 ,driver3))
    nodes=nodes[which(nodes$status>=0),]
    tipstatus=tree$rate[which(tree$edge[,2]<=length(tree$tip.label))]
    nodes$patient=PD$patient
    nodes$type="local"
    M=dim(nodes)[1]
    nodes$colony_count=sapply(1:M,function(i){n=length(which(tipstatus==i));if(i==1){n-1}else{n}})
  }else{
    fitres=PD$fit[[model]]$nullmodel
    nodes=data.frame(node=-1,driver="Global",status=1,child_count=-1,driver2="Global",driver3="Global",stringsAsFactors = FALSE)
    nodes$patient=PD$patient
    nodes$type="local"
    nodes=data.frame(node=-1,driver="",status=1,child_count=-1,driver2="",driver3="",stringsAsFactors = FALSE)
  }
  for(x in c("mean","sd","lb","ub","median")){
    nodes[[x]]=fitres$lambda[[x]]
  }
  type="local"
  for(x in c("mean","sd","lb","ub","median")){
    if(type!="snv"){
      nodes[[sprintf("%s_rescaled",x)]]=fitres$lambda[[x]]*PD[[sprintf("%sx.correction2",type)]]
      nodes[["correction"]]=PD[[sprintf("%sx.correction2",type)]]
    }else{
      stop("unsupported path")
    }
  }
  post.dist=rstan::extract(fitres$fullres,par="lambda")
  N=dim(post.dist$lambda)[1]
  nodes[["p_lt_wt"]]=sapply(1:dim(nodes)[1],function(i) mean(post.dist$lambda[,i]<post.dist$lambda[,1]))
  nodes[["p_lt_wt"]]=ifelse(nodes[["p_lt_wt"]]<1/N,1/N,nodes[["p_lt_wt"]])
  nodes[["rate_diff"]]=sapply(1:dim(nodes)[1],function(i) mean(post.dist$lambda[,i]-post.dist$lambda[,1]))
  if(!is.null(fitres$k) && length(fitres$k$mean)>0){
    for(x in c("mean","sd","lb","ub","median")){
      nodes[[sprintf("k_%s",x)]]=fitres$k[[x]]
    }
  }
  if(!is.null(fitres$p) && length(fitres$p$mean)>0){
    for(x in c("mean","sd","lb","ub","median")){
      nodes[[sprintf("p_%s",x)]]=fitres$p[[x]]
    }
  }
  nodes
}




do_paper_plot_fig2_3=function(PDD,group){
  if(group=="a"){
    patients=c("PD7271","PD5163","PD5117","PD5182","PD5179")
  }else if(group=="b"){
    patients=c("PD5847","PD9478","PD4781","PD6629","PD6646")
  }else{
    stop("group should be \"a\" or \"b\"")
  }
  nf=layout(matrix(1:6,ncol=2,byrow = FALSE),width=c(0.25,0.25))
  PDD=PDD[patients]
  colonyinfo=get_patient_timepoints()
  for(i in 1:length(PDD)){
    PD=PDD[[i]]
    tmp=colonyinfo %>% filter(patient==PD$patient)
    age1=tmp$AGE_AT_DIAGNOSIS[1]
    subtype=tmp$SUBTYPE[1]
    do_patient_paper_plot(PD,left.margin.prop = 0.2);title(sprintf("%s (%s age %d)",PD$patient,subtype,floor(age1)))
  };#plot.new();plot.new()
  plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
  add_legend2("center")
}

do_lambda_plot=function(fl,PDD,ftype="local",p.low.or.high="high"){
  res=get_lambda_plot(fl,PDD,ftype,p.low.or.high)
  res$fl
}
get_lambda_plot=function(fl,PDD,ftype="local",p.low.or.high="high",b.do.plot=TRUE,b.add.telomeres=TRUE){
  N=dim(fl)[1]
  ##Override driver column
  fl$driver=fl$driver3
  coldf=getcolor_df(data.frame(patient=PATIENTS ))
  
  ##Meta-Analysis
  tmp=fl %>% filter(driver=="WT")
  wtrma=summary(rma(tmp$mean_rescaled,tmp$sd_rescaled))
  print(wtrma)
  rate_meta_summary=sprintf("%3.2f (%3.2f-%3.2f)",wtrma$beta,wtrma$ci.lb,wtrma$ci.ub)
  
  
  #fl2=fl %>% inner_join(df2,by=c("patient","driver"))
  fl=fl %>% left_join(coldf,by="patient")
  fl$patient2=fl$patient
  fl=fl %>% mutate(patient=ifelse(driver=="WT",patient,""))
  MUSTARD="#ffdb58"
  
  fl=fl %>% mutate(col=ifelse(grepl("JAK2",driver),"red",ifelse(grepl("WT",driver),"grey",MUSTARD))) %>%
    mutate(col2=ifelse(grepl("JAK2",driver),"red",ifelse(grepl("WT",driver),"grey","black")))
  if(p.low.or.high=="high"){
    pp=ifelse(fl$patient!="","",ifelse(fl$p_lt_wt<0.025,ifelse(fl$p_lt_wt<0.025/N,"**","*"),""))
    extra="Posterior Prob. > Wild type rate"
    label2="Prob. > WT"
  }else{
    pp=ifelse(fl$patient!="","",ifelse(fl$p_lt_wt<0.001,"<0.001",sprintf("%4.3f",fl$p_lt_wt)))
    extra="Posterior Prob. <= Wild type rate"
    label2="Prob. < WT"
  }
  fl$p_lt_wt[grepl("WT",fl$driver)]=NA
  res=list(fl=fl,rmeta=wtrma,rate_meta_summary=rate_meta_summary)
  
  
  if(!b.do.plot){
    return(res)
  }
  
  bwidth=0.12
  telomere_panal_start=30
  tbs=telomere_panal_start
  if(!b.add.telomeres){
    tpl=0 ##telomere width
    tbs=35
  }else{
    tpl=12
    ##Pull in all the telomere data..
    ##Get sample per timepoint info.
    df=do.call("rbind",lapply(PDD,function(PD) PD$pdx$agedf %>% filter(tip.label!="zeros")))
    df$driver=df$driver3
    dfa=df %>% group_by(patient,age_at_sample_pcy) %>% summarise(N=n())
    dfa=dfa[order(dfa$patient,dfa$age_at_sample_pcy),]
    dfa$sample_idx=do.call("c",lapply(unique(dfa$patient),function(x) 1:length(which(dfa$patient==x))))
    dfx=df %>% inner_join(dfa[,-3],by=c("patient","age_at_sample_pcy"))
    dfx=dfx %>% left_join(dfx %>% group_by(patient) %>% summarise(offset=median(unique(sample_idx-1))),by=c( "patient"))
    df2=df %>% group_by(patient,driver) %>% summarise(mtelo=mean(telo_mean_length))
  }
  
  plot(NULL,xlim=c(-10,tbs+3*tpl),ylim=c(0,N+2),xlab="",yaxt="n",ylab="",main="",xaxt="n",bty="n")
  
  with(fl,rect(xleft = lb_rescaled,xright=ifelse(ub_rescaled>tbs,tbs,ub_rescaled),ybottom = (N:1)-bwidth,ytop =(N:1)+bwidth, col=col,border=NA))
  with(fl,points(y=N:1,x=median_rescaled,cex=1,pch=15))
  with(fl,text(y=N:1,x=rep(7,N),labels = sprintf("%s",driver),col="black",pos = 2))
  with(fl,text(y=N:1,x=rep(-10,N),labels = patient,col="black",pos = 4,offset=0.5))
  with(fl,text(y=N:1,x=rep(tbs,N),labels = pp,col="black",pos = 2))
  with(fl,text(y=N:1,x=rep(10,N),labels = colony_count,col="black",pos = 2,cex=1))
  abline(v=7)
  axis(side=1,at=c(-10,tbs+3*tpl),labels=rep("",2),lwd.ticks = 0)
  axis(side = 3,at=c(-10,tbs+3*tpl),labels=rep("",2),lwd.ticks = 0)
  abline(v=-10)
  
  axis(side = 1,at=seq(10,tbs-5,5))
  segments(x0=wtrma$b,y0=-10,y1=N+0.5,lty="dotted")
  text(x=wtrma$b,y = N+0.5,labels =rate_meta_summary ,pos = 3,cex=1,offset=0.2)
  segments(x0=tbs+tpl*0,y0=-10,y1=N+0.5)
  pchh=c(5,2,6)
  pchm=c(23:25)
  abline(v=10)
  abline(v=tbs)
  text(y=N+2,x=20,labels="Mutation Rate (SNV/Year)",pos=NULL,cex=1.2)
  text(y=N+2,x=8.5,labels="N",pos=NULL,cex=1.2)
  devnull=sapply(which(nchar(fl$patient)>1),function(i) segments(x0=-10,x1=tbs+3*tpl,y0=N-i+1+0.5,lwd=0.5))
  if(b.add.telomeres){
    abline(v=tbs+3*tpl)
    segments(x0=tbs+tpl*1,y0=-10,y1=N+0.5)
    for(k in 0:2){
      segments(x0=tbs+tpl*k,y0=-10,y1=N+0.5)
      text(y=N+1,x=tbs+tpl*k+1,labels=sprintf("Timepoint %s",k+1),col="black",pos=4)
      for(i in 1:N){
        tmp=dfx %>% filter(patient==fl$patient2[i],driver==fl$driver[i],sample_idx==k+1)
        points(x=tbs+tpl*k+(tpl/8)*tmp$telo_mean_length/1000,y=N-i+1+0.25*(runif(dim(tmp)[1])-0.5),pch=19,cex=0.5,col=fl$col[i])
        points(x=tbs+tpl*k+(tpl/8)*mean(tmp$telo_mean_length/1000),y=N-i+1,cex=1.2,pch=23,col="black",bg=fl$col[i])
      }
      tscale=seq(0,6,2)
      axis(side = 1,at=tbs+tpl*k+tscale*(tpl/8),labels=tscale)
    }
    ps=tbs+3*tpl+1
    ns=ps+tpl
    text(y=N+2,x=tbs+(tpl/2)*3,labels="Telomere Length(kb)",pos=NULL,cex=1.2)
  }
  return(res)
}

get_patient_timepoints=function(){
  res=get_sample_age_table() %>%
    group_by(patient,age_at_sample_exact,age_at_sample,DOB,DATE_OF_DIAGNOSIS,SAMPLE_DATE,SUBTYPE) %>%
    summarise(N=n()) %>%
    mutate(AGE_AT_DIAGNOSIS=time_length(as.period(interval(DOB, DATE_OF_DIAGNOSIS)), "years"),
           age_at_diagnosis_pcy=(as.integer(DATE_OF_DIAGNOSIS)-as.integer(DOB)+MEAN_AGE_AT_DELIVERY_DAYS)/365.25,
           age_at_sample_pcy=(as.integer(SAMPLE_DATE)-as.integer(DOB)+MEAN_AGE_AT_DELIVERY_DAYS)/365.25
    )
  res
}

collate_driver_info=function(PDD,treemodel="poisson_tree",b.is.null=FALSE){
  ##Does the messy collation of node timings
  type=ifelse(b.is.null,"nullmodel","altmodel")
  #cat("WARNING: Review the labelling of the output\n")
  inf1=do.call("rbind",lapply(PDD,function(x){
    out=x$fit[[treemodel]][[type]]$summary %>% mutate(patient=x$patient)
    out=out %>% left_join(add_cumulative_mutcount(x,x$fit[[treemodel]][[type]]$summary$node))
    idx=match("patient",colnames(out))
    cbind(data.frame(patient=out[,idx],stringsAsFactors = FALSE),out[,-idx])
    ##Add
  }
    ))

  ##Make driver correspond to events on the current node
  inf1$driver=gsub(":.*","",inf1$driver3)
  ##Rename 9pUPD_A etc
  inf1$driver=gsub("_[A-Z]$","",inf1$driver)
  inf1$driver=gsub("_[a-z]$","",inf1$driver)
  #Filter to non-private branches that are of interest (non-nuisance)
  inf2=inf1 %>% filter(status==1 | child_count>1)

  #The driver scheme maps specific mutations to drivers..
  ds=get_driver_scheme()
  ##The following is non-missing for copy number events
  inf2$group=ds$group[match(inf2$driver,ds$group)]
  ##Specific genes are filled in order or priority
  ## The driver column can have multiple drivers - so we prioritise
  inf2$group[grep("CBL",inf2$driver)]="CBL"
  inf2$group[grep("TET2",inf2$driver)]="TET2"
  inf2$group[grep("PPM1D",inf2$driver)]="PPM1D"
  inf2$group[grep("DNMT3A",inf2$driver)]="DNMT3A"
  inf2$group[grep("JAK2",inf2$driver)]="JAK2"
  pti2=get_patient_timepoints()
  rti=get_recap_timepoints()
  #browser()
  pti=rbind(as.data.frame(pti2[,c("patient","age_at_sample_pcy")]),as.data.frame(rti[,c("patient","age_at_sample_pcy")]))
  df=pti %>% group_by(patient) %>% summarise(N=n(),max_age_at_sample=max(age_at_sample_pcy),min_age_at_sample=min(age_at_sample_pcy))
  #browser()
  inf= inf2 %>% left_join(df,by="patient") %>% left_join(pti2 %>% group_by(patient,age_at_diagnosis_pcy) %>% summarise(totcolonies=sum(N)),by="patient")
  ##print(inf)
  inf=inf[order(inf$max_age_at_sample,inf$patient,inf$upper_median),]

  ##The following reinstate 1q+ and CBL as more specific classes (will potentially break with updated driver_scheme!)
  ##Should update driver_scheme instead..
  ds$group[grepl("1q+",ds$driver)]="1q+"
  ds$colour[grepl("1q+",ds$driver)]="brown"
  ds$group[grepl("CBL",ds$driver)]="CBL"
  inf$group[grepl("1q+",inf$driver)]="1q+"
  inf$col=ds$colour[match(inf$group,ds$group)]

  inf=inf %>% filter(!is.na(col))
  inf
}
#inf=collate_driver_info(PDD)
get_recap_timepoints=function(){
  out=Reduce(rbind,lapply(PATIENTS,get_patient_summary)) %>% filter(phase=="Recapture" & !is.na(age_at_sample_exact))
  out$age_at_sample_pcy=out$age_at_sample_exact+MEAN_AGE_AT_DELIVERY_DAYS/365.25
  out
}
do_acquisition_plot2=function(inf,##=collate_driver_info(PDD)
                              pti=get_patient_timepoints(),
                              rti=get_recap_timepoints(),
                              extra="",xmax=100,b.add.text=FALSE){
  patients=unique(inf$patient)
  N=length(patients)
  par(mar=c(5,6,2,2)+0.1)
  plot(NULL,xlim=c(0,xmax*1.15),ylim=c(0,3*length(patients)),yaxt="n",xlab="Post Conception Age (Years)",ylab="",bty="n",xaxt="n")
  title(sprintf("Age of Acquisition of Drivers%s",extra),cex=2)
  axis(side = 1,at = seq(0,90,10),labels = seq(0,90,10))
  k=1
  width=0.2
  xgap=0.3##0.5*(0.5-width)
  kk=0
  pch.sample=2
  pch.diagnosis=18
  yb=0

  for(patient in patients){
    idx=which(inf$patient==patient)
    if(length(idx)>1){
      gap=3/(length(idx)+1)
    }else{
      gap=1.5
    }
    yb=3*kk+gap
    rect(xleft=0,xright=inf$max_age_at_sample[idx[1]],ybottom=3*kk+xgap,ytop=3*(kk+1)-xgap,col="lightgray",border="black",lwd=1)
    for(j in idx){
      rect(xleft = inf$lower_median[j],xright=inf$upper_median[j],ybottom = yb-width,ytop=yb+width,border=NA,col=inf$col[j])
      segments(x0=inf$lower_lb95[j],x1=inf$lower_ub95[j],y0=yb,col="black",lwd=1)
      segments(x0=inf$upper_lb95[j],x1=inf$upper_ub95[j],y0=yb,col="black",lwd=1)
      if(b.add.text){
        text(x=inf$max_age_at_sample[idx[1]]+5,y=yb,labels = inf$driver3[j],pos = 4)
      }
      #points(x=inf$upper_ub95[j],y=yb,col="black",cex=1,pch=19)
      yb=yb+gap
      ##k=k+1
    }
    ymid=3*kk+1.5
    mtext(text = patient,side = 2,line=1,at=ymid,las=2)
    points(y=ymid,x=inf$age_at_diagnosis[j],pch=pch.diagnosis,col="black",cex=1.5)
    tmp=sort(unique(rti$age_at_sample_pcy[which(rti$patient==patient)]))
    if(!is.null(tmp)){
      segments(x0=tmp,y0=rep(ymid,length(tmp))-0.2,y1=rep(ymid,length(tmp))+0.2,col="blue",lwd=2,lend=2)
    }
    tmp=sort(unique(pti$age_at_sample_pcy[which(pti$patient==patient)]))
    segments(x0=tmp,y0=rep(ymid,length(tmp))-0.2,y1=rep(ymid,length(tmp))+0.2,col="red",lwd=2,lend=2)
    kk=kk+1
  }
  ##Add a manual legend
  if(!b.add.text){
    tmp=inf %>% group_by(group,col) %>% summarise(N=n())
    yloc=3*kk/4 +3
    #xxl=90
    xxl=0.9*xmax
    offset=5
    for(i in 1:dim(tmp)[1]){
      rect(xleft=xxl,xright=xxl+offset,ybottom=yloc+i,ytop=yloc+i+2*width,col=tmp$col[i],border=NA)
      text(x=xxl+offset,y=yloc+i+width,labels = tmp$group[i],pos=4,offset=1)
    }
    yloc=3*kk/4
    points(x=xxl+0.5*offset,y=yloc,pch=pch.diagnosis,cex=1.5)
    text(x=xxl+offset,y=yloc,pos=4,offset=1,labels="Diagnosis")
    yloc=yloc-2
    segments(x0=xxl+0.5*offset,y0=yloc-0.2,y1=yloc+0.2,col="red",lwd=2,lend=2)
    text(x=xxl+offset,y=yloc,pos=4,offset=1,labels="Colony Sampling")
    yloc=yloc-2
    segments(x0=xxl+0.5*offset,y0=yloc-0.2,y1=yloc+0.2,col="blue",lwd=2,lend=2)
    text(x=xxl+offset,y=yloc,pos=4,offset=1,labels="Recapture Sampling")
  }
  abline(v=MEAN_AGE_AT_DELIVERY_DAYS/365.25,lty="dotted")
  inf
}

do_acquisition_plot2b=function(inf,##=collate_driver_info(PDD)
                              pti=get_patient_timepoints(),
                              rti=get_recap_timepoints(),
                              extra="",xmax=100,b.add.text=FALSE){
  inf2=inf
  patients=unique(inf$patient)
  N=length(patients)
  ## Gestational Age at Birth
  ab=MEAN_AGE_AT_DELIVERY_DAYS/365.25
  ab2=(MEAN_AGE_AT_DELIVERY_DAYS+14)/365.25
  #manipulate ranges so that we can extend pre-birth period...
  browser()
  for(field in c("lower_median","lower_lb95","lower_ub95","upper_median","upper_lb95","upper_ub95","max_age_at_sample","age_at_diagnosis_pcy")){
    inf[[field]]=inf[[field]]-ab
    inf[[field]]=ifelse(inf[[field]]<0,xmax*0.1*inf[[field]]/ab2,inf[[field]])
  }
  pti$age_at_sample_pcy=pti$age_at_sample_pcy-ab
  rti$age_at_sample_pcy=rti$age_at_sample_pcy-ab
  
  par(mar=c(6,6,2,2)+0.1)
  plot(NULL,xlim=c(-xmax*0.1,xmax*1.15),ylim=c(0,3*length(patients)),yaxt="n",xlab="",ylab="",bty="n",xaxt="n")
  if(!is.null(extra)){
  title(sprintf("Age of Acquisition of Drivers%s",extra),cex=2)
  }
  axis(side = 1,at = seq(0,90,10),labels = seq(0,90,10))
  axis(side =1, at = -xmax*0.1+xmax*0.1*c(0,26/40,1),labels=c("0","26",""),cex.axis=0.6,las=1)
  axis(side =1, at = -xmax*0.1+xmax*0.1*c(13/40),labels=c("13"),cex.axis=0.6,las=1)
  mtext(side = 1,line = 4,text = "Age\n(Years)",at = xmax/2)
  mtext(side = 1,line = 4,text = "Gestational Age\n(Weeks)",at = -0.1*xmax)
  k=1
  width=0.2
  xgap=0.3##0.5*(0.5-width)
  kk=0
  pch.sample=2
  pch.diagnosis=18
  yb=0
  
  
  
  
  for(patient in patients){
    idx=which(inf$patient==patient)
    if(length(idx)>1){
      gap=3/(length(idx)+1)
    }else{
      gap=1.5
    }
    yb=3*kk+gap
    rect(xleft=0,xright=inf$max_age_at_sample[idx[1]],ybottom=3*kk+xgap,ytop=3*(kk+1)-xgap,col="lightgray",border=NA,lwd=1)
    rect(xleft=-0.1*xmax,xright=0,ybottom=3*kk+xgap,ytop=3*(kk+1)-xgap,col="lightpink",border=NA,lwd=1)
    
    for(j in idx){
      rect(xleft = inf$lower_median[j],xright=inf$upper_median[j],ybottom = yb-width,ytop=yb+width,border=NA,col=inf$col[j])
      segments(x0=inf$lower_lb95[j],x1=inf$lower_ub95[j],y0=yb,col="black",lwd=1)
      segments(x0=inf$upper_lb95[j],x1=inf$upper_ub95[j],y0=yb,col="black",lwd=1)
      if(b.add.text){
        text(x=inf$max_age_at_sample[idx[1]]+5,y=yb,labels = inf$driver3[j],pos = 4)
      }
      #points(x=inf$upper_ub95[j],y=yb,col="black",cex=1,pch=19)
      yb=yb+gap
      ##k=k+1
    }
    ymid=3*kk+1.5
    mtext(text = patient,side = 2,line=1,at=ymid,las=2)
    points(y=ymid,x=inf$age_at_diagnosis_pcy[j],pch=pch.diagnosis,col="black",cex=1.5)
    tmp=sort(unique(rti$age_at_sample_pcy[which(rti$patient==patient)]))
    if(!is.null(tmp)){
      segments(x0=tmp,y0=rep(ymid,length(tmp))-0.2,y1=rep(ymid,length(tmp))+0.2,col="blue",lwd=2,lend=2)
    }
    tmp=sort(unique(pti$age_at_sample_pcy[which(pti$patient==patient)]))
    segments(x0=tmp,y0=rep(ymid,length(tmp))-0.2,y1=rep(ymid,length(tmp))+0.2,col="red",lwd=2,lend=2)
    kk=kk+1
  }
  ##Add a manual legend
  if(!b.add.text){
    tmp=inf %>% group_by(group,col) %>% summarise(N=n())
    yloc=3*kk/4 +3
    #xxl=90
    xxl=0.85*xmax
    offset=5
    for(i in 1:dim(tmp)[1]){
      rect(xleft=xxl,xright=xxl+offset,ybottom=yloc+i,ytop=yloc+i+2*width,col=tmp$col[i],border=NA)
      text(x=xxl+offset,y=yloc+i+width,labels = tmp$group[i],pos=4,offset=1)
    }
    yloc=3*kk/4
    points(x=xxl+0.5*offset,y=yloc,pch=pch.diagnosis,cex=1.5)
    text(x=xxl+offset,y=yloc,pos=4,offset=1,labels="Diagnosis")
    yloc=yloc-2
    segments(x0=xxl+0.5*offset,y0=yloc-0.2,y1=yloc+0.2,col="red",lwd=2,lend=2)
    text(x=xxl+offset,y=yloc,pos=4,offset=1,labels="Colony Sampling")
    yloc=yloc-2
    segments(x0=xxl+0.5*offset,y0=yloc-0.2,y1=yloc+0.2,col="blue",lwd=2,lend=2)
    text(x=xxl+offset,y=yloc,pos=4,offset=1,labels="Recapture Sampling")
  }
  #abline(v=MEAN_AGE_AT_DELIVERY_DAYS/365.25,lty="dotted")
  abline(v=0,lty="dotted")
  inf2
}







get_drivergroup_colors=function(alpha=0.8){
  MUSTARD="#ffdb58"
  #data.frame(colour=c("#E41A1C", "#999999","#4DAF4A"),#   brewer.pal(3,"Set1"),
  data.frame(colour=sapply(c("red3",MUSTARD,"grey"),function(x) adjustcolor(x,alpha=alpha)),
             col2=c("red3",MUSTARD,"grey"),
             driver2=c("JAK2","Other","WT"),
             desc=c("JAK2","Other Drivers","Wild Type"),
             stringsAsFactors = FALSE)
}

##Pass in data.frame and plot the desired field by copying into feature
# e.g. via df %>% mutate(feature=mean_length)
plot_property_avg=function(df,label,jscale=0.5,yscalediv=100,b.just.show.informative=TRUE,xscalediv=5,yylim=NULL,xxlim=NULL,ccex=0.5,b.legend.inside=FALSE,rscale=0.5,asr=0.2,b.add.braces=TRUE){
  coldf=get_drivergroup_colors(alpha=0.8)
  df$colour=coldf$colour[match(df$driver2,coldf$driver2)]
  ##The following requires at WT and JAK2 to be present (see final filter)
  tmp=df %>% group_by(patient,age_at_sample,driver2) %>% summarise(N=n()) %>% group_by(patient,age_at_sample) %>% summarise(rg=length(which(.data$driver2!="Other"))) %>% filter(rg>1)
  df2=tmp %>% left_join(df,by=c("patient","age_at_sample"))
  dat2=df2 %>%
    group_by(patient,age_at_sample,driver2) %>%
    summarise(feature_mean=mean(feature)) %>%
    spread(key="driver2",value="feature_mean")
  dat2$age=dat2$age_at_sample
  if(b.just.show.informative){
    tmp=df2
  }else{
    tmp=df
  }
  tmp$age=tmp$age_at_sample+jscale*(runif(length(tmp$age_at_sample))-0.5)
  if(is.null(yylim)){
    ymin=0
    ymax=yscalediv*ceiling(max(tmp$feature)/yscalediv)
  }else{
    ymin=yylim[1]
    ymax=yylim[2]
  }
  if(is.null(xxlim)){
    xmin=xscalediv*floor(min(tmp$age)/xscalediv)
    xmax=xscalediv*ceiling(max(tmp$age)/xscalediv)
  }else{
    xmin=xxlim[1]
    xmax=xxlim[2]
  }
  par(mar=c(5,4,4,7)+0.1,xpd=2)
  with(tmp,plot(age,
                feature,
                col=colour,pch=19,cex=ccex,xlab="Age (Years)",ylab=label,main=sprintf("%s vs Driver Status",label),
                xlim=c(xmin,xmax),ylim=c(ymin,ymax),las=1
  ))
  if(TRUE){
    par(xpd=FALSE)
    grid()
  }
  JAK2col=coldf$col2[match("JAK2",coldf$driver2)]
  WTcol=coldf$col2[match("WT",coldf$driver2)]
  yscale=(asr*rscale/diff(par("usr"))[1])*par("usr")[4]

  with(dat2,points(x=age,y=WT,bg=WTcol,cex=2,pch=21,col="black"))
  with(dat2,points(x=age,y=JAK2,bg=JAK2col,cex=2,pch=21,col="black"))

  tmp$col="black"
  if(b.add.braces) add_braces(tmp,"feature")

  coldf=coldf %>% filter(driver2 %in% df$driver2)
  if(b.legend.inside){
    legend("topright",coldf$desc,col=coldf$colour,pch=19,bty="n")
  }else{
    legend(x=par("usr")[2]+1,y=0.75*ymax,coldf$desc,col=coldf$colour,pch=19,ncol = 1,cex=1,inset=c(-0.1,0.5),bty="n")
  }
  zz=dat2 %>% group_by(patient) %>% summarise(JAK2=mean(JAK2),WT=mean(WT),age_mean=mean(age))
  print(wilcox.test(zz$JAK2-zz$WT))
  zz
}
add_braces=function(dat2,yfield){
  y=par("usr")
  ygmin=y[3]
  ygmax=y[4]
  ymax=max(dat2[[yfield]])
  ymin=min(dat2[[yfield]])

  d3=dat2 %>%
    group_by(patient,age_at_sample,col) %>%
    summarise(NN=n(),mina=min(age_at_sample))
  d3=d3[order(d3$mina),]
  dat2=d3 %>% left_join(d3 %>% group_by(patient) %>% summarise(N=n()))

  dat2=dat2 %>% filter(N>1)

  voff=ygmax/20
  patients=unique(dat2$patient)
  np=length(patients)
  if(np==0){
    return(NULL)
  }
  if(ygmax-ymax>ymin-ygmin){
    cb=seq(ymax,ygmax,length.out = np+2)[2:(np+1)]
    orientation=-1
  }else{
    cb=seq(ygmin,ymin,length.out = np+2)[2:(np+1)]
    orientation=1
  }


  for(i in 1:np){
    d4=dat2 %>% filter(patient==patients[[i]])
    ages=d4$age_at_sample
    NN=length(ages)
    arrows(x0=min(ages),x1=max(ages),y0=cb[i],length=0,lwd=2,lend=1,col=d4$col)
    arrows(y0=rep(cb[i],NN),y1=rep(cb[i]+orientation*voff,NN),x0=ages,length=0,lend=2,lwd=2,col=d4$col)
  }
}

assign_genotypes_mtr=function(tree,dat,maxits=5){
  #df=reconstruct_genotype_summary(best_tree)
  mtr=dat$mtr
  dep=dat$dep
  mtr=cbind(mtr,zeros=0)
  dep=cbind(dep,zeros=10)
  #if(is.null(pdx$meta$error_rate)){
  #  cat("Assuming error rate of 0.01 - can configure this on per colony basis in inf$meta$error_rate parallel to inf$meta$clones_short\n")
  #  pdx$meta$error_rate=rep(0.01,length(pdx$meta$clones_short))}
  #browser()
  p.err=ifelse(tree$tip.label=="zeros",1e06,0.01)
  #p.err=c(pdx$meta$error_rate,1e-6)
  #p.err=p.err[match(tree$tip.label,c(pdx$meta$clones_short,"zeros"))]
  info=assign_to_tree(tree,mtr,dep,error_rate=p.err,maxits=maxits)
  dat$tree_ml=info$tree
  dat$summary=info$summary
  dat$df=info$df
  dat
}

remapvariants=function(PD,type="snv"){
  nn=get_normal(PD$inf)
  idx=with(PD$inf[[type]]$details,which(filter_bglod==1 & !FILTER_EXCEPT_GLOD & is_localx_excluded==0))
  tree=PD$pdx$tree_ml
  rr=assign_genotypes_mtr(tree,list(mtr=PD$inf[[type]]$mtr[idx,],dep=PD$inf[[type]]$dep[idx,]))
  ntip=length(tree$tip.label)
  germline.idx=which(tree$edge[,1]==ntip+1 & tree$edge[,2]>ntip)
  idx2=which(rr$summary$edge_ml !=germline.idx)
  cbind(rr$summary[idx2,],idx=idx[idx2],mtr_n=PD$inf[[type]]$mtr[idx[idx2],nn],
        dep_n=PD$inf[[type]]$dep[idx[idx2],nn],
        vaf_n=PD$inf[[type]]$mtr[idx[idx2],nn]/PD$inf[[type]]$dep[idx[idx2],nn])

}

add_back_in_variants=function(PD,rescuedvariants){
  ##Assume these have idx in PD$inf$snv
  rv=rescuedvariants
  rv=rv %>% filter(pval>1e-10)
  if(length(rv$idx)==0){
    return(PD)
  }else{
    cat("Restoring",length(rv$idx)," variants previously thought germline:")
    print(PD$inf$snv$details[rv$idx,1:10])
  }
  columns=colnames(PD$pdx$dat$details)
  PD$inf$snv$details$keep_embryonic[rv$idx]=1
  PD$inf$snv$details$FILTER[rv$idx]=0
  det=PD$inf$snv$details[rv$idx,]
  det$TYPE="SNV"
  det$node=PD$pdx$tree_ml$edge[rv$edge_ml,2]
  det=det %>% mutate(key=sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
  det=det[,colnames(PD$pdx$dat$details)]
  pdxnames=colnames(PD$pdx$dat$mtr)
  inf=PD$inf$snv
  inf$geno_colony=inf$geno
  for(field in c( "mtr","wtr","dep","ofs","minmutvaf", "geno_colony" )){
    cat(field,"\n")
    extrafields=setdiff(colnames(PD$pdx$dat[[field]]),colnames(inf[[field]]))
    if(length(extrafields)>0){
    empty=matrix(NA,nrow=length(rv$idx),ncol=length(extrafields))
    colnames(empty)=extrafields
    inf[[field]]=cbind(inf[[field]][rv$idx,],empty)
    }else{
      inf[[field]]=inf[[field]][rv$idx,]
    }
    PD$pdx$dat[[field]]=rbind(PD$pdx$dat[[field]],inf[[field]][,colnames(PD$pdx$dat[[field]])])
  }
  PD$pdx$dat$details=rbind(PD$pdx$dat$details,det)
  rv$ml=rv$edge_ml
  PD$pdx$summary=rbind(PD$pdx$summary,rv[,colnames(PD$pdx$summary)])
  tmp=PD$pdx$dat$details
  tmp$chridx=match(PD$pdx$dat$details$Chrom,c(sprintf("%s",1:22),"X","Y"))
  idx=order(tmp$chridx,tmp$Pos,tmp$Ref,tmp$Alt)
  for(field in c( "mtr","wtr","dep","ofs","minmutvaf", "geno_colony" )){
    PD$pdx$dat[[field]]=PD$pdx$dat[[field]][idx,]
  }
  PD$pdx$dat$details=PD$pdx$dat$details[idx,]
  PD$pdx$summary=PD$pdx$summary[idx,]
  ##Redo edge length

  PD$pdx$tree_ml$el.snv=sapply(PD$pdx$tree_ml$edge[,2],function(node) length(which(PD$pdx$dat$details$node==node & PD$pdx$dat$details$TYPE=="SNV")))
  PD$pdx$tree_ml$el.snv.global.filtered=sapply(PD$pdx$tree_ml$edge[,2],function(node) length(which(PD$pdx$dat$details$node==node & PD$pdx$dat$details$is_globalx_excluded==0 & PD$pdx$dat$details$TYPE=="SNV")))
  PD$pdx$tree_ml$el.snv.local.filtered=sapply(PD$pdx$tree_ml$edge[,2],function(node) length(which(PD$pdx$dat$details$node==node & PD$pdx$dat$details$is_localx_excluded==0 & PD$pdx$dat$details$TYPE=="SNV")))
  if(!is.null(PD$pdx$tree_ml$per.branch.sensitivity.reg)){
    PD=add_adjustment_models(PD)
  }
  PD
}




add_cumulative_mutcount=function(PD,nodes){
  tree=PD$pdx$tree_ml
  tree$edge.length=PD$localx.correction2*tree$el.snv.local.filtered/tree$per.branch.sensitivity.hybrid.multi
  nh=nodeHeights(tree)
  mutcount=matrix(nh[match(nodes,tree$edge[,2]),],ncol=2)
  colnames(mutcount)=c("lower_mutcount_adj","upper_mutcount_adj")
  tree$edge.length=tree$el.snv.local.filtered
  nh=nodeHeights(tree)
  mutcount2=matrix(nh[match(nodes,tree$edge[,2]),],ncol=2)
  colnames(mutcount2)=c("lower_mutcount","upper_mutcount")
  cbind(data.frame(node=nodes),as.data.frame(mutcount),as.data.frame(mutcount2))
}

get_simple_treemodel_timings=function(PDD,b.is.null){

  global_rate=do.call("rbind",lapply(PDD,function(PD){
    tmp=get_stan_result_df(PD,model = treemodel,bIsNull = TRUE)
    tmp$patient=PD$patient
    tmp[,c("patient","mean_rescaled","sd_rescaled")]
  }))

  dat=collate_driver_info(PDD,treemodel = "poisson_tree",b.is.null = b.is.null)
  ##Sensitivity adjustment...
  dat=dat %>% mutate(sensitivity_u=upper_mutcount/upper_mutcount_adj,sensitivity_l=lower_mutcount/lower_mutcount_adj)
  dat$sensitivity_l=ifelse(is.na(dat$sensitivity_l),1,dat$sensitivity_l)
  dat=dat %>% left_join(global_rate)
  #Assuming a Poisson distribution then maximum likelihood estimator for the variance is just the observed count.
  #Under the  conservative assumption that the distribution for lambda (MCMC estimate) and the estimation of lambda*s*t are independent (neg correlated?)
  #and that s is known precisely then the distribution of t is the ratio of approximately normal distributions
  #   t=lambda*/(s*lambda)  with variance lambda*
  #  if z=x/y then z_var=(x_m/y_m)**2(x_var/x_m**2)+(y_var/y_m**2)) (see https://en.wikipedia.org/wiki/Ratio_distribution)
  varfn=function(k,s,l,l_sd){((k**2)/((s*l)**2))*((1/k)+((l_sd**2)/l**2))}
  # The following simple example indicates that this works as expected..
  #TT=50;rate=rnorm(100,mean = 20,sd = 1);k=rpois(100,0.8*rate*TT)
  #vv=varfn(k,0.8,20,1)
  # Expect the following to be approximately equal (mean estimated variance and actual variance of the estimate timing)
  #mean(vv)
  #var(k/(0.8*20))


  ##Subtract 33 to take into account the developmental mutations
  dat=dat %>% mutate(t_upper_pois=(upper_mutcount-33)/(sensitivity_u*mean_rescaled),
                     t_upper_pois_var=varfn(upper_mutcount,sensitivity_u,mean_rescaled,sd_rescaled),
                     t_lower_pois=(lower_mutcount-33)/(sensitivity_l*mean_rescaled),
                     t_lower_pois_var=varfn(lower_mutcount,sensitivity_l,mean_rescaled,sd_rescaled)
                     ) %>%
      mutate(t_upper_pois_lb95=t_upper_pois-1.96*sqrt(t_upper_pois_var),
             t_upper_pois_ub95=t_upper_pois+1.96*sqrt(t_upper_pois_var),
             t_lower_pois_lb95=t_lower_pois-1.96*sqrt(t_lower_pois_var),
             t_lower_pois_ub95=t_lower_pois+1.96*sqrt(t_lower_pois_var)
             )


  dat

}

plot_timing_comparison=function(dat,xxlim=c(0,80),yylim=c(0,80),...){
  coldf=getcolor_df(dat)
  dat$col=NULL
  dat=dat %>% left_join(coldf)
  plot(NULL,xlim=xxlim,ylim=yylim,
       xlab="MCMC Tree Model",ylab="Poisson Model/MCMC Rate",main="MCMC Tree Model Timings vs Simple Poisson Model(/MCMC Rate)")
  abline(a=0,b=1)
  segments(y0=dat$t_upper_pois,x0=dat$upper_lb95,x1=dat$upper_ub95,col=dat$col,...)
  segments(y0=dat$t_upper_pois_lb95,y1=dat$t_upper_pois_ub95,x0=dat$upper_mean,col=dat$col,...)
  legend("bottomright",coldf$patient,col=coldf$col,lwd=2)

}

###  QC
extra_qc_plot_minimal=function(dat,scaling=1){
  prefix=dat$inf$meta$prefix
  dat=add_excluded_region_flag(dat)
  tmp=apply_adjustment_germline(dat)
  depth=colMeans(dat$pdx$dat$dep[,setdiff(dat$pdx$tree_ml$tip.label,"zeros")],na.rm=TRUE)
  ta1=tmp$pdx_snv_adj
  par(mfcol=c(2,1))
  ta1$tree_ml$edge.length=ta1$tree_ml$el.snv.local.filtered
  tree=plot_basic_tree(ta1,label = sprintf("%s: SNV With No Branch Length Adjustment (Autosomal - Excludes CNA/LOH Regions)",prefix),genes=GENES,cv=CV,bars=depth)
  ta1$tree_ml$edge.length=ta1$tree_ml$el.snv.local.filtered/ta1$tree_ml$sensitivity.snv.local.filtered
  tree=plot_basic_tree(ta1,label = sprintf("%s: SNV (Autosomal - Excludes CNA/LOH Regions): Scale Adjusted Assuming Somatic=Germline Sensitivity",prefix),genes=GENES,cv=CV,bars=depth)
  
  plot_chromosome_annotated_tree_pdx(dat$pdx,sprintf("%s: Chromosome composition check",prefix))
  ##
  par(mfcol=c(2,1))
  par(mar=c(5,4,4,1)+0.1)
  #plot_vaf_by_branch(dat$pdx,prefix,b_pooled=FALSE,cex.label = 0.8)
  plot_vaf_by_branch(dat$pdx,prefix,b_pooled=TRUE,cex.label = 0.8)
  tree=plot_pooled_vaf_tree(dat$pdx,sprintf("%s [Pooled VAF]",prefix),lwd=scaling)
  #if(!is.null(dat$pdx$tree_ml$BS)){
    minmuts=round(0.02*max(nodeHeights(dat$pdx$tree_ml)))
    tmp=dat$pdx
    tmp$tree_ml$color=ifelse(tmp$tree_ml$edge.length<minmuts,"grey","black")
    tmp$tree_ml$edge.length=ifelse(tmp$tree_ml$edge.length<minmuts,minmuts,tmp$tree_ml$edge.length)
    tree=plot_basic_tree(tmp,sprintf("%s: With MPBoot Bootstrap Support (Grey branches extended to %s Mutations)",prefix,minmuts))
    bs_labels(tree)
  #}
}

plot_vaf_by_branch=function(pdx,label,cex.label=1,b_pooled=T,vtype=c("SNV"),b_horizontal=FALSE,b.extra=TRUE){
  
  dd=get_vaf_info_by_branch(pdx,vtype)
  
  if(!b_pooled){
    idx.order=order(sapply(dd,function(x) x$ns))
    dd=dd[idx.order]
    ns=sapply(dd,function(x) x$ns)
    avafs=lapply(dd,function(x) x$vaf)
    if(b_horizontal){
      par(mar=c(5,8,4,2)+0.1)
      bnames=sprintf("%s:%s [%s]",ns,pdx$tree_ml$edge[idx.order,2],sapply(avafs,length))
      las=1
      xlab="VAF"
      ylab=""
    }else{
      bnames=sprintf("%s:%s\n[%s]",ns,pdx$tree_ml$edge[idx.order,2],sapply(avafs,length))
      las=1
      xlab="Edge"
      ylab="VAF"
    }
    
    idx.keep=which(ns>1)
    suppressWarnings(boxplot(avafs[idx.keep],names=bnames[idx.keep],notch=T,
                             main=sprintf("%s:VAF Distribution By Edge [Shared Branches]",label),
                             xlab=xlab,ylab=ylab,horizontal=b_horizontal,pars=list(cex.axis=cex.label),las=las))
    if(b_horizontal){
      abline(v=0.5)
    }else{
      abline(h=0.5)
    }
    return(NULL)
  }
  dcnt=unique(sapply(dd,function(x) x$ns))
  dcnt=sort(dcnt[dcnt>0])
  avafs=lapply(dcnt,function(d){idx=which(sapply(dd,function(x) x$ns)==d);do.call("c",lapply(dd[idx],function(x) x$vaf))})
  nv=sapply(avafs,length)
  bnames=sprintf("%s\n[%s]",dcnt,nv)
  suppressWarnings(boxplot(avafs,names=bnames,notch=T,
                           main=sprintf("%s:VAF Distribution By No. Colonies in Clade (Per Sample/Locus)",label),
                           xlab="No of colonies in clade",ylab="VAF",pars=list(cex.axis=cex.label)))
  abline(h=0.5)
  if(b.extra){
  mvafs=lapply(dcnt,function(d){idx=which(sapply(dd,function(x) x$ns)==d);do.call("c",lapply(dd[idx],function(x) x$mvaf))})
  bnames=sprintf("%s\n[%s]",dcnt,nv/dcnt)
  suppressWarnings(boxplot(mvafs,names=bnames,notch=T,
                           main=sprintf("%s:VAF Distribution By No. Colonies in Clade (Per Locus (Pooled over Samples))",label),
                           xlab="No of colonies in clade",ylab="VAF",pars=list(cex.axis=cex.label)))
  abline(h=0.5)
  }
  
}

get_vaf_info_by_branch=function(pdx,vtype=c("SNV","INDEL")){
  germ.idx=which(pdx$tree_ml$edge[,1]==length(pdx$tree_ml$tip.label)+1)
  lapply(pdx$tree_ml$edge[-germ.idx,2],function(i) {
    idx=get_idx_for_node(pdx$dat,node = i)
    idx=intersect(idx,which(pdx$dat$details$TYPE %in% vtype & pdx$dat$details$is_localx_excluded==0) )
    samples=setdiff(get_samples_in_clade(node=i,pdx$tree_ml),"zeros")
    if(length(samples)>1 && length(idx)>1){
      mvaf=rowSums(pdx$dat$mtr[idx,samples])/rowSums(pdx$dat$dep[idx,samples])
    }else if(length(samples)==1){
      mvaf=pdx$dat$mtr[idx,samples]/pdx$dat$dep[idx,samples]
    }else{
      ##length(idx)==1
      mvaf=sum(pdx$dat$mtr[idx,samples])/sum(pdx$dat$dep[idx,samples])
    }
    list(
      vaf=as.numeric(pdx$dat$mtr[idx,samples]/pdx$dat$dep[idx,samples]),
      ns=length(samples),
      samp=samples,
      mvaf=mvaf
    )
  }
  )
  
}
add_bs=function(dat,nexus.file){
  nexus=read_nexus(nexus.file)
  bs=data.frame(BS=sapply(nexus$splits,function(x) x$BS),
                SAMPLE=sapply(nexus$splits,function(x) paste(sort(x$taxa),collapse=":")))
  
  clades=get_clade_keys(dat$pdx)
  idx1=match(clades$c1,bs$SAMPLE)
  idx2=match(clades$c2,bs$SAMPLE)
  idx=ifelse(is.na(idx1),idx2,idx1)
  if(length(which(is.na(idx)))>0){
    cat("can't find all supporting bipartition setting BS to 0")
    # print(clades[idx,])
  }
  dat$pdx$tree_ml$BS=bs$BS[idx]
  dat$pdx$tree_ml$BS=ifelse(is.na(dat$pdx$tree_ml$BS),0,dat$pdx$tree_ml$BS)
  dat
}

bs_labels=function(tree){
  idx=which(tree$edge[,2]>length(tree$tip.label))
  text(tree$BS[idx],x=tree$coords$b1[idx],y=tree$top-tree$coords$a1[idx],col=ifelse(tree$BS[idx]<80,"red","blue"))
}

read_nexus=function(nexus.file){
  lines=readLines(nexus.file)
  idx1=which(lines=="TAXLABELS")+1
  idx2=which(lines=="END; [Taxa]")-2
  rng=idx1[1]:idx2[1]
  samples=gsub("'","",gsub("^\\[[0-9]+\\] ","",lines[rng]))
  idx1=which(lines=="MATRIX")+1
  idx2=which(lines=="END; [Splits]")-2
  rng=idx1[1]:idx2[1]
  splits=strsplit(lines[rng],"\t")
  #browser()
  splits=lapply(splits,function(x) {
    BS=as.numeric(x[2])
    taxa=unlist(strsplit(gsub(",$","",x[3])," "))
    taxa=taxa[which(nchar(taxa)>0)];list(BS=BS,taxa=samples[as.integer(taxa)])})
  list(taxa=samples,splits=splits)
}

get_clade_keys=function(pdx){
  tips=pdx$tree_ml$tip.label
  key=t(sapply(pdx$tree_ml$edge[,2],function(node){
    s=get_samples_in_clade(node = node,tree = pdx$tree_ml)
    c(c1=paste(sort(s),collapse=":"),c2=paste(sort(setdiff(tips,s)),collapse=":"))
  }))
  data.frame(key)
}


get_node_info=function(PD){
  tree=PD$pdx$tree_ml
  tmp=markup_tree(tree,PD$nodes$node[PD$nodes$status>=0])
  idx=which(tree$edge[,2] %in% which(tree$tip.label != "zeros"))
  ncolony=sapply(0:length(PD$nodes$node[PD$nodes$status>=0]),function(x) length(which(tmp$type[idx]==x)))
  descdf=data.frame(patient=rep(PD$patient,length(which(PD$nodes$status>=0))+1),id=sprintf("%s.%s",PD$patient,0:length(PD$nodes$driver3b[PD$nodes$status>=0])),
desc=c("WT",PD$nodes$driver3b[PD$nodes$status>=0]),stringsAsFactors = FALSE)
  descdf$ncolony=ncolony
  descdf
}


do_dnds=function(PDD){
  
  
  ##Add in filtering
  muts=do.call("rbind",lapply(PDD,function(PD) {
    sample=sprintf("%s_%s",PD$patient,PD$pdx$dat$details$node)
    out=PD$pdx$dat$details[,c("Chrom","Pos","Ref","Alt")]
    driver.nodes=unique(do.call("c",lapply(get_tree_drivers(PD$pdx,drivers=get_drivers("drivers_red.bed"))$node,
                                           function(x) c(x,get_all_node_children(x,tree =PD$pdx$tree_ml )))))
    out$is_in_driver_clade=PD$pdx$dat$details$node %in% driver.nodes
    out$sampleID=sample
    fields=c("VC","DBSNP","G1000_AC","EXAC_AC")
    out=cbind(out,PD$pdx$dat$details[,fields])
    
    colnames(out)=c("chr","pos","ref","mut","is_in_driver_clade","sampleID",fields)
    out[,c("sampleID",fields,"VC","chr","pos","ref","mut")]
  }
  )
  )
  ##return(muts %>% mutate(alt=mut) %>% (function(x) x[,c("sampleID","chr","pos","ref","alt")]))
  if(SPECIES_BUILD=="hg38"){
    load("../data/covariates_hg19_hg38_epigenome_pcawg.rda")
    load("../data/RefCDS_human_GRCh38_GencodeV18_recommended.rda")
    gr_genes<<-gr_genes
    dndscv(muts %>% mutate(alt=mut) %>% (function(x) x[,c("sampleID","chr","pos","ref","alt")]),refdb = RefCDS,cv=covs)
  }else{
    dndscv(muts %>% mutate(alt=mut) %>% (function(x) x[,c("sampleID","chr","pos","ref","alt")]))
  }
}


get_glod=function(mut,ref,e,f=0.05){
  lmNonGermline=ref*log10((f*e/3+(1-f)*(1-e)))+mut*log10(f*(1-e)+(1-f)*e/3)
  lm05=(mut+ref)*log10(0.5*(e/3+(1-e)))
  LOD=lmNonGermline-lm05
  LOD
}

add_timing=function(PD){
  dff=get_all_tree_drivers(PD$pdx,genes=GENES,cv = CV)
  if(dim(dff)[1]>0){
    dff$status=-1
    nodes=do.call("rbind",lapply(dff$node,function(node) {
      idx=which(dff$node==node)
      driver=paste(dff$label[idx],collapse=",")
      data.frame(node=node,driver=driver,status=-1)
    }))
    nodes$child_count=sapply(nodes$node,function(node) length(get_samples_in_clade(node=node,tree=PD$pdx$tree_ml)))
  }else{
    nodes=data.frame(node=integer(),driver=character(),status=integer(),child_count=integer())
  }
  nodes=nodes %>% filter(child_count>1)
  
  PD$nodes=nodes
  
  PD=add_driver_node_info(PD,nodes)
  PD
}

## Experimental code.

check_germline_filtering_and_contamination=function(PD,depthmin=20,depthmax=Inf,vafncutoff=0.2,xxlim=c(0,1),yylim=c(0,1)){
  ## Problematic use cases are:
  # Clonal Sweep - depends on depth of matched normal.   
  # Near clonal sweep  - don't want to eliminate variants near top of shared branch that are partly shared with the normal.   
  # Filtering process is currently permissive wrt to germline variants.  We can then remove germline branch or not.  If we see evidence of sharing
  # with normal in the top of the tree then we eliminate the trunk.
  # do visual clustering to assess true level of "contamination" in normal. 
  # Look at top of the tree. In the absence of a clonal trunk we should see some sharing at the top of the tree.
  normal=get_normal(PD$inf)
  message("Normal is ",normal)
  clones=PD$inf$meta$clones_short
  mdepth=mean(PD$inf$snv$dep[,normal])
  depthmin=max(depthmin,quantile(PD$inf$snv$dep[,normal],prob=0.1))
  depthmax=min(depthmax,quantile(PD$inf$snv$dep[,normal],prob=0.9))
  message("Mean normal depth is ",mean(PD$inf$snv$dep[,normal]))
  message("Filtering to normal depth in ",depthmin, " to ",depthmax)
  idx=which(PD$inf$snv$details$is_localx_excluded==0 & PD$inf$snv$dep[,1]>=depthmin & PD$inf$snv$dep[,1]<=depthmax)
  vn=PD$inf$snv$mtr[,normal]/PD$inf$snv$dep[,normal]
  vt=rowSums(PD$inf$snv$mtr[,clones])/rowSums(PD$inf$snv$dep[,clones])
  idx=which(PD$inf$snv$details$is_localx_excluded==0 & PD$inf$snv$dep[,1]>=depthmin & PD$inf$snv$dep[,1]<=depthmax & vn<vafncutoff)
  
  trunk=PD$pdx$tree_ml$edge[get_germline_idx(PD$pdx),2]
  
  ## highlight which ones are on the trunk.
  PD$pdx$dat$details=PD$pdx$dat$details %>% mutate(key=sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
  PD$inf$snv$details=PD$inf$snv$details  %>% mutate(key=sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
  mapped_to_trunk=PD$inf$snv$details$key %in% PD$pdx$dat$details$key[which(PD$pdx$dat$details$node==trunk)]
  ## Add in matched caveman if available.
  matchedcalls=read.table(sprintf("../data/%s_matched_calls.bed",PD$pdx$meta$prefix),head=F,stringsAsFactors = FALSE)
  colnames(matchedcalls)=c("Chrom","Pos","Ref","Alt")
  matchedcall= matchedcalls %>% mutate(key=sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
  matched=PD$inf$snv$details$key %in% matchedcall$key
  res=infer_contamination(vt[idx],vn[idx])
  cn=ifelse(res$cn>1,1,res$cn)
  cn2=res$cn2
  smoothScatter(vn,vt,nrpoints = Inf,xlim=xxlim,ylim=yylim,xlab="VAF (Normal)",ylab="VAF (Tumour [Combined Colonies])")
  browser()
  points(vn[idx],vt[idx],col="brown",cex=0.4)
  idx=which(mapped_to_trunk)
  points(vn[idx],vt[idx],col="red",cex=2,pch=1)
  idx=which(matched)
  points(vn[idx],vt[idx],col="blue",cex=1.5,pch=2)
  #points(x=0.5*cn[2],y=0.5*cn[1],pch=19,cex=2,col="red")
  #points(x=0.5*cn2[2],y=0.5*cn2[1],pch=19,cex=2,col="green")
  
  res
}


add_matched_missing_flag=function(PD){
  matchedcalls=read.table(sprintf("../data/%s_matched_calls.bed",PD$pdx$meta$prefix),head=F,stringsAsFactors = FALSE)
  colnames(matchedcalls)=c("Chrom","Pos","Ref","Alt")
  matchedcall= matchedcalls %>% mutate(key=sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
  PD$pdx$dat$details$is_in_matched_callset=FALSE
  idx=with(PD$pdx$dat$details,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt) %in% matchedcall$key)
  if(length(idx)==0){
    stop("No matching matched calls!")
  }
  PD$pdx$dat$details$is_in_matched_callset[idx]=TRUE
  ##Populate inf as well..
  PD$inf$snv$details$is_in_matched_callset=FALSE
  idx=with(PD$inf$snv$details,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt) %in% matchedcall$key)
  if(length(idx)==0){
    stop("No matching matched calls!")
  }
  PD$inf$snv$details$is_in_matched_callset[idx]=TRUE
  
  
  trunk=PD$pdx$tree_ml$edge[get_germline_idx(PD$pdx),2]
  
  ## Exclude those th
  details=PD$pdx$dat$details %>% filter(is_in_matched_callset)
  PD$pdx$tree_ml$edge.length.matched=sapply(PD$pdx$tree_ml$edge[,2],function(x) length(which(details$node==x & details$TYPE=="SNV")))
  PD$pdx$tree_ml$edge.length=sapply(PD$pdx$tree_ml$edge[,2],function(x) with(PD$pdx$dat$details,length(which(node==x & TYPE=="SNV"))))
  
  # Reglod
  PD$pdx$dat$details$NGLOD=get_glod(PD$pdx$dat$mtr[,1],PD$pdx$dat$dep[,1]-PD$pdx$dat$mtr[,1],0.001,0.005)
  PD$pdx$dat$details= PD$pdx$dat$details %>% mutate(passes_new_glod=ifelse(PD$pdx$dat$details$GSNP==0,NGLOD>1.35,NGLOD>5))
  ##PD$pdx$tree_ml$edge.length.matched=sapply(PD$pdx$tree_ml$edge[,2],function(x) length(which(details$node==x)))
  details=PD$pdx$dat$details %>% filter(passes_new_glod)
  PD$pdx$tree_ml$edge.length.newglod=sapply(PD$pdx$tree_ml$edge[,2],function(x) length(which(details$node==x )))
  PD$pdx$tree_ml$edge.length.orig=PD$pdx$tree_ml$edge.length
  PD$pdx$tree_matched=PD$pdx$tree_ml %>% (function(tree){tree$edge.length=tree$edge.length.matched;tree})
  PD$pdx$tree_glod=PD$pdx$tree_ml %>% (function(tree){tree$edge.length=tree$edge.length.newglod;tree})
  PD
}




infer_contamination=function(vt,#VAF Tumour
                             vn,#VAF Normal
                             vafn_cutoff=0.25){
  dat=cbind(vt,vn)
  dat=dat[which(vn<vafn_cutoff),]
  dat=dat[rowSums(is.na(dat))==0,]
  cls=kmeans(dat,2)
  centres=cls$centers
  
  ##Pull out the cluster with the largest VAF
  idx=which.max(centres[,1])
  message("Total of ",dim(dat),"in clustering analysis\n")
  message("Found ",length(which(cls$cluster==idx))," in largest VAF cluster\n")
  
  ##Do some kind of validation..
  list(cn=c(at=2*centres[idx,1],an=2*centres[idx,2]),cn2=c(at=2*centres[-idx,1],an=2*centres[-idx,2]))
}


plot_all_sv_on_tree=function(PD){
  testing=readBrass(PD)
  ## get rid of germline
  testing=testing %>% left_join(PD$pdx$cfg,by=c("sample"="LABEL")) %>% filter(!is.na(SHORT_LABEL))
  ## discretise to nearest 100
  testing2=testing %>% mutate(s1=ceiling(start1/100)*100,s2=ceiling(start2/100)*100)
  ## Exclude singletons
  tt=testing2 %>% group_by(chr1,s1,chr2,s2) %>% summarise(N=n()) %>% filter(N>1)
  tt=tt[order(tt$N,decreasing = TRUE),]
  cat("found events:\n")
  print(tt)
  for(i in 1:dim(tt)[1]){
    event=testing2 %>% inner_join(tt[i,])
    tree=highlight_tips(PD$pdx$tree_ml,event$SHORT_LABEL,cex.terminal=1,cex.label=0)
    title(sprintf("%s:chr%s:%s-chr%s:%s/%s-%s",PD$patient,event$chr1[1],event$start1[1],event$chr2[1],event$start2[1],event$gene1[1],event$gene2[1]))
  }
}
readBrass=function(PD){
  NST="~/volumes/nw14_network/NFS/cancer_ref01/nst_links/live/"
  cfg=PD$pdx$cfg
  cfg=cfg %>% filter(SHORT_LABEL %in% PD$pdx$tree_ml$tip.label)
  do.call("rbind",lapply(1:dim(cfg)[1],function(i) readBrassBedPE(sprintf("%s/%s/%s/%s.brass.annot.bedpe.gz",NST,cfg$PROJECT[i],cfg$LABEL[i],cfg$LABEL[i]))))
}
readBrassBedPE=function(filespec){
  if(!file.exists(filespec)){
    warning(sprintf("can't find %s",filespec))
    return(NULL)
  }else{
    cat("#")
  }
  lines=readLines(filespec)
  header=gsub("# ","",lines[grep("^# chr",lines)])
  read.table(text=c(header,lines[!grepl("^#",lines)]),sep="\t",head=TRUE,stringsAsFactors = FALSE)
}
highlight_tips=function(tree,tips,...){
  tree$tip.color=rep("black",length(tree$tip.label))
  tree$tip.color[tree$tip.label %in% tips]="red"
  plot_tree(tree,...)
}


plot_phylodyn=function(tree,label=""){
  ut=drop.tip(tree,"zeros")
  ut=multi2di(ut)
  ut=force.ultrametric(ut,method="nnls")
  class(ut)="phylo"
  fullbnpr=BNPR(ut,lengthout=250)
  #par(mfcol=c(2,1))
  #plot(ut,show.tip.label=FALSE,edge.width=0.5)
  plot_BNPR(fullbnpr,col="green",main=label)
  #abline(h=st$cfg$compartment$popsize[2]/2,col="red")
  #abline(h=st$cfg$compartment$popsize[2],col="black")
}
require("data.table")#for frollmean

plot_matrix=function(mat,details,label="",type="l",cnsamples=c()){
  ##browser()
  rng=range(mat,na.rm = TRUE)
  plot(NULL,ylim=rng,xlim=c(0,max(details$Pos)/1e6),ylab=label,xlab="Position(Mb)")
  for(i in 1:dim(mat)[2]){
    if(type=="l"){
      lines(details$Pos/1e6,mat[,i],lwd=0.5)
    }else{
      points(details$Pos/1e6,mat[,i],pch=19,cex=0.5)
    }
  }
  for(colony in cnsamples){
    if(type=="l"){
      lines(details$Pos/1e6,mat[,colony],lwd=0.5,col="red")
    }else{
      points(details$Pos/1e6,mat[,colony],pch=19,cex=0.5,col="red")
    }
  }
}

cap=function(x,x.min,x.max){
  ifelse(x>x.max,x,ifelse(x<x.min,x.min,x))
}
quickCopyNumber=function(PD,colonies=PD$pdx$meta$clones_short,
                         mutcolonies=c(),chroms=c(1:19,"X"),NSNP=20){
  dat=PD$pdx$dat
  dat$depth=dat$dep
  ##browser()
  #dat$cnsamples="PD45763b"
  dat$cnsamples=mutcolonies
  par(mfrow=c(2,2),oma = c(2, 2, 4, 2),mar=c(2,2,2,2))
  for(CHR in chroms){
    if(CHR=="Y"){
      #browser()
    }
    idx1=which(dat$details$Chrom==CHR)
    depth1=as.matrix(dat$depth[idx1,colonies])
    details1=dat$details[idx1,]
    
    #plot_matrix(depth1,details1,"Depth(Raw) Chr1",type = "p")
    
    ##NSNP=200 ## Length of rolling window
    
    rolldepth=do.call("cbind",frollmean(as.data.frame(depth1),NSNP,align="center",na.rm = TRUE));
    colnames(rolldepth)=colnames(depth1)
    
    plot_matrix(rolldepth,details1,"Rolling Depth",cnsamples = dat$cnsamples);title("Raw Rolling Depth")
    
    depthauto=dat$depth[which(dat$details$Chrom %in% 1:22),colonies]
    adjrolldepth=sapply(1:dim(depth1)[2],function(i) cap(log2(rolldepth[,i]/mean(depthauto[,i])),-3,3))
    colnames(adjrolldepth)=colnames(depth1)
    plot_matrix(adjrolldepth,details1,"Adjusted Rolling Log2-Depth",cnsamples = dat$cnsamples);title("Adjusted Rolling Log2(Depth)")
    xcentered=t(apply(adjrolldepth,1,function(x) (x-median(x))))
    colnames(xcentered)=colnames(depth1)
    plot_matrix(xcentered,details1,"Centred Log2(Depth)",cnsamples = dat$cnsamples);title("Cross-sectionally centred estimate of logR ")
    ##Robust Cross-sectional Z-Score.
    xrobustz=t(apply(adjrolldepth,1,function(x) (x-median(x))/mad(x)))
    colnames(xrobustz)=colnames(depth1)
    plot_matrix(xrobustz,details1,"Robust Z-score",cnsamples = dat$cnsamples);title("Cross-sectional Z-score")
    
    
    abline(h=4.5,lty="dotted")
    abline(h=-4.5,lty="dotted")
    mtext(text = sprintf("Chromosome %s: Rolling averages taken over %d SNPS",CHR,NSNP),side = 3,line = 1,outer = TRUE)
    res=apply(xrobustz,2,function(x) c(mean(x,na.rm=TRUE),max(abs(x),na.rm=TRUE)))
    res=t(res)
    #
    colnames(res)=c("mean","maxabs")
    print(cbind(res[order(-res[,2]),]) %>% as.data.frame %>% mutate(chr=CHR)  %>% filter(maxabs>4))
  }
}
require("VariantAnnotation")
export_vcf=function(PD){
  format="REF_COUNT:ALT_COUNT:GENOTYPE"
  header=readLines("../config/VCFHEADER.txt")
  ## Add in samples
  samples=PD$pdx$meta$clones
  clones=PD$pdx$meta$clones_short
  sample_template="##SAMPLE=<ID=%s,Description=.,Accession=.,Platform=hiseq,Protocol=WGS,SampleName=%s,Source=.>"
  sampleLines=sprintf(sample_template,samples,samples)
  mtr=as.matrix(PD$pdx$dat$mtr[,clones])
  wtr=as.matrix(PD$pdx$dat$wtr[,clones])
  geno=as.matrix(PD$pdx$dat$geno_colony[,clones])
  geno=ifelse(is.na(geno),"?",as.character(geno))
  idx=1:length(mtr)
  entries=data.frame(matrix(sprintf("%d:%d:%s",wtr[idx],mtr[idx],geno),ncol=length(clones)),stringsAsFactors = FALSE)
  colnames(entries)=samples
  #info SOMATIC;MUT_TYPE=;GENE = ;PROTEIN_CHANGE=
  N=dim(PD$pdx$dat$details)[1]
  fix=function(x,default="-"){ifelse(is.na(x),"-",x)}
  INFO=with(PD$pdx$dat$details,sprintf("SOMATIC;MUT_TYPE=%s;GENE=%s;PROTEIN_CHANGE=%s;VC=%s;IS_IN_EXCLUDED_REGION=%d",TYPE,fix(GENE),fix(HGVS_PROTEIN),fix(VC),is_localx_excluded))
  info=with(PD$pdx$dat$details,data.frame(CHROM=Chrom,POS=Pos,ID=sprintf("%s_%d",PD$patient,1:N),REF=Ref,ALT=Alt,QUAL="99",FILTER="PASS",INFO=INFO,FORMAT=format,stringsAsFactors = FALSE))
  chrom=c(1:22,"X","Y")
  idx=order(match(info$CHROM,chrom),info$POS)
  out=cbind(info,entries)[idx,]
  rownames(out)=NULL
  headerLine=paste0("#",paste(colnames(out),collapse="\t"))
  c(header,sampleLines,headerLine,apply(out,1,paste,collapse="\t"))
}

export_all=function(PDD,stub="./exports"){
  patients=names(PDD)
  system(sprintf("mkdir -p %s/VCF",stub))
  for(patient in patients){
    cat("writing VCF for ",patient,"\n")
    tmp=export_vcf(PDD[[patient]])
    fname=sprintf("%s/VCF/%s.vcf",stub,patient)
    writeLines(tmp,fname)
    ## Check recognisied as VCF
    vcf=readVcf(fname)
  }
}


