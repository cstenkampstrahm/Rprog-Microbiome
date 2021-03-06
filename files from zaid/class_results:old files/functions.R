### Prepared by Zaid Abdo/Metagenomics Class-2017
### Taxonomic analysis: pairing taxa with OTUs and recovering the higher 
##  classification levels if needed
taxa.split.combine.ftn = function(taxa, d.df, l = 7){
  ### usage taxa.split.combine.ftn(taxa = dataframe recovered from .taxonomy mothur, 
  ##                              d.df = dataframe recovered from .shared file after 
  ##                              removing first 3 columns and using 2 column to name 
  ##                              rows, l = taxonomic level [1=Domain, 2=Phylum,3=Class
  ##                              4=Order, 5=Family, 6=Genus, 7=otu/species])
  ### spliting on taxa level
  taxa.names = strsplit(as.character(taxa[,3]),split=";")
  taxa.vt = c()
  for(i in 1:length(taxa[,1])){
    temp = strsplit(taxa.names[[i]][1],split='\\(')
    taxa.vt[i] = temp[[1]][1]
  }
  ### Combining on family
  if(l==7){
    taxa.df = d.df
  }else{
    levels.vt = levels(as.factor(taxa.vt))
    taxa.mt = matrix(NA,nrow=length(d.df[,1]),ncol=length(levels.vt))
    for(i in 1:length(levels.vt)){
      col.vt = which(taxa.vt==levels.vt[i])
      if(length(col.vt)>1){
        taxa.mt[,i] = apply(d.df[,col.vt],1,sum)
      }else{
        taxa.mt[,i] = d.df[,col.vt]
      }
    }
    taxa.df = data.frame(taxa.mt,row.names = row.names(d.df))
    names(taxa.df) = levels.vt
    taxa.vt = levels.vt
    d.cs = apply(taxa.df,2,sum)
    taxa.df = taxa.df[,d.cs>0]
  }
  return(list(taxa.df,taxa.vt))
}

###### No experimental design
### species bar plot per sample using most abundant in sample
bar.taxa.sample.ftn = function(d.df,cutoff=0,xlab="Sample",llab="OTU"){
  d.rs = apply(d.df,1,sum)
  ln = length(d.rs)
  p.df = d.df/d.rs
  if(cutoff > 0){
    cut.b = apply(p.df,2,function(x)sum(x>cutoff))
    p.df = p.df[,cut.b > 0]
    p.rs = apply(p.df,1,sum)
    p.df = p.df/p.rs
  }
  nm = colnames(p.df)
  r.nm = row.names(p.df)
  sp = rep(nm,ln)
  sample = w = c()
  for(i in 1:ln){
    sample = c(sample,rep(r.nm[i],length(nm)))
    w = c(w,t(p.df[i,]))
  }
  hist.df = data.frame(sample,sp,w)
  print(qplot(data=hist.df,x=sample,fill=sp,geom="bar",weight=w)+
          labs(x=xlab,y="Proportion")+scale_fill_discrete(name=llab)+
          theme(axis.text.x=element_text(angle=90))+
          scale_color_manual(values = palette(rainbow(19))))
}

### species bar plot per sample using most abundant in total
bar.taxa.ftn = function(d.df,cutoff=0,xlab="Sample",llab="OTU"){
  d.cs = apply(d.df,2,sum)
  d.cp = d.cs/sum(d.cs)
  d.df = d.df[d.cp > cutoff]
  nm = names(d.df)
  r.nm = row.names(d.df)
  d.sc = apply(d.df,1,sum)
  p.df = as.matrix(d.df/d.sc)
  ln = length(r.nm)
  sp = rep(nm,ln)
  exp = w = c()
  for(i in 1:ln){
    exp = c(exp,rep(r.nm[i],length(nm)))
    w = c(w,p.df[i,])
  }
  hist.df = data.frame(exp,sp,w)
  print(qplot(data=hist.df,x=exp,fill=sp,geom="bar",weight=w)+labs(x=xlab,y="Proportion")+scale_fill_discrete(name=llab)+
          theme(axis.text.x=element_text(angle=90)))
  #scale_color_manual(values = palette(rainbow(19)))
}

### pheatmap no added treatments clustering function
pheatmap.1.ftn = function(d.df,cutoff=3,cellh=10,cellw=10,fonts=10){
  exp = row.names(d.df)
  exp1 = palette(rainbow(length(exp)))
  names(exp1) = exp1
  ann_col = list(exp1=exp1)
  
  print(pheatmap(as.matrix(d.df),
           cluster_rows = T,
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D",
           cluster_cols = F,
           cutree_rows = cutoff,
           cellwidth = cellw,
           cellheight = cellh,
           annotation_colors = ann_col,
           fontsize = fonts,
           height = 7,
           width = 6.5))
}

##### With experimental design
### species bar plot for an experiment
bar.ftn = function(d.df,Trt, title="Experiment",xlab="Sample",llab="OTU"){
  nm = colnames(d.df)
  ln = length(d.df[,1])
  p.df = d.df/apply(d.df,1,sum)
  hist.df = c()
  for(j in 1:length(p.df[1,])){
    sp = rep(nm[j],ln)
    w = p.df[,j]
    temp = data.frame(Trt,sp,w)
    hist.df = rbind(hist.df,temp)
  }
  if(length(hist.df[1,]) == 4){
    print(qplot(data=hist.df,x=hist.df[,2],fill=sp,geom="bar",weight=w,facets=~hist.df[,1])
          +labs(title=title,x=xlab,y="Proportion")+scale_fill_discrete(name=llab)
          +theme(axis.text.x=element_text(angle=90))
          +scale_color_manual(values = palette(rainbow(19))))
  }else{
    print(qplot(data=hist.df,x=hist.df[,3],fill=sp,geom="bar",weight=w,facets=hist.df[,1]~hist.df[,2])
          +labs(title=title,x=xlab,y="Proportion")+scale_fill_discrete(name=llab)
          +theme(axis.text.x=element_text(angle=90))
          +scale_color_manual(values = palette(rainbow(19))))
  }
}

### pheatmap no added treatments clustering function
pheatmap.2.ftn = function(d.df,Trt,cutoff=3,cellh=10,cellw=10,fonts=10){
  Trt = Trt[,1,drop=FALSE]
  exp1 = levels(Trt)
  exp2 = palette(rainbow(length(exp1)))
  names(exp2) = exp1
  ann_col = list(exp2=exp2)
  Trt = data.frame(Trt)
  
  print(pheatmap(as.matrix(d.df),
                 cluster_rows = T,
                 clustering_distance_rows = "euclidean",
                 clustering_method = "ward.D",
                 cluster_cols = F,
                 cutree_rows = cutoff,
                 cellwidth = cellw,
                 cellheight = cellh,
                 annotation_row = Trt[,1,drop=FALSE],
                 annotation_colors = ann_col,
                 fontsize = fonts,
                 height = 7,
                 width = 6.5))
}

### ordination plot
ord.plot.ftn <- function(d.df,Trt,title="Experiment",x=c(-1,1),y=c(-1,1),type=1,cex=0.5){
  l = levels(factor(Trt))
  c = seq(1,length(l))
  c1 = c()
  for(i in 1:length(l)){
    c1 = c(c1,rep(i,length(Trt[Trt==l[i]])))
  }
  s <- apply(d.df,2,sum)
  d.df <- d.df[,s>0]
  ord <- metaMDS(d.df,distance="bray",try=100,trymax=1000)
  plot(ord, type="n",xlim=x,ylim=y)
  if(type==1){
    ordiellipse(ord, Trt, col=c,lwd=1,label=TRUE)
  }else if(type == 2){
    text(ord,display="site",cex=cex,col=c1) 
  }else if(type == 3){
    ordispider(ord, Trt, col=c, cex = cex, label = TRUE)
  }else if(type == 4){
    text(ord,display="site",cex=cex,col=c1) 
    ordiellipse(ord, Trt, col=c,lwd=1,label=TRUE)
  }else if(type == 5){
    text(ord,display="site",cex=cex,col=c1) 
    ordispider(ord, Trt, col=c, cex = cex, label = TRUE)
  }else if(type==6){
    ordiellipse(ord, Trt, col=c,lwd=1,label=TRUE)
    ordispider(ord, Trt, col=c, cex = cex, label = FALSE)
  }else{
    text(ord,display="site",cex=cex,col=c1) 
    ordiellipse(ord, Trt, col=c,lwd=1,label=TRUE)
    ordispider(ord, Trt, col=c, cex = cex, label = FALSE)
  }
}

##### General function
### plotting a data frame
df.plot.ftn = function(d.df){
  rnm = row.names(d.df)
  nm = colnames(d.df)
  row.vt = col.vt = value.vt = c()
  for(i in 1:length(nm)){
    col.vt = c(col.vt,rep(nm[i],length(rnm)))
    row.vt = c(row.vt,rnm)
    value.vt = c(value.vt,d.df[,i])
  }
  result.df = cbind(col.vt,row.vt,value.vt)
  result.df = data.frame(result.df)
  print(qplot(data=result.df,x=factor(row.vt),y=as.numeric(value.vt),facets=~as.factor(col.vt))+theme(axis.text.x=element_text(angle=90)))
}
