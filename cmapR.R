#Rewrite of cmap score for publication
#NWK
#2/7/13
library(sqldf)

rank.matrix <- cmap.rank
up.genes <- tsig.up
down.genes <- tsig.down
saps.gsea <- function(rank.matrix,up.genes,down.genes){
  #Performs GSEA on a given rank.matrix, provided a vector of genes/probes up and genes/probes down
  up.genes <- match(up.genes,rownames(rank.matrix))
  down.genes <- match(down.genes,rownames(rank.matrix))
  rank.up <- rank.matrix[up.genes,]
  rank.down <- rank.matrix[down.genes,]
  samples <- ncol(rank.matrix)
  
  t.up <- nrow(trank.up)
  t.down <- nrow(trank.down)
  n <- nrow(rank.matrix)
  
  vu <- apply(rank.up,2,function(x)x[order(x)])
  vd <- apply(rank.down,2,function(x)x[order(x)])
  as.u<- (seq(nrow(vu))/t.up)-(vu/n)
  bs.u <- (vu/n)-((seq(nrow(vu))-1)/t.up)
  as.d <- (seq(nrow(vd))/t.up)-(vd/n)
  bs.d <- (vd/n)-((seq(nrow(vd))-1)/t.down)
  as.up <- apply(as.u,2,max)
  bs.up <- apply(bs.u,2,max)
  
  ksu <- ifelse(as.up>bs.up,as.up,-bs.up)
  ksd <- ifelse(as.down>bs.down,as.down,-bs.down)  
  
  nS <- numeric(ncol(trank.up))
  ns <- ksu-ksd
  np=max(ns)
  nq=min(ns)
  
  kss <- sign(ksu)!=sign(ksd)
  nS[kss&(ns>0)] <- (ns/np)[kss&(ns>0)]
  nS[kss&(ns<0)] <- -(ns/nq)[kss&(ns<0)]
 
}

