#(Total) Rewrite of GSEA code for SAPS package
#NWK
#2/7/13


#number of samples
K= 1000
#number of genes
N=3000

phenotype <- runif(K,-1,1)*10

correlated <- rnorm(N)


D <- matrix(runif(N*K),N,K)
rownames(D)<- paste0("gene_",1:N)



D <- D+correlated%*%t(phenotype)

cors <- apply(D,MARGIN=1,function(x)cor(x,phenotype))
S <- c(rownames(D)[order(cors)][1:13],tail(rownames(D)[order(cors)],n=2))

L <- rank(cors)

S.cors <- L[S]
nS.cors <- L[!names(L)%in% S]


N.r  <- sum(abs(S.cors))


phit <- numeric(N)
pmiss <- numeric()

i <- 1
for(i in 1:length(phit)){
  if(names(L)[i]%in%S){
    phit[i]<-phit[max(i-1,1)]+(L[i]/N.r)
  }else{
    phit[i]<-phit[max(i-1,1)]
  }
}

pmiss <- numeric(N)
for(i in 1:length(pmiss)){
  if(!names(L)[i]%in%S){
    pmiss[i]<- pmiss[max(i-1,1)]+(1/(N-length(S)))
  }else{
    pmiss[i] <- pmiss[max(i-1,1)]
  } 
}



