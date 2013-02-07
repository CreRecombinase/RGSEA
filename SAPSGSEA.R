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
D <- D+correlated%*%t(phenotype)

g.cors <- rank(apply(D,MARGIN=1,function(x)cor(x,phenotype)))

#cutoff
i <- 200


