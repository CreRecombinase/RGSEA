#Rewrite of GSEA code for SAPS package
#NWK
#2/7/13
library(sqldf)

cmap.rank <- read.csv.sql("C:/Users/nknoblau/Documents/R_WS/cmap_cel/nRankMat.txt",sep="\t",header=T,stringsAsFactors=F,eol="\n")
rownames(cmap.rank)<- cmap.rank[,1]
cmap.rank <- cmap.rank[,-1]


instance.names <- read.csv.sql("C:/Users/nknoblau/Documents/R_WS/cmap_cel/cmap_instances_02.txt",sep="\t",header=T,stringsAsFactors=F)


instance.names.o <- instance.names[order(instance.names$instance_id),]
rownames(instance.names.o)<-instance.names.o$instance_id

signatures <- scan("C:/Users/nknoblau/Documents/R_WS/cmap_cel/msigdb_gene_sets/msigdb_up_mapped_to_HG_U133A.gmt",what="character",sep="\n")

signatures <- strsplit(signatures,"\t")



signames <- sapply(signatures,head,n=1)

signatures <- sapply(signatures,tail,n=-2)
signatures.up <- signatures
names(signatures.up)<-signames

signatures.down <- scan("C:/Users/nknoblau/Documents/R_WS/cmap_cel/msigdb_gene_sets/msigdb_dn_mapped_to_HG_U133A.gmt",what="character",sep="\n")

signatures.down <- strsplit(signatures.down,"\t")

names(signatures.down) <- sapply(signatures.down,head,n=1)

signatures.down <- sapply(signatures.down,tail,n=-2)

tsig.up <- signatures.up[[1]]
tsig.down <- signatures.down[[1]]

trank.up <- cmap.rank[tsig.up,]

trank.down <- cmap.rank[tsig.down,]



as.up <- numeric(ncol(trank.up))
bs.up <- numeric(ncol(trank.up))
as.down <- numeric(ncol(trank.down))
bs.down <- numeric(ncol(trank.down))
t.up <- nrow(trank.up)
t.down <- nrow(trank.down)
n <- nrow(cmap.rank)
i <-1
for(i in 1:ncol(trank.up)){
  vju <- trank.up[order(trank.up[,i]),i]
  vjd <- trank.down[order(trank.down[,i]),i]
  tas.u <- (1:length(vju))/t.up-(vju/n)
  tbs.u <- (vju/n)-(((1:length(vju))-1)/t.up)
  tas.d <- (1:length(vjd))/t.down-(vjd/n)
  tbs.d <- (vjd/n)-(((1:length(vjd))-1)/t.down)
  as.up[i]<-max(tas.u)
  bs.up[i]<-max(tbs.u)
  as.down[i]<-max(tas.d)
  bs.down[i]<-max(tbs.d)
}


sum(as.up>bs.up)
ksu <- ifelse(as.up>bs.up,as.up,-bs.up)
ksd <- ifelse(as.down>bs.down,as.down,-bs.down)  


nS <- numeric(ncol(trank.up))
ns <- ksu-ksd
np=max(ns)
nq=min(ns)


kss <- sign(ksu)!=sign(ksd)
nS[kss&(ns>0)] <- (ns/np)[kss&(ns>0)]
nS[kss&(ns<0)] <- -(ns/nq)[kss&(ns<0)]

nsl <- split(nS,instance.names.o$cmap_name)

nslm <- sapply(nsl,mean)


cmaps <- split(x=nS,f=instance.names.o$cmap_name)

cmap.m <- sapply(cmaps,FUN=mean)
cmap.n <- sapply(cmaps,FUN=length)

newrank <- data.frame(instance.names.o$instance_id,instance.names.o$cmap_name,nS,ksu,ksd)
rownames(newrank)<-newrank$instance.names.o.instance_id


actual.rank <- read.table("C:/Users/nknoblau/Documents/R_WS/cmap_cel/detailedResults81929.txt",header=T,sep="\t",stringsAsFactors=T,quote="",row.names=1)
rownames(actual.rank)<-actual.rank$instance_id
actual.rank <- actual.rank[rownames(newrank),]
plot(newrank$nS,actual.rank$score)

colnames(trank.up) <- newrank$instance.names.o.instance_id
##Lets work backwards to calcualte the ksu and ksd for actual.rank[3,]
actual.rank[3,]


ttrank.u <- trank.up[,3]
ttrank.d <- trank.down[,3]

tta.u <- numeric(length(ttrank.u))
ttb.u <-  numeric(length(ttrank.u))
tta.d <- numeric(length(ttrank.d))
ttb.d <- numeric(length(ttrank.d))
tt.ut <- length(ttrank.u)
tt.dt <- length(ttrank.d)


ta.u <- max(tta.u)
tb.u <- max(ttb.u)
ta.u>tb.u
ttku <- -tb.u
actual.rank[3,]
newrank[3,]
ttku

for(i in 1:length(tta.d)){
  tta.d <- (i/tt.dt)-(ttrank.d[i]/n)
  ttb.d <- (ttrank.d[i]/n)-((i-1)/tt.dt)
}


actual.rank[3,]
###Method 2 of caluclating ksu and ksd (using ks.test)

rksd <- numeric(ncol(trank.down))

for(i in 1:length(rksd)){
  rksd[i]<-ks.test(trank.down[,i],cmap.rank[,i],alternative="less")[[1]]
}
rksd <- -rksd

rksu <- numeric(ncol(trank.up))
for(i in 1:length(rksd)){
  rksu[i] <- ks.test(trank.up[,i],cmap.rank[,i],alternative="less")[[1]]
}



rownames(trank.up)[1]
trank.up["204286_s_at","1171"]

rownames(trank.up)==tsig.up
tsig.up[1]


actual.rank$up[1]

ks.test(trank.up[,6],cmap.rank[,6],alternative="g")


actual.rank$up[5]

summary(rksu)
summary(actual.rank$up)
?ks.test
round(rksd[1],3)
plot(rksu[actual.rank$up<0],actual.rank$up[actual.rank$up<0])


prksu <- rksu[actual.rank$up>0]
paru <- actual.rank$up[actual.rank$up>0]

summary(paru)
summary(prksu)
sum(actual.rank$down>0)
sum(actual.rank$up>0)
sum(actual.rank$down>0&actual.rank$up>0)
prksd <- rksd[actual.rank$down>0]
pard <- actual.rank$down[actual.rank$down>0] 

plot(prksd,pard)
summary(prksd)

plot(prksu,paru)

summary(rksd)
summary(actual.rank$down)
length(rksd)
rksd[1]
actual.rank$down[1]











