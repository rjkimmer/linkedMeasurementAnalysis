library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gplots)
library(xlsx)
library('fgsea')

# Import raw data for GBM treated and untreated conditions 
tx.raw <- read.table('GBM_treated_rsem.txt', sep='\t', header=TRUE, row.names = 1)
tx.meta <- read.table('GBM_treated_stats.txt', sep = '\t', header = TRUE, row.names = 1)
utx.raw <- read.table('GBM_untreated_rsem.txt', sep ='\t', header = TRUE, row.names=1)
utx.meta <- read.table('GBM_untreated_stats.txt', sep = '\t', header=TRUE, row.names = 1)

full.dat <- cbind(tx.raw,utx.raw)
full.meta <- rbind(tx.meta[,1:12],utx.meta)

# Set up seurat object
gbm <- CreateSeuratObject(raw.data = full.dat,'gbm',min.cells = 5, min.genes = 1000, names.field = 1, is.expr = 0, meta.data = full.meta,do.scale = FALSE)
gbm <- NormalizeData(gbm)
gbm = FindVariableGenes(object = gbm,mean.function = ExpMean,dispersion.function=LogVMR, x.low.cutoff=0.5,y.cutoff=0.5,do.plot=TRUE)
gbm = ScaleData(gbm)

# Pare cells based on biophysical data -- (negative mass or MAR used to indicate cells that do not have matched biophysical properties)
tx.idx <- which(as.matrix(gbm@ident)=='GBMT' & gbm@meta.data$mass>0)
utx.idx <- which(as.matrix(gbm@ident)=='GBMUT' & gbm@meta.data$mass>0)
all.idx <- c(tx.idx, utx.idx)
# Include normalized MAR for single cells
gbm@meta.data$norm <- gbm@meta.data$mar/gbm@meta.data$mass


# Fig. 4b --  Biophysical response
par(pty = 's')
plot(gbm@meta.data$mass[utx.idx], gbm@meta.data$mar[utx.idx], pch = 21, bg = '#1EBDC2', cex = 1.5, xlab = 'Buoyant mass (pg)', ylab = 'MAR (pg/h)', cex.lab = 1.3, cex.axis = 1.25)
points(gbm@meta.data$mass[tx.idx], gbm@meta.data$mar[tx.idx], pch = 24, bg = '#F3766E', cex = 1.5 )
legend(x= 'topleft', legend = c('Treated', 'DMSO'), pch = c(24,21), pt.bg = c('#F3766E','#1EBDC2'), pt.cex = 1.5, box.lwd = 0)

# Supp. Fig 8 -- Normalized MAR boxplot and jitter plot
par(pty = 's')
x1 = rep(1,length(utx.idx))
x2 = rep(2,length(tx.idx))
par(pty = "s")
ddf <- data.frame(norm = gbm@meta.data$norm[c(utx.idx,tx.idx)], id = gbm@ident[c(utx.idx, tx.idx)])
boxplot(norm~id, data = ddf, outline = FALSE, lwd = 2, ylim = c(-0.025,0.075), staplewex = 0.1)
plot(jitter(x1,10),gbm@meta.data$norm[utx.idx],xlim= c(.5,2.5), ylim = c(-.025,0.075),bg='#1EBDC2',pch=21,cex = 1.75, ylab = "Normalized MAR (1/h)",cex.lab = 1.5, cex.axis=1.5)
points(jitter(x2,5),gbm@meta.data$norm[tx.idx],pch=24,cex=1.75,bg ='#F3766E')


# Differentially expressed genes between DMSO and RG7388 treated cells -- Supp. Table 8  
markers = FindMarkers(gbm,ident.1 = 'GBMT', ident.2 = 'GBMUT', logfc.threshold = 0)

# Mass correlation with gene expression 
# Determine mass correlation for the untreated cells
mass.corr.gbm.utx = apply(t(as.matrix(gbm@data[,utx.idx])),2,function(x) cor(x,gbm@meta.data$mass[utx.idx],method='spearman'))
mass.corr.gbm.utx[is.na(mass.corr.gbm.utx)]=0


# Figure 4d -- Mitosis scoring
mass.corr.gbm.utx = apply(t(as.matrix(gbm@data[,utx.idx])),2,function(x) cor(x,gbm@meta.data$mass[utx.idx],method='spearman'))
mass.corr.gbm.utx[is.na(mass.corr.gbm.utx)]=0
mitotic.genes <- as.matrix(read.table('mitosis.txt', sep = '\t', header = TRUE))
mitotic.corr <- mass.corr.gbm.utx[mitotic.genes]
mitotic.corr.sig <- mitotic.corr[which(mitotic.corr>0.2147)] #2*sd of null distribution based on untreated mass data

# Scatter between mito score and cell mass
utx.mito <- gbm@scale.data[names(mitotic.corr.sig), utx.idx]
utx.mito.score <- colMeans(utx.mito)

tx.mito <- gbm@scale.data[names(mitotic.corr.sig), tx.idx]
tx.mito.score <- colMeans(tx.mito)

par(pty = 's')
plot(gbm@meta.data$mass[utx.idx], utx.mito.score, pch = 21, cex = 1.5, bg = '#1EBDC2', xlab = 'Buoyant mass (pg)', ylab = 'Mitotic score', cex.lab=1.3, cex.axis = 1.2 )
points(gbm@meta.data$mass[tx.idx], tx.mito.score, pch = 24, cex = 1.5, bg ='#F3766E')

## Supp. Tables. 1/2
mass.corr.gbm.tx = apply(t(as.matrix(gbm@data[,tx.idx])),2,function(x) cor(x,gbm@meta.data$mass[tx.idx],method='spearman'))
mass.corr.gbm.tx[is.na(mass.corr.gbm.tx)]=0
mass.corr.gbm.tx.rank = sort(mass.corr.gbm.tx, decreasing=TRUE,index.return=TRUE)$ix
mass.ranked.list = mass.corr.gbm.tx[mass.corr.gbm.tx.rank]
names(mass.ranked.list) = toupper(names(mass.corr.gbm.tx)[mass.corr.gbm.tx.rank])
mass.ranked.enrichment = fgsea(go.set,mass.ranked.list,nperm=15000,maxSize = 500)

mass.corr.gbm.utx = apply(t(as.matrix(gbm@data[,utx.idx])),2,function(x) cor(x,gbm@meta.data$mass[utx.idx],method='spearman'))
mass.corr.gbm.utx[is.na(mass.corr.gbm.utx)]=0
mass.corr.gbm.utx.rank = sort(mass.corr.gbm.utx, decreasing=TRUE,index.return=TRUE)$ix
mass.ranked.list = mass.corr.gbm.utx[mass.corr.gbm.utx.rank]
names(mass.ranked.list) = toupper(names(mass.corr.gbm.utx)[mass.corr.gbm.utx.rank])
mass.ranked.enrichment = fgsea(go.set,mass.ranked.list,nperm=15000,maxSize = 1500)

norm.corr.gbm.tx = apply(t(as.matrix(gbm@data[,tx.idx])),2,function(x) cor(x,gbm@meta.data$norm[tx.idx],method='spearman'))
norm.corr.gbm.tx[is.na(norm.corr.gbm.tx)]=0
norm.corr.gbm.tx.rank = sort(norm.corr.gbm.tx, decreasing=TRUE,index.return=TRUE)$ix
norm.ranked.list = norm.corr.gbm.tx[norm.corr.gbm.tx.rank]
names(norm.ranked.list) = toupper(names(norm.corr.gbm.tx)[norm.corr.gbm.tx.rank])
norm.ranked.enrichment = fgsea(go.set,norm.ranked.list,nperm=15000,maxSize = 500)

norm.corr.gbm.utx = apply(t(as.matrix(gbm@data[,utx.idx])),2,function(x) cor(x,gbm@meta.data$norm[utx.idx],method='spearman'))
norm.corr.gbm.utx[is.na(norm.corr.gbm.utx)]=0
norm.corr.gbm.utx.rank = sort(norm.corr.gbm.utx, decreasing=TRUE,index.return=TRUE)$ix
norm.ranked.list = norm.corr.gbm.utx[norm.corr.gbm.utx.rank]
names(norm.ranked.list) = toupper(names(norm.corr.gbm.utx)[norm.corr.gbm.utx.rank])
norm.ranked.enrichment = fgsea(go.set,norm.ranked.list,nperm=15000,maxSize = 500)

# Supp. Fig. 3/4
# Treated -- normalized MAR
gbm.tx.norm.null <- sapply(1:10, function(x) sample(tx.idx, length(tx.idx), replace =FALSE))
gbm.tx.norm.corr.null <- sapply(1:10, function(x) (apply(gbm@data[,tx.idx], 1, function(y) cor(y,gbm@meta.data$norm[gbm.tx.norm.null[,x]], method = 'spearman'))))
gbm.tx.norm.corr.null[is.na(gbm.tx.norm.corr.null)]<-0
gbm.tx.norm.null.mean <- mean(colMeans(gbm.tx.norm.corr.null))
gbm.tx.norm.null.sd <- mean(apply(gbm.tx.norm.corr.null,2, function(x) sd(x)))
norm.corr.gbm.tx = apply(t(as.matrix(gbm@data[,tx.idx])),2,function(x) cor(x,gbm@meta.data$norm[tx.idx],method='spearman'))
norm.corr.gbm.tx[is.na(norm.corr.gbm.tx)]=0
norm.corr.gbm.tx.rank = sort(norm.corr.gbm.tx, decreasing=TRUE,index.return=TRUE)$ix

cols<- c("blue","red")[(abs(norm.corr.gbm.tx[norm.corr.gbm.tx.rank])>2*gbm.tx.norm.null.sd) + 1]
barplot(norm.corr.gbm.tx[norm.corr.gbm.tx.rank], col=cols,border=NA,ylim = c(-0.75,.75),xaxt='n', main = c('GBM tx norm',as.character(length(which(norm.corr.gbm.tx>2*gbm.tx.norm.null.sd))),as.character(length(which(norm.corr.gbm.tx< -2*gbm.tx.norm.null.sd)))))
box()
abline(h = 2*gbm.tx.norm.null.sd,lwd=2,col=1,lty=2)
abline(h = -2*gbm.tx.norm.null.sd,lwd=2,col=1,lty=2)
abline(h=0,lwd=1)

# untreated -- normalized MAR
gbm.utx.norm.null <- sapply(1:10, function(x) sample(utx.idx, length(utx.idx), replace =FALSE))
gbm.utx.norm.corr.null <- sapply(1:10, function(x) (apply(gbm@data[,utx.idx], 1, function(y) cor(y,gbm@meta.data$norm[gbm.utx.norm.null[,x]], method = 'spearman'))))
gbm.utx.norm.corr.null[is.na(gbm.utx.norm.corr.null)]<-0
gbm.utx.norm.null.mean <- mean(colMeans(gbm.utx.norm.corr.null))
gbm.utx.norm.null.sd <- mean(apply(gbm.utx.norm.corr.null,2, function(x) sd(x)))
norm.corr.gbm.utx = apply(t(as.matrix(gbm@data[,utx.idx])),2,function(x) cor(x,gbm@meta.data$norm[utx.idx],method='spearman'))
norm.corr.gbm.utx[is.na(norm.corr.gbm.utx)]=0
norm.corr.gbm.utx.rank = sort(norm.corr.gbm.utx, decreasing=TRUE,index.return=TRUE)$ix

cols<- c("blue","red")[(abs(norm.corr.gbm.utx[norm.corr.gbm.utx.rank])>2*gbm.utx.norm.null.sd) + 1]
barplot(norm.corr.gbm.utx[norm.corr.gbm.utx.rank], col=cols,border=NA,ylim = c(-0.75,.75),xaxt='n', main = c('GBM utx norm',as.character(length(which(norm.corr.gbm.utx>2*gbm.utx.norm.null.sd))),as.character(length(which(norm.corr.gbm.utx< -2*gbm.utx.norm.null.sd)))))
box()
abline(h = 2*gbm.utx.norm.null.sd,lwd=2,col=1,lty=2)
abline(h = -2*gbm.utx.norm.null.sd,lwd=2,col=1,lty=2)
abline(h=0,lwd=1)

# Treated -- mass
gbm.tx.mass.null <- sapply(1:10, function(x) sample(tx.idx, length(tx.idx), replace =FALSE))
gbm.tx.mass.corr.null <- sapply(1:10, function(x) (apply(gbm@data[,tx.idx], 1, function(y) cor(y,gbm@meta.data$mass[gbm.tx.mass.null[,x]], method = 'spearman'))))
gbm.tx.mass.corr.null[is.na(gbm.tx.mass.corr.null)]<-0
gbm.tx.mass.null.mean <- mean(colMeans(gbm.tx.mass.corr.null))
gbm.tx.mass.null.sd <- mean(apply(gbm.tx.mass.corr.null,2, function(x) sd(x)))
mass.corr.gbm.tx = apply(t(as.matrix(gbm@data[,tx.idx])),2,function(x) cor(x,gbm@meta.data$mass[tx.idx],method='spearman'))
mass.corr.gbm.tx[is.na(mass.corr.gbm.tx)]=0
mass.corr.gbm.tx.rank = sort(mass.corr.gbm.tx, decreasing=TRUE,index.return=TRUE)$ix

cols<- c("blue","red")[(abs(mass.corr.gbm.tx[mass.corr.gbm.tx.rank])>2*gbm.tx.mass.null.sd) + 1]
barplot(mass.corr.gbm.tx[mass.corr.gbm.tx.rank], col=cols,border=NA,ylim = c(-0.75,.75),xaxt='n', main = c('GBM tx mass',as.character(length(which(mass.corr.gbm.tx>2*gbm.tx.mass.null.sd))),as.character(length(which(mass.corr.gbm.tx< -2*gbm.tx.mass.null.sd)))))
box()
abline(h = 2*gbm.tx.mass.null.sd,lwd=2,col=1,lty=2)
abline(h = -2*gbm.tx.mass.null.sd,lwd=2,col=1,lty=2)
abline(h=0,lwd=1)

# Untreated -- mass
gbm.utx.mass.null <- sapply(1:10, function(x) sample(utx.idx, length(utx.idx), replace =FALSE))
gbm.utx.mass.corr.null <- sapply(1:10, function(x) (apply(gbm@data[,utx.idx], 1, function(y) cor(y,gbm@meta.data$mass[gbm.utx.mass.null[,x]], method = 'spearman'))))
gbm.utx.mass.corr.null[is.na(gbm.utx.mass.corr.null)]<-0
gbm.utx.mass.null.mean <- mean(colMeans(gbm.utx.mass.corr.null))
gbm.utx.mass.null.sd <- mean(apply(gbm.utx.mass.corr.null,2, function(x) sd(x)))
mass.corr.gbm.utx = apply(t(as.matrix(gbm@data[,utx.idx])),2,function(x) cor(x,gbm@meta.data$mass[utx.idx],method='spearman'))
mass.corr.gbm.utx[is.na(mass.corr.gbm.utx)]=0
mass.corr.gbm.utx.rank = sort(mass.corr.gbm.utx, decreasing=TRUE,index.return=TRUE)$ix

cols<- c("blue","red")[(abs(mass.corr.gbm.utx[mass.corr.gbm.utx.rank])>2*gbm.utx.mass.null.sd) + 1]
barplot(mass.corr.gbm.utx[mass.corr.gbm.utx.rank], col=cols,border=NA,ylim = c(-0.75,.75),xaxt='n', main = c('GBM utx mass',as.character(length(which(mass.corr.gbm.utx>2*gbm.utx.mass.null.sd))),as.character(length(which(mass.corr.gbm.utx< -2*gbm.utx.mass.null.sd)))))
box()
abline(h = 2*gbm.utx.mass.null.sd,lwd=2,col=1,lty=2)
abline(h = -2*gbm.utx.mass.null.sd,lwd=2,col=1,lty=2)
abline(h=0,lwd=1)

## Compare spearman and pearson coefficients for mass/mar treated/untreated
mass.corr.gbm.utx.spearman = apply(t(as.matrix(gbm@data[,utx.idx])),2,function(x) cor(x,gbm@meta.data$mass[utx.idx],method='spearman'))
mass.corr.gbm.utx.spearman[is.na(mass.corr.gbm.utx.spearman)] <- 0
mass.corr.gbm.utx.pearson = apply(t(as.matrix(gbm@data[,utx.idx])),2,function(x) cor(x,gbm@meta.data$mass[utx.idx],method='pearson'))
mass.corr.gbm.utx.pearson[is.na(mass.corr.gbm.utx.pearson)] <- 0

norm.corr.gbm.utx.spearman = apply(t(as.matrix(gbm@data[,utx.idx])),2,function(x) cor(x,gbm@meta.data$norm[utx.idx],method='spearman'))
norm.corr.gbm.utx.spearman[is.na(norm.corr.gbm.utx.spearman)] <- 0
norm.corr.gbm.utx.pearson = apply(t(as.matrix(gbm@data[,utx.idx])),2,function(x) cor(x,gbm@meta.data$norm[utx.idx],method='pearson'))
norm.corr.gbm.utx.pearson[is.na(norm.corr.gbm.utx.pearson)] <- 0

mass.corr.gbm.tx.spearman = apply(t(as.matrix(gbm@data[,tx.idx])),2,function(x) cor(x,gbm@meta.data$mass[tx.idx],method='spearman'))
mass.corr.gbm.tx.spearman[is.na(mass.corr.gbm.tx.spearman)] <- 0
mass.corr.gbm.tx.pearson = apply(t(as.matrix(gbm@data[,tx.idx])),2,function(x) cor(x,gbm@meta.data$mass[tx.idx],method='pearson'))
mass.corr.gbm.tx.pearson[is.na(mass.corr.gbm.tx.pearson)] <- 0

norm.corr.gbm.tx.spearman = apply(t(as.matrix(gbm@data[,tx.idx])),2,function(x) cor(x,gbm@meta.data$norm[tx.idx],method='spearman'))
norm.corr.gbm.tx.spearman[is.na(norm.corr.gbm.tx.spearman)] <- 0
norm.corr.gbm.tx.pearson = apply(t(as.matrix(gbm@data[,tx.idx])),2,function(x) cor(x,gbm@meta.data$norm[tx.idx],method='pearson'))
norm.corr.gbm.tx.pearson[is.na(norm.corr.gbm.tx.pearson)] <- 0

# mass utx
par(pty = 's')
plot(mass.corr.gbm.utx.spearman, mass.corr.gbm.utx.pearson, pch = 21, bg = rgb(0,0,1,0.3), xlab = 'spearman Coefficient', ylab = 'Pearson Coefficeint', xlim = c(-0.8,0.8), ylim = c(-0.8,.8), cex.lab =1.2, cex.axis = 1.2,cex=0.75,col=rgb(0,0,0,0))
title(c('mass utx', as.character(cor(mass.corr.gbm.utx.spearman, mass.corr.gbm.utx.pearson))))

# mass tx
par(pty = 's')
plot(mass.corr.gbm.tx.spearman, mass.corr.gbm.tx.pearson, pch = 21, bg = rgb(0,0,1,0.3), xlab = 'spearman Coefficient', ylab = 'Pearson Coefficeint', xlim = c(-0.8,0.8), ylim = c(-0.8,.8), cex.lab =1.2, cex.axis = 1.2,cex=0.75,col=rgb(0,0,0,0))
title(c('mass tx', as.character(cor(mass.corr.gbm.tx.spearman, mass.corr.gbm.tx.pearson))))

# norm utx
par(pty = 's')
plot(norm.corr.gbm.utx.spearman, norm.corr.gbm.utx.pearson, pch = 21, bg = rgb(0,0,1,0.3), xlab = 'spearman Coefficient', ylab = 'Pearson Coefficeint', xlim = c(-0.8,0.8), ylim = c(-0.8,.8), cex.lab =1.2, cex.axis = 1.2,cex=0.75,col=rgb(0,0,0,0))
title(c('norm utx', as.character(cor(norm.corr.gbm.utx.spearman, norm.corr.gbm.utx.pearson))))

# norm tx
par(pty = 's')
plot(norm.corr.gbm.tx.spearman, norm.corr.gbm.tx.pearson, pch = 21, bg = rgb(0,0,1,0.3), xlab = 'spearman Coefficient', ylab = 'Pearson Coefficeint', xlim = c(-0.8,0.8), ylim = c(-0.8,.8), cex.lab =1.2, cex.axis = 1.2,cex=0.75,col=rgb(0,0,0,0))
title(c('norm tx', as.character(cor(norm.corr.gbm.tx.spearman, norm.corr.gbm.tx.pearson))))


## Supp. Fig. 9
gbm1 <- gbm

####1: TSNEs of TREATED AND UNTREATED CELLS#####
gbm1<-SetAllIdent(gbm1,id = "orig.ident")
gbm1.t<-SubsetData(gbm1,ident.use = "GBMT")
gbm1.ut<-SubsetData(gbm1,ident.use = "GBMUT")

# Treated 
var.genes.tx<-gbm1.t@var.genes[-c(grep("MT-",gbm1.t@var.genes),
                                                    grep("MTR",gbm1.t@var.genes),
                                                    grep("RPL",gbm1.t@var.genes),
                                                    grep("RPS",gbm1.t@var.genes),
                                                    grep("RNA",gbm1.t@var.genes))]
var.genes.tx
gbm1.t<-RunPCA(gbm1.t,pc.genes = var.genes.tx)
gbm1.t<-SetAllIdent(gbm1.t,id = "orig.ident")
gbm1.t<-RunTSNE(gbm1.t,dims.use=c(1,2,3,4),dim_embed=2,max_iter = 2500,perplexity=15,do.fast = T)
TSNEPlot(gbm1.t,dim.1=1,dim.2=2,group.by = "orig.ident",pt.size = 2)

# Untreated
var.genes.utx<-gbm1.ut@var.genes[-c(grep("MT-",gbm1.ut@var.genes),
                                  grep("MTR",gbm1.ut@var.genes),
                                  grep("RPL",gbm1.ut@var.genes),
                                  grep("RPS",gbm1.ut@var.genes),
                                  grep("RNA",gbm1.ut@var.genes))]
var.genes.utx
gbm1.ut<-RunPCA(gbm1.ut,pc.genes = var.genes.utx)
gbm1.ut<-SetAllIdent(gbm1.ut,id = "orig.ident")
gbm1.ut<-RunTSNE(gbm1.ut,dims.use=c(1,2,3,4),dim_embed=2,max_iter = 2500,perplexity=15,do.fast = T)
TSNEPlot(gbm1.ut,dim.1=1,dim.2=2,group.by = "orig.ident",pt.size = 2)


# Both together
var.genes<-gbm1@var.genes[-c(grep("MT-",gbm1@var.genes),
                                    grep("MTR",gbm1@var.genes),
                                    grep("RPL",gbm1@var.genes),
                                    grep("RPS",gbm1@var.genes),
                                    grep("RNA",gbm1@var.genes))]
var.genes
gbm1<-RunPCA(gbm1,pc.genes = var.genes.utx)
gbm1<-SetAllIdent(gbm1,id = "orig.ident")
gbm1<-RunTSNE(gbm1,dims.use=c(1,2,3,4),dim_embed=2,max_iter = 2500,perplexity=15,do.fast = T)
TSNEPlot(gbm1,dim.1=1,dim.2=2,group.by = "orig.ident",pt.size = 2)