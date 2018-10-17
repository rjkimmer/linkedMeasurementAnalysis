library(Seurat)
library(ggplot2)
library(xlsx)
library(fgsea)

# Enrichment reference set 
go.set = gmtPathways('c5.all.v6.1.symbols.gmt')

## READ IN RAW DATA AND META DATA
raw.data = read.table("CD8_rsem.txt",sep = "\t", header = TRUE, row.names=1)
#raw.data=log(raw.data+1)
meta = read.table("cd8_qc2.txt",sep = "\t", header = TRUE, row.names = 1)

## SET UP SEURAT OBJECT
cd8 = CreateSeuratObject(raw.data=raw.data,project = "cd8",min.cells = 19,names.field=1,min.genes=1000,is.expr=0,meta.data=meta, do.scale=FALSE)
cd8 = NormalizeData(object= cd8, normalization.method="LogNormalize")
cd8 = FindVariableGenes(object = cd8,mean.function = ExpMean,dispersion.function=LogVMR, x.low.cutoff=0.5,y.cutoff=0,do.plot=TRUE)
cd8 = ScaleData(object=cd8)
cd8 = RunPCA(object=cd8,pcs.compute=15,do.print=TRUE,pcs.print=1:5,genes.print=1:5)

cd8.24 = which(cd8@meta.data$time==24 & cd8@meta.data$mar>-100)
cd8.48 = which(cd8@meta.data$time==48 & cd8@meta.data$mar>-100)
cd8.tot = c(cd8.24,cd8.48)

# Determine differentially expressed genes -- Supp. Tables 5/6
markers = FindMarkers(cd8,ident.1 = 'T24hr', ident.2 = 'T48hr')
# Sort based on average diff. exp. fold change
markers.idx = sort(markers[,2], decreasing=FALSE, index.return=TRUE)$ix
markers.sort = as.numeric(markers[markers.idx,2])
names(markers.sort) = toupper(rownames(markers)[markers.idx])
marker.ranked.enrichment = fgsea(go.set,markers.sort,nperm=15000,maxSize = 500)
marker.ranked.enrichment = marker.ranked.enrichment[sort(marker.ranked.enrichment$NES,index.return=TRUE, decreasing = FALSE)$ix,1:6]


# Figure 4a --> Mass and Mar for both time points
par(pty = "s")
plot(cd8@meta.data$mass[cd8.24],cd8@meta.data$mar[cd8.24], xlim = c(20,80), ylim = c(0,8.5), pch=21, cex=1.75, bg = rgb(0,0,1,1),xlab = "Mass (pg)",ylab = "MAR (pg/h)", cex.axis = 1.5, cex.lab=1.5)
points(cd8@meta.data$mass[cd8.48],cd8@meta.data$mar[cd8.48], xlim = c(20,80), ylim = c(0,8.5), pch=24, cex=1.75, bg = rgb(1,0,0,1))

# Figure 4b --> Plot normalized growth rates 
norm.24 = cd8@meta.data$mar[cd8.24]/cd8@meta.data$mass[cd8.24]
norm.48 = cd8@meta.data$mar[cd8.48]/cd8@meta.data$mass[cd8.48]
x1 = rep(1,length(norm.24))
x2 = rep(2,length(norm.48))
par(pty = "s")
plot(jitter(x1,10),norm.24,xlim= c(.5,2.5), ylim = c(0,0.15),bg=rgb(0,0,1,1),pch=21,cex = 1.75, ylab = "Normalized MAR (1/h)",cex.lab = 1.5, cex.axis=1.5)


# Generate supplementary gene lists and enrichment analyses
# Mass -- Pearson
mass.corr.pears = apply(t(as.matrix(cd8@data[,cd8.24])),2,function(x) cor(x,cd8@meta.data$mass[cd8.24]))
mass.corr.pears[is.na(mass.corr.pears)]=0


# Mass -- Spearman
mass.corr.spear.24 = apply(t(as.matrix(cd8@data[,cd8.24])),2,function(x) cor(x,cd8@meta.data$mass[cd8.24],method='spearman'))
mass.corr.spear.24[is.na(mass.corr.spear.24)]=0
mass.corr.spear.24.rank = sort(mass.corr.spear.24, decreasing=TRUE,index.return=TRUE)$ix
mass.ranked.list = mass.corr.spear.24[mass.corr.spear.24.rank]
names(mass.ranked.list) = toupper(names(mass.corr.spear.24)[mass.corr.spear.24.rank])
mass.ranked.enrichment = fgsea(go.set,mass.ranked.list,nperm=10000,maxSize = 500)
mass.ranked.enrichment <- mass.ranked.enrichment[mass.ranked.enrichment$padj< 0.1]
mass.ranked.enrichment$NES[is.na(mass.ranked.enrichment$NES)]=0

# Null Distribution for mass 
cd8.mass.null = sapply(1:10, function(x) sample(cd8.24,length(cd8.24),replace=FALSE))
cd8.mass.corr.null = sapply(1:10,function(x) (apply(cd8@data[,cd8.24],1,function(y) cor(y,cd8@meta.data$mass[cd8.mass.null[,x]],method='spearman'))))
cd8.mass.corr.null[is.na(cd8.mass.corr.null)]=0
cd8.mass.null.mean = mean(colMeans(cd8.mass.corr.null))
cd8.mass.null.sd = mean(apply(cd8.mass.corr.null,2,function(x) sd(x)))


cols<- c("blue","red")[(abs(mass.corr.spear.24[mass.corr.spear.24.rank])>2*cd8.mass.null.sd) + 1]
barplot(mass.corr.spear.24[mass.corr.spear.24.rank], col=cols,border=NA,ylim = c(-0.75,.75),xaxt='n', main = c('cd8 24 mass',as.character(length(which(mass.corr.spear.24>2*cd8.mass.null.sd))),as.character(length(which(mass.corr.spear.24< -2*cd8.mass.null.sd)))))
box()
abline(h = 2*cd8.mass.null.sd,lwd=2,col=1,lty=2)
abline(h = -2*cd8.mass.null.sd,lwd=2,col=1,lty=2)
abline(h=0,lwd=1)

# Plot Spearman v. Pearson for cd8 mass correlations 
plot(mass.corr.spear.24, mass.corr.pears, pch = 21, bg = rgb(0,0,1,0.3), xlab = 'Spearman Coefficient', ylab = 'Pearson Coefficeint', xlim = c(-0.8,0.8), ylim = c(-0.8,.8), cex.lab =1.2, cex.axis = 1.2,cex=0.75,col=rgb(0,0,0,0))
title(c('cd8 24 Mass', as.character(cor(mass.corr.spear.24,mass.corr.pears))))

# Normalized growth rate 24h
cd8@meta.data$norm = cd8@meta.data$mar/cd8@meta.data$mass

# norm -- Pearson
norm.corr.pears = apply(t(as.matrix(cd8@data[,cd8.24])),2,function(x) cor(x,cd8@meta.data$norm[cd8.24]))
norm.corr.pears[is.na(norm.corr.pears)]=0

# norm -- Spearman
norm.corr.spear.24 = apply(t(as.matrix(cd8@data[,cd8.24])),2,function(x) cor(x,cd8@meta.data$norm[cd8.24],method='spearman'))
norm.corr.spear.24[is.na(norm.corr.spear.24)]=0
norm.corr.spear.24.rank = sort(norm.corr.spear.24, decreasing=TRUE,index.return=TRUE)$ix
norm.ranked.list = norm.corr.spear.24[norm.corr.spear.24.rank]
names(norm.ranked.list) = toupper(names(norm.corr.spear.24)[norm.corr.spear.24.rank])
norm.ranked.enrichment = fgsea(go.set,norm.ranked.list,nperm=10000,maxSize = 500)
norm.ranked.enrichment <- norm.ranked.enrichment[norm.ranked.enrichment$padj< 0.1]
norm.ranked.enrichment$NES[is.na(norm.ranked.enrichment$NES)]=0
norm.ranked.enrichment = norm.ranked.enrichment[sort(norm.ranked.enrichment$NES,index.return=TRUE, decreasing = FALSE)$ix,1:6]

# Null Distribution for norm 
cd8.norm.null = sapply(1:10, function(x) sample(cd8.24,length(cd8.24),replace=FALSE))
cd8.norm.corr.null = sapply(1:10,function(x) (apply(cd8@data[,cd8.24],1,function(y) cor(y,cd8@meta.data$norm[cd8.norm.null[,x]],method='spearman'))))
cd8.norm.corr.null[is.na(cd8.norm.corr.null)]=0
cd8.norm.null.mean = mean(colMeans(cd8.norm.corr.null))
cd8.norm.null.sd = mean(apply(cd8.norm.corr.null,2,function(x) sd(x)))


cols<- c("blue","red")[(abs(norm.corr.spear.24[norm.corr.spear.24.rank])>2*cd8.norm.null.sd) + 1]
barplot(norm.corr.spear.24[norm.corr.spear.24.rank], col=cols,border=NA,ylim = c(-0.75,.75),xaxt='n', main = c('cd8 24 norm',as.character(length(which(norm.corr.spear.24>2*cd8.norm.null.sd))),as.character(length(which(norm.corr.spear.24< -2*cd8.norm.null.sd)))))
box()
abline(h = 2*cd8.norm.null.sd,lwd=2,col=1,lty=2)
abline(h = -2*cd8.norm.null.sd,lwd=2,col=1,lty=2)
abline(h=0,lwd=1)

# Plot Spearman v. Pearson for cd8 norm correlations 
plot(norm.corr.spear.24, norm.corr.pears, pch = 21, bg = rgb(0,0,1,0.3), xlab = 'Spearman Coefficient', ylab = 'Pearson Coefficeint', xlim = c(-0.8,0.8), ylim = c(-0.8,.8), cex.lab =1.2, cex.axis = 1.2,cex=0.75,col=rgb(0,0,0,0))
title(c('cd8 24 norm', as.character(cor(norm.corr.spear.24,norm.corr.pears))))

# 48h analysis
# Mass -- Pearson
mass.corr.pears = apply(t(as.matrix(cd8@data[,cd8.48])),2,function(x) cor(x,cd8@meta.data$mass[cd8.48]))
mass.corr.pears[is.na(mass.corr.pears)]=0

# Mass -- Spearman
mass.corr.spear.48 = apply(t(as.matrix(cd8@data[,cd8.48])),2,function(x) cor(x,cd8@meta.data$mass[cd8.48],method='spearman'))
mass.corr.spear.48[is.na(mass.corr.spear.48)]=0
mass.corr.spear.48.rank = sort(mass.corr.spear.48, decreasing=TRUE,index.return=TRUE)$ix
mass.ranked.list = mass.corr.spear.48[mass.corr.spear.48.rank]
names(mass.ranked.list) = toupper(names(mass.corr.spear.48)[mass.corr.spear.48.rank])
mass.ranked.enrichment = fgsea(go.set,mass.ranked.list,nperm=10000,maxSize = 500)
mass.ranked.enrichment <- mass.ranked.enrichment[mass.ranked.enrichment$padj< 0.1]
mass.ranked.enrichment$NES[is.na(mass.ranked.enrichment$NES)]=0
mass.ranked.enrichment = mass.ranked.enrichment[sort(mass.ranked.enrichment$NES,index.return=TRUE, decreasing = FALSE)$ix,1:6]


# Null Distribution for mass 
cd8.mass.null = sapply(1:10, function(x) sample(cd8.48,length(cd8.48),replace=FALSE))
cd8.mass.corr.null = sapply(1:10,function(x) (apply(cd8@data[,cd8.48],1,function(y) cor(y,cd8@meta.data$mass[cd8.mass.null[,x]],method='spearman'))))
cd8.mass.corr.null[is.na(cd8.mass.corr.null)]=0
cd8.mass.null.mean = mean(colMeans(cd8.mass.corr.null))
cd8.mass.null.sd = mean(apply(cd8.mass.corr.null,2,function(x) sd(x)))


cols<- c("blue","red")[(abs(mass.corr.spear.48[mass.corr.spear.48.rank])>2*cd8.mass.null.sd) + 1]
barplot(mass.corr.spear.48[mass.corr.spear.48.rank], col=cols,border=NA,ylim = c(-0.75,.75),xaxt='n', main = c('cd8 48 mass',as.character(length(which(mass.corr.spear.48>2*cd8.mass.null.sd))),as.character(length(which(mass.corr.spear.48< -2*cd8.mass.null.sd)))))
box()
abline(h = 2*cd8.mass.null.sd,lwd=2,col=1,lty=2)
abline(h = -2*cd8.mass.null.sd,lwd=2,col=1,lty=2)
abline(h=0,lwd=1)

# Plot Spearman v. Pearson for cd8 mass correlations 
plot(mass.corr.spear.48, mass.corr.pears, pch = 21, bg = rgb(0,0,1,0.3), xlab = 'Spearman Coefficient', ylab = 'Pearson Coefficeint', xlim = c(-0.8,0.8), ylim = c(-0.8,.8), cex.lab =1.2, cex.axis = 1.2,cex=0.75,col=rgb(0,0,0,0))
title(c('cd8 48 Mass', as.character(cor(mass.corr.spear.48,mass.corr.pears))))

# Normalized growth rate 48h
cd8@meta.data$norm = cd8@meta.data$mar/cd8@meta.data$mass

# norm -- Pearson
norm.corr.pears = apply(t(as.matrix(cd8@data[,cd8.48])),2,function(x) cor(x,cd8@meta.data$norm[cd8.48]))
norm.corr.pears[is.na(norm.corr.pears)]=0

# norm -- Spearman
norm.corr.spear.48 = apply(t(as.matrix(cd8@data[,cd8.48])),2,function(x) cor(x,cd8@meta.data$norm[cd8.48],method='spearman'))
norm.corr.spear.48[is.na(norm.corr.spear.48)]=0
norm.corr.spear.48.rank = sort(norm.corr.spear.48, decreasing=TRUE,index.return=TRUE)$ix
norm.ranked.list = norm.corr.spear.48[norm.corr.spear.48.rank]
names(norm.ranked.list) = toupper(names(norm.corr.spear.48)[norm.corr.spear.48.rank])
norm.ranked.enrichment = fgsea(go.set,norm.ranked.list,nperm=10000,maxSize = 500)
norm.ranked.enrichment <- norm.ranked.enrichment[norm.ranked.enrichment$padj< 0.1]
norm.ranked.enrichment$NES[is.na(norm.ranked.enrichment$NES)]=0
norm.ranked.enrichment = norm.ranked.enrichment[sort(norm.ranked.enrichment$NES,index.return=TRUE, decreasing = FALSE)$ix,1:6]

# Null Distribution for norm 
cd8.norm.null = sapply(1:10, function(x) sample(cd8.48,length(cd8.48),replace=FALSE))
cd8.norm.corr.null = sapply(1:10,function(x) (apply(cd8@data[,cd8.48],1,function(y) cor(y,cd8@meta.data$norm[cd8.norm.null[,x]],method='spearman'))))
cd8.norm.corr.null[is.na(cd8.norm.corr.null)]=0
cd8.norm.null.mean = mean(colMeans(cd8.norm.corr.null))
cd8.norm.null.sd = mean(apply(cd8.norm.corr.null,2,function(x) sd(x)))


cols<- c("blue","red")[(abs(norm.corr.spear.48[norm.corr.spear.48.rank])>2*cd8.norm.null.sd) + 1]
barplot(norm.corr.spear.48[norm.corr.spear.48.rank], col=cols,border=NA,ylim = c(-0.75,.75),xaxt='n', main = c('cd8 48 norm',as.character(length(which(norm.corr.spear.48>2*cd8.norm.null.sd))),as.character(length(which(norm.corr.spear.48< -2*cd8.norm.null.sd)))))
box()
abline(h = 2*cd8.norm.null.sd,lwd=2,col=1,lty=2)
abline(h = -2*cd8.norm.null.sd,lwd=2,col=1,lty=2)
abline(h=0,lwd=1)

# Plot Spearman v. Pearson for cd8 norm correlations 
plot(norm.corr.spear.48, norm.corr.pears, pch = 21, bg = rgb(0,0,1,0.3), xlab = 'Spearman Coefficient', ylab = 'Pearson Coefficeint', xlim = c(-0.8,0.8), ylim = c(-0.8,.8), cex.lab =1.2, cex.axis = 1.2,cex=0.75,col=rgb(0,0,0,0))
title(c('cd8 48 norm', as.character(cor(norm.corr.spear.48,norm.corr.pears))))

# Figure 3c -- Correlation with genes from Kimmerling et al., Nat. Comm. 2016
cc.genes = as.matrix(read.table("Cell_cycle_CD8.txt",sep="\t",header = TRUE))
cc.idx = match(cc.genes,rownames(cd8@data))
cc.idx = cc.idx[!is.na(cc.idx)]

# Create null distributions for correlation coefficients (limited to 10 iterations to reduce run time)
cd8.24.null.idx = sapply(1:10,function(x) sample(cd8.24,length(cd8.24),replace=FALSE))
null.mat.24 = sapply(1:10,function(x) (apply(cd8@data[cc.idx,cd8.24],1,function(y) cor(y,cd8@meta.data$mass[cd8.24.null.idx[,x]]))))
null.mean.24 = mean(colMeans(null.mat.24))
null.sd.24 = mean(apply(null.mat.24,2,function(x) sd(x)))
null.vec.24 = rnorm(10,mean=null.mean.24,sd=null.sd.24)

# 48h 
cd8.48.null.idx = sapply(1:10,function(x) sample(cd8.48,length(cd8.48),replace=FALSE))
null.mat.48 = sapply(1:10,function(x) (apply(cd8@data[cc.idx,cd8.48],1,function(y) cor(y,cd8@meta.data$mass[cd8.48.null.idx[,x]]))))
null.mean.48 = mean(colMeans(null.mat.48))
null.sd.48 = mean(apply(null.mat.48,2,function(x) sd(x)))
null.vec.48 = rnorm(10,mean=null.mean.48,sd=null.sd.48)

# Find correlations between mass and cluster genes for both time points
cl3.corr.24 = apply(cd8@data[cc.idx,cd8.24],1,function(x) cor(x,cd8@meta.data$mass[cd8.24],method='spearman'))
cl3.corr.48 = apply(cd8@data[cc.idx,cd8.48],1,function(x) cor(x,cd8@meta.data$mass[cd8.48],method = 'spearman'))
boxplot(null.vec.24,cl3.corr.24,null.vec.48,cl3.corr.48, ylim = c(-0.7,0.7), outline=FALSE,col=c(rgb(0,0,0,0.5),rgb(0,0,1,1),rgb(0,0,0,0.5),rgb(1,0,0,1)), names = c("24h Null","24h Mass","48h Null","48 Mass"), ylab = "Pearson Coefficient",lwd=2,cex.lab = 1.5, cex.axis=1.5)

