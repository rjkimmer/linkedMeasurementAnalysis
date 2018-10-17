library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library('fgsea')
proper<-function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))


# Figure 2 -- Single-cell mass/MAR plot -----------------------------------

fl5.growth1 = read.table("coefs_fl5_serial.txt", sep = "\t",header = FALSE)
plot.idx1 = which(fl5.growth1[,5]<0.75 & fl5.growth1[,8]>10 & fl5.growth1[,2]<85)
plot(fl5.growth1[plot.idx1,2],fl5.growth1[plot.idx1,3], ylim = c(0,8), xlim = c(30,85), pch =21, bg = rgb(0,0,1,0.65), xlab = "Mass (pg)", ylab = "MAR (pg/h)", cex = 1.5,cex.axis=2.2,cex.lab=2.2)


# Setup Seurat object -----------------------------------------------------


## READ IN RAW DATA AND META DATA
raw.data = read.table("fl5_serial_rsem3.txt",sep = "\t", header = TRUE, row.names=1)
raw.data = log(raw.data+1)
meta = read.table("qc_fl5_serial3.txt",sep = "\t", header = TRUE, row.names = 1)

## SET UP SEURAT OBJECT
fl5s = CreateSeuratObject(raw.data=raw.data,project = "fl5s",min.cells = 13,names.field=1,min.genes=4000,is.expr=0,meta.data=meta, do.scale=FALSE)
#fl5s = NormalizeData(object= fl5s, normalization.method="LogNormalize")
fl5s = FindVariableGenes(object = fl5s,mean.function = ExpMean,dispersion.function=LogVMR, x.low.cutoff=0.5,y.cutoff=0.5,do.plot=FALSE)
fl5s = ScaleData(object=fl5s)

# PARE DOWN TO CELLS WITH HIGH QUALITY GROWTH MEASUREMENTS
fl5s.phys1 = which(fl5s@meta.data$mass>-100 & fl5s@meta.data$mass<80 & fl5s@meta.data$day ==1)
# Update cell identities 
fl5s@ident <- factor(levels = c(levels(fl5s@ident), 'FL51'))
fl5s@ident[which(fl5s@ident =='FL51a')] <- 'FL51'
fl5s@ident[which(fl5s@ident =='FL51b')] <- 'FL51'


# Supplementary tables 1 & 2 ----------------------------------------------
go.set = gmtPathways('c5.all.v6.1.symbols.gmt')

# Mass -- Spearman
mass.corr.spear = apply(t(as.matrix(fl5s@data[,fl5s.phys1])),2,function(x) cor(x,fl5s@meta.data$mass[fl5s.phys1],method='spearman'))
mass.corr.spear[is.na(mass.corr.spear)]=0
mass.corr.spear.rank = sort(mass.corr.spear, decreasing=TRUE,index.return=TRUE)$ix
mass.ranked.list = mass.corr.spear[mass.corr.spear.rank]
names(mass.ranked.list) = toupper(names(mass.corr.spear)[mass.corr.spear.rank])
mass.ranked.enrichment = fgsea(go.set,mass.ranked.list,nperm=15000,maxSize = 500)
mass.ranked.enrichment$NES[is.na(mass.ranked.enrichment$NES)]=0
mass.ranked.enrichment <- mass.ranked.enrichment[which(mass.ranked.enrichment$padj<0.1)]

# Normalized growth rate
fl5s@meta.data$norm = fl5s@meta.data$mar/fl5s@meta.data$mass

# norm -- Spearman
norm.corr.spear = apply(t(as.matrix(fl5s@data[,fl5s.phys1])),2,function(x) cor(x,fl5s@meta.data$norm[fl5s.phys1],method='spearman'))
norm.corr.spear[is.na(norm.corr.spear)]=0
norm.corr.spear.rank = sort(norm.corr.spear, decreasing=TRUE,index.return=TRUE)$ix
norm.ranked.list = norm.corr.spear[norm.corr.spear.rank]
names(norm.ranked.list) = toupper(names(norm.corr.spear)[norm.corr.spear.rank])
norm.ranked.enrichment = fgsea(go.set,norm.ranked.list,nperm=15000,maxSize = 500)
norm.ranked.enrichment$NES[is.na(norm.ranked.enrichment$NES)]=0
norm.ranked.enrichment <- norm.ranked.enrichment[norm.ranked.enrichment$padj<0.1]


# Supplementary figure 4 & 5 ----------------------------------------------
# Null Distribution for mass 
fl5.mass.null = sapply(1:10, function(x) sample(fl5s.phys1,length(fl5s.phys1),replace=FALSE))
fl5.mass.corr.null = sapply(1:10,function(x) (apply(fl5s@data[,fl5s.phys1],1,function(y) cor(y,fl5s@meta.data$mass[fl5.mass.null[,x]],method='spearman'))))
fl5.mass.corr.null[is.na(fl5.mass.corr.null)]=0
fl5.mass.null.mean = mean(colMeans(fl5.mass.corr.null))
fl5.mass.null.sd = mean(apply(fl5.mass.corr.null,2,function(x) sd(x)))

cols<- c("blue","red")[(abs(mass.corr.spear[mass.corr.spear.rank])>2*fl5.mass.null.sd) + 1]
barplot(mass.corr.spear[mass.corr.spear.rank], col=cols,border=NA,ylim = c(-0.75,.75),xaxt='n', main = c('fl5 mass',as.character(length(which(mass.corr.spear>2*fl5.mass.null.sd))),as.character(length(which(mass.corr.spear< -2*fl5.mass.null.sd)))))
box()
abline(h = 2*fl5.mass.null.sd,lwd=2,col=1,lty=2)
abline(h = -2*fl5.mass.null.sd,lwd=2,col=1,lty=2)
abline(h=0,lwd=1)

# Plot Spearman v. Pearson for fl5 mass correlations 
mass.corr.pears = apply(t(as.matrix(fl5s@data[,fl5s.phys1])),2,function(x) cor(x,fl5s@meta.data$mass[fl5s.phys1]))
plot(mass.corr.spear, mass.corr.pears, pch = 21, bg = rgb(0,0,1,0.3), xlab = 'Spearman Coefficient', ylab = 'Pearson Coefficeint', xlim = c(-0.8,0.8), ylim = c(-0.8,.8), cex.lab =1.2, cex.axis = 1.2,cex=0.75,col=rgb(0,0,0,0))
title(c('fl5 Mass', as.character(cor(mass.corr.spear,mass.corr.pears))))

# Null Distribution for norm 
fl5.norm.null = sapply(1:10, function(x) sample(fl5s.phys1,length(fl5s.phys1),replace=FALSE))
fl5.norm.corr.null = sapply(1:10,function(x) (apply(fl5s@data[,fl5s.phys1],1,function(y) cor(y,fl5s@meta.data$norm[fl5.norm.null[,x]],method='spearman'))))
fl5.norm.corr.null[is.na(fl5.norm.corr.null)]=0
fl5.norm.null.mean = mean(colMeans(fl5.norm.corr.null))
fl5.norm.null.sd = mean(apply(fl5.norm.corr.null,2,function(x) sd(x)))


cols<- c("blue","red")[(abs(norm.corr.spear[norm.corr.spear.rank])>2*fl5.norm.null.sd) + 1]
barplot(norm.corr.spear[norm.corr.spear.rank], col=cols,border=NA,ylim = c(-0.75,.75),xaxt='n', main = c('fl5 norm',as.character(length(which(norm.corr.spear>2*fl5.norm.null.sd))),as.character(length(which(norm.corr.spear< -2*fl5.norm.null.sd)))))
box()
abline(h = 2*fl5.norm.null.sd,lwd=2,col=1,lty=2)
abline(h = -2*fl5.norm.null.sd,lwd=2,col=1,lty=2)
abline(h=0,lwd=1)

# Plot Spearman v. Pearson for fl5 norm correlations 
norm.corr.pears = apply(t(as.matrix(fl5s@data[,fl5s.phys1])),2,function(x) cor(x,fl5s@meta.data$norm[fl5s.phys1]))
plot(norm.corr.spear, norm.corr.pears, pch = 21, bg = rgb(0,0,1,0.3), xlab = 'Spearman Coefficient', ylab = 'Pearson Coefficeint', xlim = c(-0.8,0.8), ylim = c(-0.8,.8), cex.lab =1.2, cex.axis = 1.2,cex=0.75,col=rgb(0,0,0,0))
title(c('fl5 norm', as.character(cor(norm.corr.spear,norm.corr.pears))))



# Figure 2b -- Cell cycle heat map ----------------------------------------

# Mass ranked heat map

mass.corr.spear = apply(t(as.matrix(fl5s@data[,fl5s.phys1])),2,function(x) cor(x,fl5s@meta.data$mass[fl5s.phys1],method='spearman'))
mass.corr.spear[is.na(mass.corr.spear)]=0
mass.corr.spear.rank = sort(mass.corr.spear, decreasing=TRUE,index.return=TRUE)$ix
mass.ranked.list = mass.corr.spear[mass.corr.spear.rank]

# Find genes in chromosome segregation and DNA replication GO terms with significant correlation with mass
chr.list <- as.matrix(read.table('chromosome_segregation.txt', header = TRUE, sep ='\t'))
# Find those with significant positive correlation with gene expression based on null SD
chr.list.keep <- which(mass.ranked.list[(chr.list)]>2*0.1178759)
chr.list <- proper(chr.list[chr.list.keep])

# Now the negatively correlated 
#dna.list <- mass.ranked.enrichment[mass.ranked.enrichment$pathway == "GO_DNA_REPLICATION"]$leadingEdge[[1]]
dna.list <- as.matrix(read.table('dna_replication.txt', header = TRUE, sep ='\t'))
dna.keep <- which(mass.ranked.list[(dna.list)]< -2*0.1178759)
dna.list <- proper(dna.list[dna.keep])

all.genes <- c(chr.list, dna.list)
# Order cells by mass
mass.sort = sort(fl5s@meta.data$mass[fl5s.phys1],decreasing=FALSE,index.return=TRUE)$ix
c1 <- colnames(fl5s@data)[fl5s.phys1][mass.sort]
DoHeatmap(fl5s, genes.use = all.genes, cells.use = c1, disp.min = -1.5, disp.max = 1.5, cex.col = 0)


# Bar plot of masses for heat map
barplot(fl5s@meta.data$mass[fl5s.phys1][mass.sort],space = 1,col=1, ylim = c(20,90))

# Supp. Fig. 7
# G1S scoring 

norm.corr.spear = apply(t(as.matrix(fl5s@data[,fl5s.phys1])),2,function(x) cor(x,fl5s@meta.data$norm[fl5s.phys1],method='spearman'))
norm.corr.spear[is.na(norm.corr.spear)]=0
norm.corr.spear.rank = sort(norm.corr.spear, decreasing=TRUE,index.return=TRUE)$ix
norm.ranked.list = norm.corr.spear[norm.corr.spear.rank]


# Re-run null distribution 
fl5.norm.null = sapply(1:10, function(x) sample(fl5s.phys1,length(fl5s.phys1),replace=FALSE))
fl5.norm.corr.null = sapply(1:10,function(x) (apply(fl5s@data[,fl5s.phys1],1,function(y) cor(y,fl5s@meta.data$norm[fl5.norm.null[,x]],method='spearman'))))
fl5.norm.corr.null[is.na(fl5.norm.corr.null)]=0
fl5.norm.null.mean = mean(colMeans(fl5.norm.corr.null))
fl5.norm.null.sd = mean(apply(fl5.norm.corr.null,2,function(x) sd(x)))

g1s.list <- as.matrix(read.table('g1_s_transition.txt', header = TRUE, sep = '\t'))
g1s.keep <- which(norm.ranked.list[(g1s.list)] > 2*fl5.norm.null.sd)
g1s.list <- g1s.list[g1s.keep]

g1s.score <- colMeans(as.matrix(fl5s@scale.data)[g1s.list, fl5s.phys1])
g1s.score <- MinMax(g1s.score, -1.0,1.0)
colfunc = colorRampPalette(c("blue","white","red"))
col.score = colfunc(length(g1s.score))[as.numeric(cut((g1s.score),breaks=length(g1s.score)))]
plot(fl5s@meta.data$mass[fl5s.phys1],fl5s@meta.data$norm[fl5s.phys1],pch=21,cex=1.5, xlab = "Mass (pg)", ylab = "MAR (pg/h)", bg=col.score, xlim = c(30,80), ylim = c(0,.12))


