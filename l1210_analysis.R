library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library('fgsea')


# Figure 2a -- biophysical plot --------------------------------------------
par(pty='s')
# Includes additional L1210 cells that were not collected for scRNA-seq but had mass/MAR measured
l1210.growth1 = read.table("coefs_l1210_serial1.txt", sep = "\t",header = FALSE)
l1210.growth2 = read.table("coefs_l1210_serial2.txt",sep = "\t",header = FALSE)
l1210.growth = rbind(l1210.growth1,l1210.growth2)
plot.idx1 = which(l1210.growth[,5]<0.65 & l1210.growth[,8]>10&l1210.growth[,2]<90)
plot(l1210.growth[plot.idx1,2],l1210.growth[plot.idx1,3], ylim = c(0,8), xlim = c(30,90), pch =21, bg = rgb(0,0,1,0.65), xlab = "Mass (pg)", ylab = "MAR (pg/h)", cex = 1.5,cex.axis=2.2,cex.lab=2.2)



# Set up Seurat object ----------------------------------------------------
## READ IN RAW DATA AND META DATA
raw.data = read.table("l1210_trap_serial_rsem.txt",sep = "\t", header = TRUE, row.names=1)
raw.data = log(raw.data+1)
meta = read.table("qc_trap_serial.txt",sep = "\t", header = TRUE, row.names = 1)


## SET UP SEURAT OBJECT
l1210s = CreateSeuratObject(raw.data=raw.data,project = "l1210s",min.cells = 5,names.field=1,min.genes=4000,is.expr=0,meta.data=meta, do.scale=FALSE)
l1210s = FindVariableGenes(object = l1210s,mean.function = ExpMean,dispersion.function=LogVMR, x.low.cutoff=0.5,y.cutoff=0.5,do.plot=FALSE)
l1210s = ScaleData(l1210s)
# PARE DOWN TO CELLS WITH HIGH QUALITY GROWTH MEASUREMENTS
l1210s.phys = which(l1210s@meta.data$mass>-100 & l1210s@meta.data$mass<85 & l1210s@meta.data$system == 2)




# Supplementary tables 1 & 2 ----------------------------------------------
# Enrichment reference set 
go.set = gmtPathways('c5.all.v6.1.symbols.gmt')

# Mass -- spearman
mass.corr.spear = apply(t(as.matrix(l1210s@data[,l1210s.phys])),2,function(x) cor(x,l1210s@meta.data$mass[l1210s.phys],method='spearman'))
mass.corr.spear[is.na(mass.corr.spear)]=0
mass.corr.spear.rank = sort(mass.corr.spear, decreasing=TRUE,index.return=TRUE)$ix
mass.ranked.list = mass.corr.spear[mass.corr.spear.rank]
names(mass.ranked.list) = toupper(names(mass.corr.spear)[mass.corr.spear.rank])
mass.ranked.enrichment = fgsea(go.set,mass.ranked.list,nperm=15000,maxSize = 500)
mass.ranked.enrichment$NES[is.na(mass.ranked.enrichment$NES)]=0
mass.ranked.enrichment <- mass.ranked.enrichment[which(mass.ranked.enrichment$padj<0.1)]

# Normalized growth rate
l1210s@meta.data$norm = l1210s@meta.data$mar/l1210s@meta.data$mass
# norm -- spearman
norm.corr.spear = apply(t(as.matrix(l1210s@data[,l1210s.phys])),2,function(x) cor(x,l1210s@meta.data$norm[l1210s.phys],method='spearman'))
norm.corr.spear[is.na(norm.corr.spear)]=0
norm.corr.spear.rank = sort(norm.corr.spear, decreasing=TRUE,index.return=TRUE)$ix
norm.ranked.list = norm.corr.spear[norm.corr.spear.rank]
names(norm.ranked.list) = toupper(names(norm.corr.spear)[norm.corr.spear.rank])
norm.ranked.enrichment = fgsea(go.set,norm.ranked.list,nperm=10000,maxSize = 500)
norm.ranked.enrichment$NES[is.na(norm.ranked.enrichment$NES)]=0
norm.ranked.enrichment <- norm.ranked.enrichment[which(norm.ranked.enrichment$padj<0.1)]




# Supplementary Figure 4 & 5 ----------------------------------------------

# Null Distribution for mass 
l1210.mass.null = sapply(1:10, function(x) sample(l1210s.phys,length(l1210s.phys),replace=FALSE))
l1210.mass.corr.null = sapply(1:10,function(x) (apply(l1210s@data[,l1210s.phys],1,function(y) cor(y,l1210s@meta.data$mass[l1210.mass.null[,x]],method='spearman'))))
l1210.mass.null.mean = mean(colMeans(l1210.mass.corr.null))
l1210.mass.null.sd = mean(apply(l1210.mass.corr.null,2,function(x) sd(x)))

# Setting significance threshold for mass correlations
cols<- c("blue","red")[(abs(mass.corr.spear[mass.corr.spear.rank])>2*l1210.mass.null.sd) + 1]
barplot(mass.corr.spear[mass.corr.spear.rank], col=cols,border=NA,ylim = c(-0.75,.75),xaxt='n', main = c('L1210 mass',as.character(length(which(mass.corr.spear>2*l1210.mass.null.sd))),as.character(length(which(mass.corr.spear< -2*l1210.mass.null.sd)))))
box()
abline(h = 2*l1210.mass.null.sd,lwd=2,col=1,lty=2)
abline(h = -2*l1210.mass.null.sd,lwd=2,col=1,lty=2)
abline(h=0,lwd=1)

# Plot spearman v. Pearson for L1210 mass correlations 
mass.corr.pears = apply(t(as.matrix(l1210s@data[,l1210s.phys])),2,function(x) cor(x,l1210s@meta.data$mass[l1210s.phys]))
plot(mass.corr.spear, mass.corr.pears, pch = 21, bg = rgb(0,0,1,0.3), xlab = 'spearman Coefficient', ylab = 'Pearson Coefficeint', xlim = c(-0.8,0.8), ylim = c(-0.8,.8), cex.lab =1.2, cex.axis = 1.2,cex=0.75,col=rgb(0,0,0,0))
title(c('L1210 Mass', as.character(cor(mass.corr.spear,mass.corr.pears))))

# Null Distribution for norm 
l1210.norm.null = sapply(1:10, function(x) sample(l1210s.phys,length(l1210s.phys),replace=FALSE))
l1210.norm.corr.null = sapply(1:10,function(x) (apply(l1210s@data[,l1210s.phys],1,function(y) cor(y,l1210s@meta.data$norm[l1210.norm.null[,x]],method='spearman'))))
l1210.norm.null.mean = mean(colMeans(l1210.norm.corr.null))
l1210.norm.null.sd = mean(apply(l1210.norm.corr.null,2,function(x) sd(x)))

# Setting significance threshold for normalized MAR correlation
cols<- c("blue","red")[(abs(norm.corr.spear[norm.corr.spear.rank])>2*l1210.norm.null.sd) + 1]
barplot(norm.corr.spear[norm.corr.spear.rank], col=cols,border=NA,ylim = c(-0.75,.75),xaxt='n', main = c('L1210 norm',as.character(length(which(norm.corr.spear>2*l1210.norm.null.sd))),as.character(length(which(norm.corr.spear< -2*l1210.norm.null.sd)))))
box()
abline(h = 2*l1210.norm.null.sd,lwd=2,col=1,lty=2)
abline(h = -2*l1210.norm.null.sd,lwd=2,col=1,lty=2)
abline(h=0,lwd=1)

# Plot spearman v. Pearson for L1210 norm correlations 
norm.corr.pears = apply(t(as.matrix(l1210s@data[,l1210s.phys])),2,function(x) cor(x,l1210s@meta.data$norm[l1210s.phys]))
plot(norm.corr.spear, norm.corr.pears, pch = 21, bg = rgb(0,0,1,0.3), xlab = 'spearman Coefficient', ylab = 'Pearson Coefficeint', xlim = c(-0.8,0.8), ylim = c(-0.8,.8), cex.lab =1.2, cex.axis = 1.2,cex=0.75,col=rgb(0,0,0,0))
title(c('L1210 norm', as.character(cor(norm.corr.spear,norm.corr.pears))))





# Figure 2b -- Cell cycle heat map -----------------------------------------
proper<-function(x) paste(toupper(substr(x, 1, 1)), tolower(substring(x, 2)), sep ='')
# Find genes in chromosome segregation and DNA replication terms with signficant correlation with mass
chr.list <- as.matrix(read.table('chromosome_segregation.txt', header = TRUE, sep ='\t'))
# Find those with significant positive correlation with gene expression based on null SD
chr.list.keep <- which(mass.ranked.list[toupper(chr.list)]>2*l1210.mass.null.sd)
chr.list <- proper(chr.list[chr.list.keep])

# Now the negatively correlated 
dna.list <- as.matrix(read.table('dna_replication.txt', header = TRUE, sep ='\t'))
dna.keep <- which(mass.ranked.list[toupper(dna.list)]< -2*l1210.mass.null.sd)
dna.list <- proper(dna.list[dna.keep])

all.genes <- c(chr.list, dna.list)
# Order cells by mass
mass.sort = sort(l1210s@meta.data$mass[l1210s.phys],decreasing=FALSE,index.return=TRUE)$ix
c1 <- colnames(l1210s@data)[l1210s.phys][mass.sort]
DoHeatmap(l1210s, genes.use = all.genes, cells.use = c1, disp.min = -1.5, disp.max = 1.5, cex.col = 0  )
# Bar plot of masses for heat map
barplot(l1210s@meta.data$mass[l1210s.phys][mass.sort],space = 1,col=1, ylim = c(20,90))

