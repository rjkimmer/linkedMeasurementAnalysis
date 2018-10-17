# Analysis of various cell data sets based on Mackosko 2015 cell cycle analysis (Supp. Fig. 5)

# Load in the cell-cycle gene list 
g1s <- as.matrix(read.table('G1S.txt', header = TRUE))
g2m <- as.matrix(read.table('g2m.txt', header = TRUE))


# Load sequencing data 
library(Seurat)
proper<-function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))

## READ IN RAW DATA AND META DATA
raw.data = read.table("CD8_rsem.txt",sep = "\t", header = TRUE, row.names=1)
meta = read.table("cd8_qc2.txt",sep = "\t", header = TRUE, row.names = 1)

## SET UP SEURAT OBJECT
cd8 = CreateSeuratObject(raw.data=raw.data,project = "cd8",min.cells = 19,names.field=1,min.genes=1000,is.expr=0,meta.data=meta, do.scale=FALSE)
cd8 = NormalizeData(object= cd8, normalization.method="LogNormalize")
cd8 = FindVariableGenes(object = cd8,mean.function = ExpMean,dispersion.function=LogVMR, x.low.cutoff=0.5,y.cutoff=0,do.plot=FALSE)
cd8 = ScaleData(object=cd8)
cd8.24 = which(cd8@meta.data$time==24 & cd8@meta.data$mar>-100)
cd8.48 = which(cd8@meta.data$time==48 & cd8@meta.data$mar>-100)
cd8.tot = c(cd8.24,cd8.48)


# Determine scoring for each cell cycle phase
genes <- toupper(row.names(cd8@data))

#  G1S
g1s.present <- intersect(g1s,genes)
g1s.data <- cd8@data[proper(g1s.present),cd8.24]
g1s.rough.scores <- colMeans(g1s.data)
g1s.cors <- sapply(1:length(g1s.present), function(x) cor(g1s.rough.scores, (g1s.data[x,])))
g1s.keep <- g1s.present[which(g1s.cors>0.3)]
g1s.score <- colMeans(cd8@data[proper(g1s.keep),cd8.24])
g1s.score.normal <- scale(g1s.score, center = TRUE, scale = TRUE)

# G2M
g2m.present <- intersect(g2m,genes)
g2m.data <- cd8@data[proper(g2m.present), cd8.24]
g2m.rough.scores <- colMeans(g2m.data)
g2m.cors <- sapply(1:length(g2m.present), function(x) cor(g2m.rough.scores, (g2m.data[x,])))
g2m.keep <- g2m.present[which(g2m.cors>0.3)]
g2m.score <- colMeans(cd8@data[proper(g2m.keep),cd8.24])
g2m.score.normal <- scale(g2m.score, center = TRUE, scale = TRUE)

# Append scores and find whether a cell is more well-defined as G1/S or G2/M
app.scores <- cbind(g1s.score.normal, g2m.score.normal)
sort.idx <- sort(cd8@meta.data$mass[cd8.24], index.return = TRUE)$ix
g1s.assign <- sapply(1:dim(app.scores)[1], function(x) app.scores[x,1]>app.scores[x,2])
app.scores2 <- MinMax(app.scores,-2,2)
my_palette = colorRampPalette(c("#FF00FF","#000000","#FFFF00"))(n=length(cd8.24))
heatmap.2(t(app.scores2[sort.idx,]), Rowv = NA, Colv = NA, dendrogram = 'none', trace = 'none', col = my_palette, key = FALSE)

# Look at masses based on assignment
mass.sub <- cd8@meta.data$mass[cd8.24]
par(pty = 's')
boxplot(mass.sub[g1s.assign], mass.sub[!g1s.assign], outline= FALSE, names = c('G1S', 'G2M'), ylab = 'Buoyant Mass (pg)', ylim= c(20,90), col = c('blue', 'red'), cex.axis = 2.25, cex.lab = 2.5, lwd = 2.5)


#### And again for the 48h time point 

# Determine scoring for each cell cycle phase
genes <- toupper(row.names(cd8@data))

#  G1S
g1s.present <- intersect(g1s,genes)
g1s.data <- cd8@data[proper(g1s.present),cd8.48]
g1s.rough.scores <- colMeans(g1s.data)
g1s.cors <- sapply(1:length(g1s.present), function(x) cor(g1s.rough.scores, (g1s.data[x,])))
g1s.keep <- g1s.present[which(g1s.cors>0.3)]
g1s.score <- colMeans(cd8@data[proper(g1s.keep),cd8.48])
g1s.score.normal <- scale(g1s.score, center = TRUE, scale = TRUE)


# G2M
g2m.present <- intersect(g2m,genes)
g2m.data <- cd8@data[proper(g2m.present), cd8.48]
g2m.rough.scores <- colMeans(g2m.data)
g2m.cors <- sapply(1:length(g2m.present), function(x) cor(g2m.rough.scores, (g2m.data[x,])))
g2m.keep <- g2m.present[which(g2m.cors>0.3)]
g2m.score <- colMeans(cd8@data[proper(g2m.keep),cd8.48])
g2m.score.normal <- scale(g2m.score, center = TRUE, scale = TRUE)


# Append scores and find whether a cell is more well-defined as G1/S or G2/M
app.scores <- cbind(g1s.score.normal,g2m.score.normal)
sort.idx <- sort(cd8@meta.data$mass[cd8.48], index.return = TRUE)$ix
g1s.assign <- sapply(1:dim(app.scores)[1], function(x) app.scores[x,1]>app.scores[x,2])
app.scores2 <- MinMax(app.scores,-2,2)
my_palette = colorRampPalette(c("#FF00FF","#000000","#FFFF00"))(n=length(cd8.48))
heatmap.2(t(app.scores2[sort.idx,]), Rowv = NA, Colv = NA, dendrogram = 'none', trace = 'none', col = my_palette, key = FALSE)

# Look at masses based on assignment
mass.sub <- cd8@meta.data$mass[cd8.48]
boxplot(mass.sub[g1s.assign], mass.sub[!g1s.assign], outline= FALSE, names = c('G1S', 'G2M'), ylab = 'Buoyant Mass (pg)', ylim= c(20,90), col = c('blue', 'red'), cex.axis = 2.25, cex.lab = 2.5, lwd = 2.5)


