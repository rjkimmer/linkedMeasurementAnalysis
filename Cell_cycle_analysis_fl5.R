# Analysis of various cell data sets based on Mackosko 2015 cell cycle analysis (Supp. Fig. 5)

# Load in the cell-cycle gene list 
g1s <- as.matrix(read.table('G1S.txt', header = TRUE))
g2m <- as.matrix(read.table('g2m.txt', header = TRUE))


# Load sequencing data -- fl5 cells 
library(Seurat)
proper<-function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
## READ IN RAW DATA AND META DATA
raw.data = read.table("fl5_serial_rsem3.txt",sep = "\t", header = TRUE, row.names=1)
raw.data = log(raw.data+1)
meta = read.table("qc_fl5_serial3.txt",sep = "\t", header = TRUE, row.names = 1)


## SET UP SEURAT OBJECT
fl5s = CreateSeuratObject(raw.data=raw.data,project = "fl5s",min.cells = 13,names.field=1,min.genes=4000,is.expr=0,meta.data=meta, do.scale=FALSE)
fl5s = FindVariableGenes(object = fl5s,mean.function = ExpMean,dispersion.function=LogVMR, x.low.cutoff=0.5,y.cutoff=0.5,do.plot=FALSE)
fl5s = ScaleData(object=fl5s)

# PARE DOWN TO CELLS WITH HIGH QUALITY GROWTH MEASUREMENTS
fl5s.phys1 = which(fl5s@meta.data$mass>-100 & fl5s@meta.data$mass<80 & fl5s@meta.data$day ==1)

# Determine scoring for each cell cycle phase
genes <- toupper(row.names(fl5s@data))

#  G1S
g1s.present <- intersect(g1s,genes)
g1s.data <- fl5s@data[proper(g1s.present),fl5s.phys1]
g1s.rough.scores <- colMeans(g1s.data)
g1s.cors <- sapply(1:length(g1s.present), function(x) cor(g1s.rough.scores, t(g1s.data[x,])))
g1s.keep <- g1s.present[which(g1s.cors>0.3)]
g1s.score <- colMeans(fl5s@data[proper(g1s.keep),fl5s.phys1])
g1s.score.normal <- scale(g1s.score, center = TRUE, scale = TRUE)

# G2M
g2m.present <- intersect(g2m,genes)
g2m.data <- fl5s@data[proper(g2m.present), fl5s.phys1]
g2m.rough.scores <- colMeans(g2m.data)
g2m.cors <- sapply(1:length(g2m.present), function(x) cor(g2m.rough.scores, t(g2m.data[x,])))
g2m.keep <- g2m.present[which(g2m.cors>0.3)]
g2m.score <- colMeans(fl5s@data[proper(g2m.keep),fl5s.phys1])
g2m.score.normal <- scale(g2m.score, center = TRUE, scale = TRUE)


# Append scores and find whether a cell is more well-defined as G1/S or G2/M
app.scores <- cbind(g1s.score.normal, g2m.score.normal)
sort.idx <- sort(fl5s@meta.data$mass[fl5s.phys1], index.return = TRUE)$ix
g1s.assign <- sapply(1:dim(app.scores)[1], function(x) app.scores[x,1]>app.scores[x,2])
app.scores2 <- MinMax(app.scores,-2,2)

my_palette = colorRampPalette(c("#FF00FF","#000000","#FFFF00"))(n=length(fl5s.phys1))
heatmap.2(t(app.scores2[sort.idx,]), Rowv = NA, Colv = NA, dendrogram = 'none', trace = 'none', col = my_palette, key = FALSE)

# Look at masses based on assignmen
mass.sub <- fl5s@meta.data$mass[fl5s.phys1]
par(pty = 's')
boxplot(mass.sub[g1s.assign], mass.sub[!g1s.assign], outline= FALSE, names = c('G1S', 'G2M'), ylab = 'Buoyant Mass (pg)', ylim= c(30,90), col = c('blue', 'red'), cex.axis = 2.25, cex.lab = 2.5, lwd = 2.5)




