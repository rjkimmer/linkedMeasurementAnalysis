# Analysis of various cell data sets based on Mackosko 2015 cell cycle analysis (Supp. Fig. 5)

# Load in the cell-cycle gene list 
g1s <- as.matrix(read.table('G1S.txt', header = TRUE))
g2m <- as.matrix(read.table('g2m.txt', header = TRUE))

# Load sequencing data -- L1210 cells 
library(Seurat)
proper<-function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))

## READ IN RAW DATA AND META DATA
raw.data = read.table("l1210_trap_serial_rsem.txt",sep = "\t", header = TRUE, row.names=1)
raw.data = raw.data[,146:241]
#raw.data = log(raw.data+1)
meta = read.table("qc_trap_serial.txt",sep = "\t", header = TRUE, row.names = 1)
meta = meta[147:242,]


## SET UP SEURAT OBJECT

l1210s = CreateSeuratObject(raw.data=raw.data,project = "l1210s",min.cells = 5,names.field=1,min.genes=4000,is.expr=0,meta.data=meta, do.scale=FALSE)
l1210s = NormalizeData(object= l1210s, normalization.method="LogNormalize")
l1210s = FindVariableGenes(object = l1210s,mean.function = ExpMean,dispersion.function=LogVMR, x.low.cutoff=0.5,y.cutoff=0.5,do.plot=FALSE)
l1210s = ScaleData(l1210s)


# PARE DOWN TO CELLS WITH HIGH QUALITY GROWTH MEASUREMENTS
l1210s.phys = which(l1210s@meta.data$mass>-100 & l1210s@meta.data$mass<85 & l1210s@meta.data$system == 2)

# Determine scoring for each cell cycle phase
genes <- toupper(row.names(l1210s@data))

#  G1S
g1s.present <- intersect(g1s,genes)
g1s.data <- l1210s@data[proper(g1s.present),l1210s.phys]
g1s.rough.scores <- colMeans(g1s.data)
g1s.cors <- sapply(1:length(g1s.present), function(x) cor(g1s.rough.scores, (g1s.data[x,])))
g1s.keep <- g1s.present[which(g1s.cors>0.3)]
g1s.score <- colMeans(l1210s@data[proper(g1s.keep),l1210s.phys])
g1s.score.normal <- scale(g1s.score, center = TRUE, scale = TRUE)

# G2M
g2m.present <- intersect(g2m,genes)
g2m.data <- l1210s@data[proper(g2m.present), l1210s.phys]
g2m.rough.scores <- colMeans(g2m.data)
g2m.cors <- sapply(1:length(g2m.present), function(x) cor(g2m.rough.scores, (g2m.data[x,])))
g2m.keep <- g2m.present[which(g2m.cors>0.3)]
g2m.score <- colMeans(l1210s@data[proper(g2m.keep),l1210s.phys])
g2m.score.normal <- scale(g2m.score, center = TRUE, scale = TRUE)

s1 <- c(0,1)
s2 <- c(1,0)

# Append scores and find whether a cell is more well-defined as G1/S or G2/M
app.scores <- cbind(g1s.score.normal, g2m.score.normal)
sort.idx <- sort(l1210s@meta.data$mass[l1210s.phys], index.return = TRUE)$ix
g1s.assign <- sapply(1:dim(app.scores)[1], function(x) app.scores[x,1]>app.scores[x,2])

app.scores2 <- MinMax(app.scores, -2,2)
my_palette = colorRampPalette(c("#FF00FF","#000000","#FFFF00"))(n=length(l1210s.phys))
heatmap.2(t(app.scores2[sort.idx,]), Rowv = NA, Colv = NA, dendrogram = 'none', trace = 'none', col = my_palette, key = FALSE)

# Look at masses based on assignmen
mass.sub <- l1210s@meta.data$mass[l1210s.phys]
par(pty = 's')
boxplot(mass.sub[g1s.assign], mass.sub[!g1s.assign], outline= FALSE, names = c('G1S', 'G2M'), ylab = 'Buoyant Mass (pg)', ylim= c(30,90), col = c('blue', 'red'), cex.axis = 2.25, cex.lab = 2.5, lwd = 2.5)


