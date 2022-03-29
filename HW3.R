#load libraries
library(Rtsne)
library(ggplot2)
library(scattermore)
library(gridExtra) 

#load MERFISH data
data <- read.csv('~/Dropbox/JHU/Courses/genomic-data-visualization/data/MERFISH_Slice2Replicate2_halfcortex.csv.gz')

#downsample to 5000
vi <- sample(data[,1],5000)
ds <- data[data$X %in% vi,]

#make matrix of position
pos <- ds[, c('x','y')]
rownames(pos) <- ds[,1]
plot(pos, pch='.') #change geometric primitive to point

#make matrix of gene expression data without position
gexp <- ds[, 4:ncol(ds)]
rownames(gexp) <- ds[,1]

# for each gene, count how many cells express that gene
cellexp <- colSums(gexp != 0)

cdata <- data.frame(cellexp)
ggplot(data = cdata, mapping = aes(x=log10(cellexp))) +
  geom_histogram(mapping = aes(y=stat(count/sum(count))), 
                 binwidth = 1, bins = 4, boundary=0, closed="right", 
                 color="black", fill="white") +
  labs(title="Prevalence of Genes by Orders of Magnitude", x="log10(# of cells expressing a gene)", y = "proportion of genes")

head(cellexp[order(cellexp, decreasing = TRUE)])

#CPM normalize
numgenes <- rowSums(gexp)
normgexp <- gexp/numgenes*1e6

#add pseudocounts to do log
mat <- log10(normgexp+1)

# histogram of Gad1 expression among cells
mydf <- data.frame(mat)
#colors <- c("blue", rep("blue",19))
phist <- ggplot(data = mydf,
       mapping = aes(x = Gad1)) +
  geom_histogram(mapping = aes(fill = mydf['Gad1'] > 0), bins = 20) +
  scale_fill_manual("Gad1>0",values = c("blue", "red"))		 
phist


###############
#PCA
pcs <-prcomp(mat)
save(pcs, file="MERFISH_pcs_5000.RData")
df <- data.frame(x=c(1:30), y=pcs$sdev[1:30])
p <- ggplot(data = df,mapping = aes(x=x,y=y) ) + geom_line() +
  labs(title="Principle Components", x="index", y = "standard deviation")
p

###############
# tSNE
set.seed(0) 
emb <- Rtsne(pcs$x[,1:2], dims=2, perplexity = 30)$Y
rownames(emb) <- rownames(mat)
head(emb)
dim(emb)
df1 <- data.frame(x=emb[,1],
                 y=emb[,2],
                 col = mat[,'Gad1']) 

p1 <- ggplot(data=df1, mapping = aes(x=x, y=y)) +
  geom_point(mapping = aes(col = col), size=1) + 
  theme_classic() + 
  scale_color_gradientn("Gad1", colours =
                          c("blue", "white","red"), values = c(0,0.01,1)) +
  labs(title="tSNE of PC1 to PC2")
p1

emb2 <- Rtsne(pcs$x[,1:5], dims=2, perplexity = 30)$Y
rownames(emb2) <- rownames(mat)
head(emb2)
dim(emb2)
df2 <- data.frame(x=emb2[,1],
                 y=emb2[,2],
                 col = mat[,'Gad1']) 

p2 <- ggplot(data=df2, mapping = aes(x=x, y=y)) +
  geom_point(mapping = aes(col = col), size=1) + 
  theme_classic() + 
  scale_color_gradientn("Gad1", colours =
                          c("blue", "white","red"), values = c(0,0.01,1)) +
  labs(title="tSNE of PC1 to PC5")

p2


emb3 <- Rtsne(pcs$x[,1:10], dims=2, perplexity = 30)$Y
rownames(emb3) <- rownames(mat)
head(emb3)
dim(emb3)
df3 <- data.frame(x=emb3[,1],
                  y=emb3[,2],
                  col = mat[,'Gad1']) 

p3 <- ggplot(data=df3, mapping = aes(x=x, y=y)) +
  geom_point(mapping = aes(col = col), size=1) + 
  theme_classic() + 
  scale_color_gradientn("Gad1", colours =
                          c("blue", "white","red"), values = c(0,0.01,1)) +
  labs(title="tSNE of PC1 to PC10")

p3

emb4 <- Rtsne(pcs$x[,1:20], dims=2, perplexity = 30)$Y
rownames(emb4) <- rownames(mat)
head(emb4)
dim(emb4)
df4 <- data.frame(x=emb4[,1],
                  y=emb4[,2],
                  col = mat[,'Gad1']) 

p4 <- ggplot(data=df4, mapping = aes(x=x, y=y)) +
  geom_point(mapping = aes(col = col), size=1) + 
  theme_classic() + 
  scale_color_gradientn("Gad1", colours =
                          c("blue", "white","red"), values = c(0,0.01,1)) +
  labs(title="tSNE of PC1 to PC20")

p4

emb5 <- Rtsne(pcs$x[,1:30], dims=2, perplexity = 30)$Y
rownames(emb5) <- rownames(mat)
head(emb5)
dim(emb5)
df5 <- data.frame(x=emb5[,1],
                  y=emb5[,2],
                  col = mat[,'Gad1']) 

p5 <- ggplot(data=df5, mapping = aes(x=x, y=y)) +
  geom_point(mapping = aes(col = col), size=1) + 
  theme_classic() + 
  scale_color_gradientn("Gad1", colours =
    c("blue", "white","red"), values = c(0,0.01,1)) +
  labs(title="tSNE of PC1 to PC30")

p5

## for arranging plots side by side
grid.arrange(p, phist, p1, p2, p3, p4, ncol=2)
ggsave("HW3.png")
