#What I want to make salient with this visualization is:
#what is the chance that if I pick a gene at random, it will be expressed in many cells?

#Load MERFISH data
data <- read.csv('~/Dropbox/JHU/Courses/genomic-data-visualization/data/MERFISH_Slice2Replicate2_halfcortex.csv.gz')

#make matrix of gene expression data without position
gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]

# for each gene, count how many cells express that gene
cellexp <- colSums(gexp != 0)

# Make histogram of prevalence of genes
# to make more salient visualization, 
# plot x-axis in log10 to use orders of magnitude
# plot y-axis in proportion 
library(ggplot2)
df <- data.frame(cellexp)
ggplot(data = df, mapping = aes(x=log10(cellexp))) +
  geom_histogram(mapping = aes(y=stat(count/sum(count))), 
                 binwidth = 1, bins = 4, boundary=0, closed="right", 
                 color="black", fill="white") +
  labs(title="Prevalence of Genes by Orders of Magnitude", x="log10(# of cells expressing a gene)", y = "proportion of genes")

#save as png 
ggsave("gene_prevalence.png")
