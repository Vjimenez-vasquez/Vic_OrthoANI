setwd("/home/hp/Documentos/NGS_bartonella/rochalimae/ortho-ANI/rochalimae_genomes_new")
dir()

data <- read.csv("distances_1.tsv", header=T, sep="\t")
data$distance <- as.numeric(data$distance)

head(data)
genomes <- unique(c(data$genome_1,data$genome_2))

data1 <- data[!data$genome_1 %in% genomes[c(5,14,15)], ]
data2 <- data1[!data1$genome_2 %in% genomes[c(5,14,15)], ]
dim(data2)
data <- data2

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# install.packages("igraph") #
library(igraph)

pairwise <- data
head(pairwise)
g1 <- graph.data.frame(pairwise, directed=FALSE)
matrix2 <- get.adjacency(g1, attr="distance", sparse=FALSE)
head(matrix2)

reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

cormat <- reorder_cormat(matrix2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

colors <- c("#33ccff","#00ff00","#ff0000")
p2 <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white") +
  geom_text(aes(label = round(value,2)), color = "#333333", size = 3) +
  scale_fill_gradientn(colours = colors, 
                       values = rescale(c(89,95,100)),
                       guide = "colorbar", limits=c(89,100),
                       name="OrthoANI\ndistance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name="") + scale_y_discrete(name="") + coord_fixed() 
p2