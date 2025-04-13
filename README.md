# Vic_OrthoANI
A collection of codes for OrthoANI distances estimations for fasta files
ORHTO-ANI as been adapted from https://github.com/althonos/orthoani

![distances plot](https://github.com/user-attachments/assets/24995460-984e-4584-8e4c-b6ac3df79851)

# step 1 : generate a file containing pairwise combinantions
```r
#!/bin/bash
input=$1
ls *.fasta > sequences.txt
ruby -0777 -F'\n' -lane '$F.combination(2) { |c| puts c.join("@")}' sequences.txt > ${1}.txt
cat ${1}.txt
```
# step 2 : estimate OrthoANI distances 
```r
#!/bin/bash
input=$1
output=$2
for line in $(cat $input)
do
genome1=$(cut -f1 -d'@' <<< "$line")
genome2=$(cut -f2 -d'@' <<< "$line")
res=$(orthoani -q $genome1 -r $genome2 -j 20)
echo $genome1 $genome2 $res
done > a.txt
sed '1i genome_1 genome_2 distance' a.txt | tr " " "\t" > ${output}.tsv ;
rm a.txt ; 
cat ${output}.tsv ;
ls -lh ;
```
# usage: 
```r
./command0.sh input_name_prefix
./command1.sh input_name_prefix.txt output_name_prefix
```

# step3: Very nice R-plot:
```r
![hextable_colors](https://github.com/user-attachments/assets/f80f07d3-453f-4d1b-81d6-e2800200c345)

# 1. read the table 
data <- read.csv("distances_1.tsv", header=T, sep="\t")
data$distance <- as.numeric(data$distance)

# 2. load functions 
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

# 3. transform the table into matrix 

# install.packages("igraph") #
library(igraph)

pairwise <- data
head(pairwise)
g1 <- graph.data.frame(pairwise, directed=FALSE)
matrix2 <- get.adjacency(g1, attr="distance", sparse=FALSE)
head(matrix2)

# 4. load function to reorder the comparisons and observe patterns 

reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# 5. generate the final matrix
 
cormat <- reorder_cormat(matrix2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# 6. PLOT
 
p2 <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white") +
  geom_text(aes(label = round(value,2)), color = "#333333", size = 3) +
  scale_fill_gradientn(colours = c("blue","green","red"), 
                       values = rescale(c(89,95,100)),
                       guide = "colorbar", limits=c(89,100),
                       name="Ortho-ANI\ndistance") +

  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name="") + scale_y_discrete(name="") + coord_fixed() 
p2
```


