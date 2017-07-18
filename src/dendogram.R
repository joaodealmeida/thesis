#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(vegan)
library(ape)
library(geiger)

if (length(args) < 1) {
  stop("./Rscript files_directory", call.=FALSE)
} else if(length(args) > 2) {
  stop("./Rscript files_directory (ONLY 1 ARGUMENTS SHOULD BE PASSED)", call.=FALSE)
}
directory <- args[1]


setwd(directory)
teste <- read.csv("heatmap.csv")
clust_teste <- vegdist(teste, method = "jaccard")

png("dendogram.png",    # create PNG for the heat map
    width = 10*300,        # 5 x 300 pixels
    height = 10*200,
    res = 300)            # 300 pixels per inch

hc <- hclust(clust_teste, method="average")

neighbor_joining <- nj(clust_teste)


#dend <- as.dendrogram(hc)
#Apply Colors to different groups
#cores <- read.csv("grupos_tecidos.csv")
#colors <- as.character(cores[,2])
#Sort colors based on dendrogram
#colors <- colors[order.dendrogram(dend)]

#Apply to dendrogram
#labels_colors(dend, labels = FALSE) <- colors


#hang.dendrogram(dend, hang = -1)
#par(mar=c(16.1,4.1,0,2.1))
#plot(dend)
#plot(tr, edge.color = $)
#boot <- boot.phylo(tr, teste, function(xx) nj(clust_teste))

#plot(tr)

#nodelabels()
#dev.off()

png("dendogram_NJ.png",    # create PNG for the heat map
    width = 10*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300)            # 300 pixels per inch
exo <- neighbor_joining$tip.label[c(1,2,3,20,30,33,36,38,39,45,46)]
cardio <- neighbor_joining$tip.label[c(4,5,6,28,29)]
nervous <- neighbor_joining$tip.label[c(7,8,9,10,11,12,13,14,15,16,17,18,19,35)]
integu <- neighbor_joining$tip.label[c(21,22,40,41)]
digestive <- neighbor_joining$tip.label[c(23,24,25,26,27,31,37,42,44)]
urogenital <- neighbor_joining$tip.label[c(47,48)]
hemic <- neighbor_joining$tip.label[c(43,49)]
respiratory <- neighbor_joining$tip.label[c(32)]

plot(neighbor_joining, tip.color = ifelse(neighbor_joining$tip.label %in% exo, "chartreuse3",
                            ifelse(neighbor_joining$tip.label %in% cardio, "red",
                                   ifelse(neighbor_joining$tip.label %in% nervous, "purple",
                                          ifelse(neighbor_joining$tip.label %in% integu, "blue",
                                                 ifelse(neighbor_joining$tip.label %in% digestive, "brown",
                                                        ifelse(neighbor_joining$tip.label %in% urogenital, "orange",
                                                               ifelse(neighbor_joining$tip.label %in% hemic, "peachpuff4",
                                                                      ifelse(neighbor_joining$tip.label %in% respiratory, "cyan4","black")))))))))

nodelabels(frame = "circle", cex = 0.7 , text = , bg = "white")
dev.off()

print(neighbor_joining);
print("EXPORTING NEIGHBOR JOINING TREE...")
dir.create(file.path(getwd(), "trees"))
write.tree(neighbor_joining, file ="trees/neighbor_joining.newick")
upgma <- as.phylo(hc)


png("dendogram_UPGMA.png",    # create PNG for the heat map        
    width = 10*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300)            # 300 pixels per inch

plot(upgma, tip.color = ifelse(upgma$tip.label %in% exo, "chartreuse3",
                            ifelse(upgma$tip.label %in% cardio, "red",
                                   ifelse(upgma$tip.label %in% nervous, "purple",
                                          ifelse(upgma$tip.label %in% integu, "blue",
                                                 ifelse(upgma$tip.label %in% digestive, "brown",
                                                        ifelse(upgma$tip.label %in% urogenital, "orange",
                                                               ifelse(upgma$tip.label %in% hemic, "peachpuff4",
                                                                      ifelse(upgma$tip.label %in% respiratory, "cyan4","black")))))))))

nodelabels(frame = "circle", cex = 0.7 , text = , bg = "white")
dev.off()


print("EXPORTING UPGMA TREE...")
write.tree(upgma, file= "trees/upgma.newick")
