#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ape)
library(geiger)

if (length(args) < 3) {
  stop("./Rscript directory tree_name node_number", call.=FALSE)
} else if(length(args) > 4) {
  stop("./Rscript directory tree_name node_numbes (ONLY 3 ARGUMENTS SHOULD BE PASSED)", call.=FALSE)
}

setwd(args[1])
tree_name <- paste(getwd(),"/",args[2], sep ="")
nodeN <- args[3]

print("READING TREE INFO...")

dendrogram <- read.tree(file = tree_name)

print("GETTING NODE DESCENDANTS...")
descendants <- tips(dendrogram, node = nodeN)

#REMOVING "_" FROM SOME TISSUES
for(i in 1:length(descendants)){
  ifelse(substr(descendants[i],nchar(descendants[i]),nchar(descendants[i])) == "_", descendants[i] <- substr(descendants[i], 1, nchar(descendants[i])-1), descendants[i])
}

print("EXPORTING TISSUES to /clusters FOLDER...")

dir.create(file.path(getwd(), "clusters"))
cat(descendants, sep = "\n", file = paste("clusters/cluster", nodeN, sep = "-"))