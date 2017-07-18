#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

setwd("~/Desktop/Dissertacao/sol/data/mitocondrias/feup-pp/project/data/output/camacho-ficheiros")

#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

data <- read.csv("occurrences.csv")
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
mat_data[lower.tri(mat_data)] <- NA
rownames(mat_data) <- rnames                  # assign row names


#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)



# creates a 5 x 5 inch image
png("heatmap.png",    # create PNG for the heat map        
    width = 15*300,        # 5 x 300 pixels
    height = 15*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size



heatmap.2(mat_data,
          cellnote = mat_data,  # same data set for cell labels
          main = "GO_BP 08 Pathway Ocurrences", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(26,26),     # widens margins around plot
          cexRow = 1.5,
          cexCol = 1.5,
          col=my_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Colv="NA",
          Rowv = "NA",
          lwid = c(0.1,4),
          lhei = c(0.1,4))            # turn off column clustering

dev.off()               # close the PNG device

# creates a 5 x 5 inch image
png("heatmap_clu.png",    # create PNG for the heat map        
    width = 15*300,        # 5 x 300 pixels
    height = 15*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(mat_data,
          cellnote = mat_data,  # same data set for cell labels
          main = "GO_BP 08 Pathway Ocurrences", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(26,26),     # widens margins around plot
          cexRow = 1.5,
          cexCol = 1.5,
          col=my_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Colv="NA",
          lwid = c(0.1,4),
          lhei = c(0.1,4))            # turn off column clustering

dev.off()               # close the PNG device
