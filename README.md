# FEUP - Thesis

## Co-expression networks between protein encoding mitochondrial genes and all the remaining genes in human tissues

#### Compile

```
g++ parser.cpp -o parser.out -Wall -std=c++11
g++ correlations_prot_filter.cpp -o correlations_prot_filter.out -Wall -std=c++11
g++ count_pathways.cpp -o count_pathways.out -Wall -std=c++11
g++ filter_camacho.cpp -o filter_camacho.out -Wall -std=c++11
g++ filter_data.cpp -o filter_data.out -Wall -std=c++11
g++ goterms.cpp -o goterms.out -Wall -std=c++11
g++ pathways.cpp -o pathways.out -Wall -std=c++11
g++ prepend_tissue_name.cpp -o prepend_tissue_name.out -Wall -std=c++11
```

#### Run

```
./parser.out gtex/rna-seq-data/All_Tissue_Site_Details_Analysis.combined.rpkm.gct gtex/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt
sh runAllCorrelactionsParallel.sh
sh runAllEnrichment.sh
```
```
Import each "-enriched" tissue to cytoscape. (Settings: Source Node | Target Node | Edge Attribute | Source Node Attribute (List of string) | Target Node Attribute (List of string)

HOW TO IMPORT ? - https://giphy.com/gifs/3ohz6wh5m4Huz5rkn6/fullscreen


Go to Tools -> Network Analyzer -> Network Analysis -> Analyze Network -> Treat the network as undirected.
(Wait)

HOW TO ANALYZE NETWORK ? - https://giphy.com/gifs/3ohz6Bp0yR5mS0XtxC/fullscreen


Export (default edge) to data/output/cytoscape-networks with the name tissuename-edges.csv
Export (default node) to data/output/cytoscape-networks with the name tissuename-nodes.csv

HOW TO EXPORT THE NETWORK ? - https://giphy.com/gifs/3oAkalAYrlnx1Nffva

Example of file names: Adipose - Visceral (Omentum)-edges.csv and Adipose - Visceral (Omentum)-nodes.csv
```
```
When every tissue has its own -edges/nodes cytoscape files in the folder, run:
sh makeDendograms.sh

Use the /data/output/camacho-ficheiros/tissueInfo.txt and the tree in newick format (/data/output/camacho-ficheiros/trees) and import them into BioTree Viewer for better analysis.
(https://joaoalmeida.me/dissertacao/viewer/) or https://github.com/joaodealmeida/biotreeviewer if you want to run it on your local server.

If you have need help, contact me: contact at joaoalmeida dot me.
```
