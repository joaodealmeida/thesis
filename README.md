# FEUP - Thesis

## Identificação de grupos de co-expressão entre todos os genes de proteínas mitocondriais

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

Import each "-enriched" tissue to cytoscape. (Settings: Source Node | Target Node | Edge Attribute | Source Node Attribute (List of string) | Target Node Attribute (List of string)

(http://i.imgur.com/IOTRWTi.png)


Go to Tools -> Network Analyzer -> Network Analysis -> Analyze Network -> Treat the network as undirected.
(Wait)
Export (default edge) to data/output/cytoscape-networks with the name tissuename-edges.csv
Export (default node) to data/output/cytoscape-networks with the name tissuename-nodes.csv

Example of file names: Adipose - Visceral (Omentum)-edges.csv and Adipose - Visceral (Omentum)-nodes.csv

When every tissue has its own -edges/nodes cytoscape files in the folder, run:
sh makeDendograms.sh

Use the /data/output/camacho-ficheiros/tissueInfo.txt and the tree in newick format (/data/output/camacho-ficheiros/trees) in the BioTree Viewer for better analysis.
(https://joaoalmeida.me/dissertacao/viewer/)
```
