# FEUP - Thesis

## Identificação de grupos de co-expressão entre todos os genes de proteínas mitocondriais

#### Compile

```
g++ parser.cpp -o parser.out -Wall -std=c++11
g++ correlations.cpp -o correlations.out -Wall -std=c++11
g++ pathways.cpp -o pathways.out -Wall -std=c++11
g++ filter_data.cpp -o filter_data.out -Wall -std=c++11
```

#### Run

```
./parser.out gtex/rna-seq-data/All_Tissue_Site_Details_Analysis.combined.rpkm.gct gtex/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt
sh runallCorrelactionsParallel.sh
sh runallFilter.sh
sh mergeOverlays.sh
sh runAllEnrichment.sh
./filter_data "TissueName" 1/2/3/4/5/6
```
