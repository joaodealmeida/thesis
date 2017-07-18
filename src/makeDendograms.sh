./runAllCamacho.sh
./runAllPrependTissue.sh
cd ../data/output/camacho-ficheiros/
cat *_*.txt > tissueInfo.txt
cd ../../../src/
./count_pathways.o "../data/output/camacho-ficheiros/tissueInfo.txt"
./dendogram.R "../data/output/camacho-ficheiros"
echo "Make Dendograms Finished, check /tree folders and graphic images generated."
