#!/bin/bash
dir="/Users/joaoalmeida/Desktop/Dissertacao/sol/data/mitocondrias/feup-pp/project/data/output/camacho-ficheiros/trees/"
init=50
end=96
for ((i = $init; i <= $end; i++)); do
   ./dendogram_tissues_by_node.R $dir "neighbor_joining.newick" $i
done
