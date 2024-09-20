cd ../../accuracy/

algs=("dm" "lg" "meta" "mnc" "res")

# algs=("mnc" "res")

t=${1}

for alg in ${algs[@]}; do
    g++ -std=c++17 -O3 -g -fopenmp ./${alg}.cpp -o ./${alg}
    
    ./${alg} \
    ${t} \
    1 \
    ../materials/input/graph/IMDB/graph.txt \
    ../materials/input/graph/IMDB/vertex.txt \
    ../materials/input/graph/IMDB/edge.txt \
    ../materials/input/metapath/IMDB/IMDB.txt \
    ../materials/output/accuracy/IMDB/${alg}_acc.txt
done