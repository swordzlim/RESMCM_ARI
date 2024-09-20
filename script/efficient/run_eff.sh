cd ../../efficient/

t=${1}

algs=("l2r" "naive" "dm" "lg" "meta" "mnc" "rescsr" "resmcm")
# algs=("resmcm")


for alg in ${algs[@]}; do
    g++ -std=c++17 -O3 -g -fopenmp ./${alg}.cpp -o ./${alg}
    
    ./${alg} \
    ${t} \
    1 \
    ../materials/input/graph/IMDB/graph.txt \
    ../materials/input/graph/IMDB/vertex.txt \
    ../materials/input/graph/IMDB/edge.txt \
    ../materials/input/metapath/IMDB/IMDB.txt \
    ../materials/output/efficient/IMDB/${alg}_eff.txt
done