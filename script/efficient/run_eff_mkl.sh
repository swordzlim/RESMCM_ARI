cd ../../efficient/

t=${1}
mkl_root=${2}

algs=("mkl_csr" "mkl_csc" "mkl_bsr" "rescsr" "resmcm")
# algs=("mkl_csr" "mkl_csc" "mkl_bsr")

. ${mkl_root}/env/vars.sh

for alg in ${algs[@]}; do
    g++ -std=c++17 -O3 -g -fopenmp ./${alg}.cpp -o ./${alg} \
    -I${mkl_root}/include \
    -L${mkl_root}/lib/intel64 \
    -Wl,--no-as-needed \
    -lmkl_intel_lp64 \
    -lmkl_gnu_thread \
    -lmkl_sequential \
    -lmkl_core \
    -lgomp \
    -lpthread \
    -lm \
    -ldl
    
    ./${alg} \
    ${t} \
    1 \
    ../materials/input/graph/IMDB/graph.txt \
    ../materials/input/graph/IMDB/vertex.txt \
    ../materials/input/graph/IMDB/edge.txt \
    ../materials/input/metapath/IMDB/IMDB.txt \
    ../materials/output/efficient/IMDB/${alg}_eff.txt
done