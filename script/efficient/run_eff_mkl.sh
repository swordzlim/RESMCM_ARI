cd ../../efficient/

mkl_root=${1}
t=64
len=5
DS=(IMDB FourSquare DBLP)

algs=("mkl_csr" "mkl_csc" "mkl_bsr")
# algs=("mkl_csr" "mkl_csc" "mkl_bsr")

. ${mkl_root}/env/vars.sh

for ds in ${DS[@]}; do
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
        ../materials/input/graph/$ds/graph.txt \
        ../materials/input/graph/$ds/vertex.txt \
        ../materials/input/graph/$ds/edge.txt \
        ../materials/input/metapath/$ds/${ds}-l${len}.txt \
        ../materials/output/efficient/$ds/${alg}_eff-l${len}-t${t}.txt
    done
done