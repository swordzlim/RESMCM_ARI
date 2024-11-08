cd ../../efficient/


algs=("l2r" "naive" "dm" "lg" "meta" "mnc" "rescsr" "resmcm")
DS=(IMDB FourSquare DBLP)
T=(64)
L=(2 3 4 5 6 7)

for ds in ${DS[@]}; do
    for t in ${T[@]}; do
        for len in ${L[@]}; do
            for alg in ${algs[@]}; do
                g++ -std=c++17 -O3 -g -fopenmp ./${alg}.cpp -o ./${alg}
                ./${alg} \
                ${t} \
                1 \
                ../materials/input/graph/${ds}/graph.txt \
                ../materials/input/graph/${ds}/vertex.txt \
                ../materials/input/graph/${ds}/edge.txt \
                ../materials/input/metapath/${ds}/${ds}-l${len}.txt \
                ../materials/output/efficient/${ds}/${alg}_eff-l${len}-t${t}.txt
            done
        done
    done
done

T=(8 16 32)
L=(7)

for ds in ${DS[@]}; do
    for t in ${T[@]}; do
        for len in ${L[@]}; do
            for alg in ${algs[@]}; do
                g++ -std=c++17 -O3 -g -fopenmp ./${alg}.cpp -o ./${alg}
                ./${alg} \
                ${t} \
                1 \
                ../materials/input/graph/${ds}/graph.txt \
                ../materials/input/graph/${ds}/vertex.txt \
                ../materials/input/graph/${ds}/edge.txt \
                ../materials/input/metapath/${ds}/${ds}-l${len}.txt \
                ../materials/output/efficient/${ds}/${alg}_eff-l${len}-t${t}.txt
            done
        done
    done
done