cd ../../tool/



DS=(IMDB DBLP FourSquare)
# DS=(DBLP FourSquare)
L=(2 3 4 5 6 7)

for ds in ${DS[@]}; do
    g++ -std=c++17 -O3 -g -fopenmp ./generate_${ds}.cpp -o ./generate_${ds}
    for len in ${L[@]}; do
        ./generate_${ds} ${len}
    done
done