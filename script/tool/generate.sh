len=${1}

cd ../../tool/

g++ -std=c++17 -O3 -g -fopenmp ./generate.cpp -o ./generate

./generate ${len}