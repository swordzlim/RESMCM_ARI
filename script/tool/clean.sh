#!/bin/bash

cur_dir=${PWD}

cd ${cur_dir}
pwd

DS=(IMDB DBLP FourSquare)

for ds in ${DS[@]}; do
    cd ${cur_dir}
    cd ../../materials/input/metapath/$ds/
    rm -f *

    cd ${cur_dir}
    cd ../../materials/output/efficient/$ds/
    rm -f *

    cd ${cur_dir}
    cd ../../materials/output/accuracy/$ds/
    rm -f *
done
