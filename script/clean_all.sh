#!/bin/bash

cur_dir=${PWD}

cd ${cur_dir}

DS=(IMDB FourSquare DBLP)

for ds in ${DS[@]}; do
    pwd
    cd ${cur_dir}
    cd ../materials/input/metapath/${ds}/
    rm -f *
done

for ds in ${DS[@]}; do
    pwd
    cd ${cur_dir}
    cd ../materials/output/efficient/$ds/
    rm -f *

    cd ${cur_dir}
    cd ../materials/output/accuracy/$ds/
    rm -f *
done
