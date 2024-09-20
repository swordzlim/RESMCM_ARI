#!/bin/bash

cur_dir=${PWD}

cd ${cur_dir}
pwd

cd ${cur_dir}
cd ../materials/input/metapath/IMDB/
rm -f *

cd ${cur_dir}
cd ../materials/output/efficient/IMDB/
rm -f *

cd ${cur_dir}
cd ../materials/output/accuracy/IMDB/
rm -f *

# cd ${cur_dir}
# cd ../materials/output/fig/IMDB/
# rm -f *