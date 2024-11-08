cur_dir=${PWD}

cd $cur_dir/tool
bash ./clean.sh

cd $cur_dir/tool
bash ./generate_metapath.sh

cd $cur_dir
bash ./clean_all.sh

# run efficient first
cd $cur_dir/efficient
bash ./run_eff.sh

# run accuracy second
cd $cur_dir/accuracy
bash ./run_acc.sh
