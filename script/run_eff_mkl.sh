for arg in "$@"
do
  case $arg in
    -thread_num=*)
    thread_num="${arg#*=}"
    shift
    ;;
    -length=*)
    len="${arg#*=}"
    shift
    ;;
    -mkl_root=*)
    mkl_root="${arg#*=}"
    shift
    ;;
    *)
    ;;
  esac
done

echo "length: $len"
echo "thread_num: $thread_num"

cur_dir=${PWD}

cd $cur_dir/tool
bash ./clean.sh

cd $cur_dir/tool
bash ./generate.sh ${len}


# run efficient mkl
cd $cur_dir/efficient
bash ./run_eff_mkl.sh ${thread_num} ${mkl_root}



# draw
cd $cur_dir/efficient
python ./draw_mkl.py ${len}

