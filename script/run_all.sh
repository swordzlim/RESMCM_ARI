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

# run efficient first
cd $cur_dir/efficient
bash ./run_eff.sh ${thread_num}
# draw
cd $cur_dir/efficient
python ./draw.py ${len}

# run accuracy second
cd $cur_dir/accuracy
bash ./run_acc.sh ${thread_num}
# draw
cd $cur_dir/accuracy
python ./draw.py ${len}