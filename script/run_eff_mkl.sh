for arg in "$@"
do
  case $arg in
    -mkl_root=*)
    mkl_root="${arg#*=}"
    shift
    ;;
    *)
    ;;
  esac
done

cur_dir=${PWD}

# run efficient mkl
cd $cur_dir/efficient
bash ./run_eff_mkl.sh ${mkl_root}


