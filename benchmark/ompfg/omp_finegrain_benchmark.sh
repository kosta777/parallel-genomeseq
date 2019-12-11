#!/bin/sh

benchmark_dir=$(dirname `realpath $0`)
project_dir=$(dirname $benchmark_dir)

#args:
code_section="solve_small"
n_reads="10"

finegrain_type="1"
#timing_file_path=$project_dir"/data/timings/timing_ompfg_0612.csv"
timing_file_path=$project_dir"/data/timings/ompfg_timing_results_1112.csv"
ref_file_path=$project_dir"/data/custom_ref_1.fa"
reads_file_path=$project_dir"/data/custom_reads_1.csv"

while true
do
  for n_threads in 1 2 4 5 6 7 8 10 16 32 64
  do
    eval "export OMP_NUM_THREADS="$n_threads
    echo "sh updated OMP_NUM_THREADS:"
    printenv OMP_NUM_THREADS
    echo "executing binary..."
  #  bin/omp_sw_solve_small solve_small 10 n_threads 1
    cmd=$project_dir"/bin/omp_sw_solve_small "$code_section" "$n_reads" "$n_threads" "$finegrain_type" "$timing_file_path" "$ref_file_path" "$reads_file_path
    echo $cmd
    eval $cmd
  done
done
