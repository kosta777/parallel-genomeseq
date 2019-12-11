#!/bin/sh

benchmark_dir=$(dirname `realpath $0`)
project_dir=$(dirname $benchmark_dir)
data_dir=$project_dir"/data"

ucsc_chr22_wget_cmd="wget \"http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz\" -P "$data_dir
eval $ucsc_chr22_wget_cmd
gunzip $data_dir"/chr22.fa.gz"


