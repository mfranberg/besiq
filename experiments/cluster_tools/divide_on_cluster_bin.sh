#!/bin/bash

if [[ $# != 4 ]]; then
    echo "Usage: divide_on_cluster.sh num_cores pairs plink_prefix output_dir"
    exit 1
fi

num_cores=$1
pairs_file=$2
plink_file=$3
output_dir=$4

script_root=$(dirname $0)

tmp_dir="$output_dir/tmp"
mkdir -p $tmp_dir
split_prefix="$tmp_dir/split_pair"

num_pairs=$(wc -l $pairs_file | awk '{ print $1}')
pairs_per_node=$((num_cores + num_pairs / num_cores))

split -l $pairs_per_node -d $pairs_file $split_prefix
for split_pair_file in `ls $split_prefix*`;
do
    sbatch -A b2012114 -t 8:00:00 -p core -n 1 $script_root/run_bayesic_bin.sh $split_pair_file $plink_file $output_dir
done
