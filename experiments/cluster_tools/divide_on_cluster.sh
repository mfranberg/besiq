#!/bin/bash

if [[ $# != 5 ]]; then
    echo "Usage: divide_on_cluster.sh num_cores pairs plink_prefix pheno_file output_dir"
    exit 1
fi

num_cores=$1
pairs_file=$2
plink_file=$3
pheno_file=$4
output_dir=$5

script_root=$(dirname $0)

tmp_dir="$output_dir/tmp"
mkdir -p $tmp_dir
split_prefix="$tmp_dir/split_pair"

num_pairs=$(wc -l $pairs_file | awk '{ print $1}')
pairs_per_node=$((num_cores + num_pairs / num_cores))

split -l $pairs_per_node -d $pairs_file $split_prefix
for split_pair_file in `ls $split_prefix*`;
do
    sbatch -A b2013088 -t 8:00:00 -p core -n 1 $script_root/run_bayesic.sh $split_pair_file $plink_file $pheno_file $output_dir
done
