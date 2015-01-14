#!/bin/bash

if [[ $# != 3 ]]; then
    echo "Usage: run_bayesic.sh pair_file plink_prefix output_dir"
    exit 1
fi

pair_file=$1
plink_prefix=$2
output_dir=$3

pair_name=$(basename $pair_file)

#bayesic -m stepwise $pair_file $plink_prefix > $output_dir/$pair_name.res
bayesic -m wald $pair_file $plink_prefix > $output_dir/$pair_name.res
