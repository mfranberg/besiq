#!/bin/bash

if [[ $# != 4 ]]; then
    echo "Usage: run_bayesic.sh pair_file plink_prefix pheno_file output_dir"
    exit 1
fi

pair_file=$1
plink_prefix=$2
pheno_file=$3
output_dir=$4

pair_name=$(basename $pair_file)

bayesic -p $pheno_file -m lm-stepwise $pair_file $plink_prefix > $output_dir/$pair_name.res
