#!/bin/bash -l

module add python/2.7.4

if [[ $# != 4 ]]; then
    echo "Usage: env_single model std num_snps output_dir"
    exit 1
fi

model=$1
std=$2
num_snps=$3
output_dir=$4

script_dir=$(dirname $0)

python $script_dir/../env/env.py --model $model --std $std -n 1 --num-variants 1000 --sample-size 4000 --num-tests $num_snps 1 --bonferroni $output_dir
