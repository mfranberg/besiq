#!/bin/bash -l

module add python/2.7.4

if [[ $# != 3 ]]; then
    echo "Usage: env_single model std output_dir"
    exit 1
fi

model=$1
std=$2
output_dir=$3

python env.py --model $model --std $std -n 100 $output_dir
