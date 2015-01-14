#!/bin/bash -l

module add python/2.7.4

if [[ $# != 2 ]]; then
    echo "Usage: env_power allocation output_dir"
    exit 1
fi

allocation=$1
output_dir=$2

script_dir=$(dirname $0)

model_id=0
for model in {"1.1 1.2 1.3 1.1 1.0 0.9","1.0 1.0 1.0 1.0 1.2 1.3"};
do
    for num_snps in {10000,1000000,10000000};
    do
        batch_output_dir=$output_dir/model$model_id/snps$num_snps/
        mkdir -p $batch_output_dir

        sbatch -A $allocation -t 2:00:00 -p core -n 1 $script_dir/env_power_single.sh "$model" 1.0 $num_snps $batch_output_dir
    done

    model_id=$((model_id+1))
done
