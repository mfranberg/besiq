#!/bin/bash -l

module add python/2.7.4

if [[ $# != 2 ]]; then
    echo "Usage: env_multiple allocation output_dir"
    exit 1
fi

allocation=$1
output_dir=$2

script_dir=$(dirname $0)

for i in {1..10};
do
    batch_output_dir=$output_dir/batch$i/
    mkdir -p $batch_output_dir

    sbatch -A $allocation -t 8:00:00 -p core -n 1 $script_dir/env_single.sh "1.1 1.2 1.3 1.2 1.3 1.4" 1.0 $batch_output_dir
done
