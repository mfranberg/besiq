#!/bin/bash -l

module add python/2.7.4
module add R/2.15.2

if [ $# -ne 2 ]
then
    echo "Usage: run_all.sh num_tests output_dir"
    exit 1
fi

num_tests=$1
root_dir=$2

python ../run_experiment/run_experiment.py --heritability 0.015 --closed-num-tests 1000000000000 2000000 2000000 1 -m 0.3 0.3 -n 2000 2000 -s 200 -i 1000000000000 ../run_experiment/types/all.example.json $root_dir/h10s22 &
python ../run_experiment/run_experiment.py --heritability 0.015 --closed-num-tests 1000000000000 2000000 2000000 1 -m 0.3 0.3 -n 3000 3000 -s 200 -i 1000000000000 ../run_experiment/types/all.example.json $root_dir/h10s33 &
python ../run_experiment/run_experiment.py --heritability 0.015 --closed-num-tests 1000000000000 2000000 2000000 1 -m 0.3 0.3 -n 4000 4000 -s 200 -i 1000000000000 ../run_experiment/types/all.example.json $root_dir/h10s44 &
python ../run_experiment/run_experiment.py --heritability 0.020 --closed-num-tests 1000000000000 2000000 2000000 1 -m 0.3 0.3 -n 2000 2000 -s 200 -i 1000000000000 ../run_experiment/types/all.example.json $root_dir/h20s22 &
python ../run_experiment/run_experiment.py --heritability 0.020 --closed-num-tests 1000000000000 2000000 2000000 1 -m 0.3 0.3 -n 3000 3000 -s 200 -i 1000000000000 ../run_experiment/types/all.example.json $root_dir/h20s33 &
python ../run_experiment/run_experiment.py --heritability 0.020 --closed-num-tests 1000000000000 2000000 2000000 1 -m 0.3 0.3 -n 4000 4000 -s 200 -i 1000000000000 ../run_experiment/types/all.example.json $root_dir/h20s44 &
python ../run_experiment/run_experiment.py --heritability 0.025 --closed-num-tests 1000000000000 2000000 2000000 1 -m 0.3 0.3 -n 2000 2000 -s 200 -i 1000000000000 ../run_experiment/types/all.example.json $root_dir/h30s22 &
python ../run_experiment/run_experiment.py --heritability 0.025 --closed-num-tests 1000000000000 2000000 2000000 1 -m 0.3 0.3 -n 3000 3000 -s 200 -i 1000000000000 ../run_experiment/types/all.example.json $root_dir/h30s33 &
python ../run_experiment/run_experiment.py --heritability 0.025 --closed-num-tests 1000000000000 2000000 2000000 1 -m 0.3 0.3 -n 4000 4000 -s 200 -i 1000000000000 ../run_experiment/types/all.example.json $root_dir/h30s44 &

wait

touch $root_dir/combined_power.out
rm $root_dir/combined_power.out

cat $root_dir/h10s22/power/power_all.out | awk '{ print $0, 0.3, "2000" }' >> $root_dir/combined_power.out
cat $root_dir/h10s33/power/power_all.out | awk '{ print $0, 0.3, "3000" }' >> $root_dir/combined_power.out
cat $root_dir/h10s44/power/power_all.out | awk '{ print $0, 0.3, "4000" }' >> $root_dir/combined_power.out
cat $root_dir/h20s22/power/power_all.out | awk '{ print $0, 0.3, "2000" }' >> $root_dir/combined_power.out
cat $root_dir/h20s33/power/power_all.out | awk '{ print $0, 0.3, "3000" }' >> $root_dir/combined_power.out
cat $root_dir/h20s44/power/power_all.out | awk '{ print $0, 0.3, "4000" }' >> $root_dir/combined_power.out
cat $root_dir/h30s22/power/power_all.out | awk '{ print $0, 0.3, "2000" }' >> $root_dir/combined_power.out
cat $root_dir/h30s33/power/power_all.out | awk '{ print $0, 0.3, "3000" }' >> $root_dir/combined_power.out
cat $root_dir/h30s44/power/power_all.out | awk '{ print $0, 0.3, "4000" }' >> $root_dir/combined_power.out

Rscript ../run_experiment/explib/external/plot_all_joint.r "Power (t)" "Fraction of models with power >= t" X $root_dir/combined_power.out $root_dir/combined_all.pdf
