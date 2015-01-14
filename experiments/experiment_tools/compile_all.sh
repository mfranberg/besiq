#!/bin/bash

root_dir=~/proj/private/data/bayesic_all/1000000/

touch $root_dir/combined_power.out
rm $root_dir/combined_power.out

cat $root_dir/maf22s22/power/power_all.out | awk '{ print $0, 0.2, "2000" }' >> $root_dir/combined_power.out
cat $root_dir/maf22s33/power/power_all.out | awk '{ print $0, 0.2, "3000" }' >> $root_dir/combined_power.out
cat $root_dir/maf22s44/power/power_all.out | awk '{ print $0, 0.2, "4000" }' >> $root_dir/combined_power.out
cat $root_dir/maf33s22/power/power_all.out | awk '{ print $0, 0.3, "2000" }' >> $root_dir/combined_power.out
cat $root_dir/maf33s33/power/power_all.out | awk '{ print $0, 0.3, "3000" }' >> $root_dir/combined_power.out
cat $root_dir/maf33s44/power/power_all.out | awk '{ print $0, 0.3, "4000" }' >> $root_dir/combined_power.out
cat $root_dir/maf44s22/power/power_all.out | awk '{ print $0, 0.4, "2000" }' >> $root_dir/combined_power.out
cat $root_dir/maf44s33/power/power_all.out | awk '{ print $0, 0.4, "3000" }' >> $root_dir/combined_power.out
cat $root_dir/maf44s44/power/power_all.out | awk '{ print $0, 0.4, "4000" }' >> $root_dir/combined_power.out

Rscript ../run_experiment/explib/external/plot_all_joint.r "Power (t)" "Fraction of models with power >= t" X $root_dir/combined_power.out $root_dir/combined_all.pdf
