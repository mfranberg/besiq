#!/bin/bash -l

module add R/2.15.2

script_dir=$(dirname $0)

for power_path in `ls paper/maf20s22/power/*.out`;
do
    power_type=$(basename $power_path);
    
    touch paper/$power_type
    rm paper/$power_type

    cat paper/maf20s22/power/$power_type | awk '{ print $0, "0.2", "2000" }' >> paper/$power_type
    cat paper/maf20s33/power/$power_type | awk '{ print $0, "0.2", "3000" }' >> paper/$power_type
    cat paper/maf20s44/power/$power_type | awk '{ print $0, "0.2", "4000" }' >> paper/$power_type
    cat paper/maf30s22/power/$power_type | awk '{ print $0, "0.3", "2000" }' >> paper/$power_type
    cat paper/maf30s33/power/$power_type | awk '{ print $0, "0.3", "3000" }' >> paper/$power_type
    cat paper/maf30s44/power/$power_type | awk '{ print $0, "0.3", "4000" }' >> paper/$power_type
    cat paper/maf40s22/power/$power_type | awk '{ print $0, "0.4", "2000" }' >> paper/$power_type
    cat paper/maf40s33/power/$power_type | awk '{ print $0, "0.4", "3000" }' >> paper/$power_type
    cat paper/maf40s44/power/$power_type | awk '{ print $0, "0.4", "4000" }' >> paper/$power_type

    Rscript $script_dir/../run_experiment/explib/external/plot_joint_power.r "Heritability" "Power" X paper/$power_type paper/$power_type.pdf
done

