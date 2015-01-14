#!/bin/bash -l

module add python/2.7.4

python $script_dir/../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 $script_dir/../run_experiment/types/null/flat.example.json paper/null_flat/ &
python $script_dir/../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 $script_dir/../run_experiment/types/null/single.example.json paper/null_single/ &
python $script_dir/../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 $script_dir/../run_experiment/types/null/oddsadd.example.json paper/null_oddsadd/ &
python $script_dir/../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 $script_dir/../run_experiment/types/null/penadd.example.json paper/null_penadd/ &
python $script_dir/../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 $script_dir/../run_experiment/types/null/penmul.example.json paper/null_penmul/ &
python $script_dir/../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 $script_dir/../run_experiment/types/null/oddsmul.example.json paper/null_oddsmul/ &
python $script_dir/../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 $script_dir/../run_experiment/types/null/penhet.example.json paper/null_penhet/ &

wait
