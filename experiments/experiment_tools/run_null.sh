#!/bin/bash -l

module add python/2.7.4

python ../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 ../run_experiment/types/null/flat.example.json paper/null_flat/ &
python ../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 ../run_experiment/types/null/single.example.json paper/null_single/ &
python ../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 ../run_experiment/types/null/oddsadd.example.json paper/null_oddsadd/ &
python ../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 ../run_experiment/types/null/penadd.example.json paper/null_penadd/ &
python ../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 ../run_experiment/types/null/penmul.example.json paper/null_penmul/ &
python ../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 ../run_experiment/types/null/oddsmul.example.json paper/null_oddsmul/ &
python ../run_experiment/run_experiment.py -m 0.1 0.5 --sample-maf -n 3000 3000 -s 1000 -i 1000 ../run_experiment/types/null/penhet.example.json paper/null_penhet/ &

wait
