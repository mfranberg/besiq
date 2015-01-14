#!/bin/bash -l

module add python/2.7.4

python ../run_experiment/run_experiment.py -m 0.2 0.2 -n 2000 2000 -s 200 -i 1000000000000 --closed-num-tests 1000000000000 2000000 2000000 1 ../run_experiment/types/models.example.json paper/maf20s22 &
python ../run_experiment/run_experiment.py -m 0.2 0.2 -n 3000 3000 -s 200 -i 1000000000000 --closed-num-tests 1000000000000 2000000 2000000 1 ../run_experiment/types/models.example.json paper/maf20s33 &
python ../run_experiment/run_experiment.py -m 0.2 0.2 -n 4000 4000 -s 200 -i 1000000000000 --closed-num-tests 1000000000000 2000000 2000000 1 ../run_experiment/types/models.example.json paper/maf20s44 &
python ../run_experiment/run_experiment.py -m 0.3 0.3 -n 2000 2000 -s 200 -i 1000000000000 --closed-num-tests 1000000000000 2000000 2000000 1 ../run_experiment/types/models.example.json paper/maf30s22 &
python ../run_experiment/run_experiment.py -m 0.3 0.3 -n 3000 3000 -s 200 -i 1000000000000 --closed-num-tests 1000000000000 2000000 2000000 1 ../run_experiment/types/models.example.json paper/maf30s33 &
python ../run_experiment/run_experiment.py -m 0.3 0.3 -n 4000 4000 -s 200 -i 1000000000000 --closed-num-tests 1000000000000 2000000 2000000 1 ../run_experiment/types/models.example.json paper/maf30s44 &
python ../run_experiment/run_experiment.py -m 0.4 0.4 -n 2000 2000 -s 200 -i 1000000000000 --closed-num-tests 1000000000000 2000000 2000000 1 ../run_experiment/types/models.example.json paper/maf40s22 &
python ../run_experiment/run_experiment.py -m 0.4 0.4 -n 3000 3000 -s 200 -i 1000000000000 --closed-num-tests 1000000000000 2000000 2000000 1 ../run_experiment/types/models.example.json paper/maf40s33 &
python ../run_experiment/run_experiment.py -m 0.4 0.4 -n 4000 4000 -s 200 -i 1000000000000 --closed-num-tests 1000000000000 2000000 2000000 1 ../run_experiment/types/models.example.json paper/maf40s44 &

wait
