#!/bin/bash


START_SEED=$1
N_WORKERS=$2
N_SEEDS_EACH=$3 # Seeds each

for i in $(seq 0 $((N_WORKERS-1)))
do
    FROM=$((i * N_SEEDS_EACH + START_SEED))
    TO=$((i * N_SEEDS_EACH + N_SEEDS_EACH -1 + START_SEED))
    echo $FROM
    echo $TO
    Rscript run_opruner.R --seeds_from $FROM --seeds_to $TO --output_dir $4 --filename $5 &

done
