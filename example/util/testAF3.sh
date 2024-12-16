#!/bin/bash

#SBATCH -J testAF3
#SBATCH -o log/testAF3.out
#SBATCH -p gpu
#SBATCH -w devrtx
#SBATCH --gres=gpu:1


GPUS=1
MODEL_PARAMETERS_DIR=/fast/databases/alphafold3/weights
DATABASES_DIR=/fast/databases/alphafold3

INPUTDIR=/home/seamustard52/repository/colabfold-rachelse/example/input_json/yoshi
OUTPUTDIR=/home/seamustard52/repository/colabfold-rachelse/example/output/yoshi

docker run -it \
    --volume $INPUTDIR:/root/af_input \
    --volume $OUTPUTDIR:/root/af_output \
    --volume $MODEL_PARAMETERS_DIR:/root/models \
    --volume $DATABASES_DIR:/root/public_databases \
    --gpus $GPUS \
    alphafold3 \
    python run_alphafold.py \
    --json_path=/root/af_input \
    --model_dir=/root/models \
    --output_dir=/root/af_output

#docker run -it --volume /home/seamustard52/repository/colabfold-rachelse/example/input_json/yoshi/:/root/af_input --volume /home/seamustard52/repository/colabfold-rachelse/example/output/yoshi/:/root/af_output --volume /fast/databases/alphafold3/weights/:/root/models --volume /fast/databases/alphafold3/:/root/public_databases --gpus 1 alphafold3 python run_alphafold.py --json_path=/root/af_input/cf_hetero2mer_P01241_P10912_711a7.json  --model_dir=/root/models --output_dir=/root/af_output
# devrtx: needs cuda update
# devbox002: worked

#TODO: how to do batch process?
#TODO: how to use docker in slurm?
#TODO: update cuda on devrtx