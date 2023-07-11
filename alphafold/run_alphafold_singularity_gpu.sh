#!/bin/bash

#SBATCH --job-name=alphafold_2023-04-29_CasY15
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=30:00:00
#SBATCH --partition=v100_12
#SBATCH --gres=gpu:v100:1
#SBATCH --mail-user=dnnguy23@colby.edu
#SBATCH --mail-type=ALL
#SBATCH --output=alphafold-log/%x.%j.out     
#SBATCH --error=alphafold-log/%x.%j.err 

# View hostname and state of assigned GPUs
hostname
nvidia-smi

# clean modules
module purge

# Move to lab folder
cd /dartfs/rc/lab/H/HowellA

# set env 
export SINGULARITYENV_TF_FORCE_UNIFIED_MEMORY=1
export SINGULARITYENV_XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
export SINGULARITYENV_LD_LIBRARY_PATH=/optnfs/el7/cuda/11.2/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

# Run with GPU
python3 run_alphafold_singularity_gpu.py \
    --fasta_paths /dartfs/rc/lab/H/HowellA/faa_files/CasY15.faa \
    --output_dir /dartfs/rc/lab/H/HowellA/alphafold-output/2023-04-29_CasY15

# Clean up
rm /dartfs/rc/lab/H/HowellA/ld.so.cache