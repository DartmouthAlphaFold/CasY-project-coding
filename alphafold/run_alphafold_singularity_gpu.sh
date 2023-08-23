#!/bin/bash

#SBATCH --job-name=alphafold_2023-08-21_test232_params_v3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --partition=v100_12
#SBATCH --gres=gpu:v100:1
#SBATCH --mail-user=dnnguy23@colby.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/dartfs/rc/lab/H/HowellA/alphafold-log/%x.%j.out     
#SBATCH --error=/dartfs/rc/lab/H/HowellA/alphafold-log/%x.%j.err 

# View hostname and state of assigned GPUs
hostname
nvidia-smi

# Unload all modules
module purge

# Enable Multi-Process Service for GPU-based relaxation
# nvidia-cuda-mps-control -d

# Activate conda
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate alphafold-python

# Set environmental variables
export SINGULARITYENV_TF_FORCE_UNIFIED_MEMORY=1
export SINGULARITYENV_XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
export SINGULARITYENV_LD_LIBRARY_PATH=/optnfs/el7/cuda/11.2/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

# Output directory
OUTPUT_DIR=/dartfs/rc/lab/H/HowellA/alphafold-output/alphafold_2023-08-21_test232_params_v3

# Run with GPU
python3 /dartfs/rc/lab/H/HowellA/run_alphafold_singularity_gpu.py \
    --fasta_paths /dartfs/rc/lab/H/HowellA/faa_files/IFGSC_6mer.fasta \
    --model_preset multimer \
    --max_template_date 2018-01-01 \
    --output_dir $OUTPUT_DIR

# Visualize result
python3 /dartfs/rc/lab/H/HowellA/visualize_alphafold_results.py $OUTPUT_DIR

# Clean up
conda deactivate
# echo quit | nvidia-cuda-mps-control