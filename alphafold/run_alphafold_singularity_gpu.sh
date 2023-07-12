#!/bin/bash

#SBATCH --job-name=alphafold_2023-07-11_mgnify-CasY
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=30:00:00
#SBATCH --partition=v100_12
#SBATCH --gres=gpu:v100:1
#SBATCH --mail-user=dnnguy23@colby.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/dartfs/rc/lab/H/HowellA/alphafold-log/%x.%j.out     
#SBATCH --error=/dartfs/rc/lab/H/HowellA/alphafold-log/%x.%j.err 

# View hostname and state of assigned GPUs
hostname
nvidia-smi

# clean modules
module purge

# activate conda
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate alphafold-python

# Move to submit folder
cd $SLURM_SUBMIT_DIR

# set env 
export SINGULARITYENV_TF_FORCE_UNIFIED_MEMORY=1
export SINGULARITYENV_XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
export SINGULARITYENV_LD_LIBRARY_PATH=/optnfs/el7/cuda/11.2/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

# output directory
OUTPUT_DIR=/dartfs/rc/lab/H/HowellA/alphafold-output/2023-07-11_mgnify-CasY

# Run with GPU
python3 /dartfs/rc/lab/H/HowellA/run_alphafold_singularity_gpu.py \
    --fasta_paths /dartfs/rc/lab/H/HowellA/faa_files/CasY_MGYP000178263316.fa,/dartfs/rc/lab/H/HowellA/faa_files/CasY_MGYP000626431764.fa\
    --output_dir $OUTPUT_DIR

# Visualize result
python3 /dartfs/rc/lab/H/HowellA/visualize_alphafold_results.py $OUTPUT_DIR

# Clean up
rm $SLURM_SUBMIT_DIR/ld.so.cache
conda deactivate