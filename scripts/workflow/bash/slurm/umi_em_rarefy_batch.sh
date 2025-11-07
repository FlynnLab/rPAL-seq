#!/bin/bash
#SBATCH --job-name=umi_em_rarefy_batch
#SBATCH --output=logs/umi_em_rarefy_batch_%j.out
#SBATCH --error=logs/umi_em_rarefy_batch_%j.err
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --chdir=/path/to/workdir

# Define working directory
WORKDIR="/path/to/workdir"
INPUT_DIR="$WORKDIR/bowtie_alignment/aligned"
LOG_DIR="$WORKDIR/logs"
SCRIPTS_DIR="/path/to/rPAL-seq/scripts"

# Create necessary directories
mkdir -p "$LOG_DIR"

# Find all input BAM files
BAM_LIST=($(find "$INPUT_DIR" -type f -name "*.bam"))

# Export WORKDIR so worker scripts inherit it
export WORKDIR

# Run in parallel
parallel --jobs 4 "bash $SCRIPTS_DIR/umi_em_rarefy.sh {}" ::: "${BAM_LIST[@]}"