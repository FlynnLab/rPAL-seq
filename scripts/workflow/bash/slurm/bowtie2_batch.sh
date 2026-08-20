#!/bin/bash
#SBATCH --job-name=bowtie2_align
#SBATCH --output=logs/bowtie2_align_%A.out
#SBATCH --error=logs/bowtie2_align_%A.err
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --chdir=/path/to/workdir

# Define working directory
WORKDIR="/path/to/workdir"
INPUT_DIR="$WORKDIR/cutadapt/final"
LOG_DIR="$WORKDIR/logs"
SCRIPTS_DIR="/path/to/rPAL-seq/scripts/workflow"

# Create necessary directories
mkdir -p "$LOG_DIR"

# Find all input FASTQ files
FASTQ_LIST=($(find "$INPUT_DIR" -name "*.fastq.gz"))

# Export WORKDIR so worker scripts inherit it
export WORKDIR

# Run in parallel
parallel --jobs 4 "bash $SCRIPTS_DIR/bash/worker/bowtie2.sh {}" ::: "${FASTQ_LIST[@]}"
