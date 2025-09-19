#!/bin/bash
#SBATCH --job-name=cutadapt_trim
#SBATCH --output=logs/cutadapt_trim_%A.out
#SBATCH --error=logs/cutadapt_trim_%A.err
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --chdir=/path/to/workdir

# Define working directory
WORKDIR="/path/to/workdir"
INPUT_DIR="$WORKDIR/raw_fastq"
LOG_DIR="$WORKDIR/logs"
SCRIPTS_DIR="/path/to/rPAL-seq/scripts"

# Create output directories
mkdir -p "$LOG_DIR"

# Find all input FASTQ files
FASTQ_LIST=($(find "$INPUT_DIR" -name "*.fastq.gz"))

# Export WORKDIR so worker script inherits it
export WORKDIR

# Run in parallel
parallel --jobs 6 "bash $SCRIPTS_DIR/cutadapt.sh {}" ::: "${FASTQ_LIST[@]}"
