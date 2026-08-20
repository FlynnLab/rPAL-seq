#!/bin/bash
#SBATCH --job-name=cpup
#SBATCH --output=logs/cpup_%A.out
#SBATCH --error=logs/cpup_%A.err
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --chdir=/path/to/workdir

# Define working directory
WORKDIR="/path/to/workdir"
INPUT_DIR="$WORKDIR/umi_em/bam"
OUTPUT_DIR="$WORKDIR/pileup_results"
LOG_DIR="$WORKDIR/logs"
SCRIPTS_DIR="/path/to/rPAL-seq/scripts/workflow"

# Create necessary directories
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Find all input BAM files
BAM_LIST=($(find "$INPUT_DIR" -name "*.emhard.coord.bam"))

# Run in parallel
parallel --jobs 4 "bash $SCRIPTS_DIR/bash/worker/cpup.sh {} $OUTPUT_DIR/{/.}.cpup.txt" ::: "${BAM_LIST[@]}"
