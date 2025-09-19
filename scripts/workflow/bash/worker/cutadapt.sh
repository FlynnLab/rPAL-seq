#!/bin/bash
# Usage: bash cutadapt.sh input.fastq.gz
set -euo pipefail

# Define WORKDIR fallback if not passed
WORKDIR="${WORKDIR:-/path/to/workdir}"

OUTPUT_DIR="$WORKDIR/cutadapt"
STEP1_DIR="$OUTPUT_DIR/step1"
FINAL_DIR="$OUTPUT_DIR/final"
STATS_DIR="$OUTPUT_DIR/stats"

# Create output directories
mkdir -p "$STEP1_DIR" "$FINAL_DIR" "$STATS_DIR"

# Adapters
POLYA_ADAPTER="A{10}"
POLYG_ADAPTER="G{10}"
ILLUMINA_ADAPTER="AGATCGGAAGAG"

# CPUs per sample
CPUS_PER_SAMPLE=8

# Input file passed as argument
FASTQ_FILE="$1"

# Base name without extensions
base=$(basename "$FASTQ_FILE" .fastq.gz)

# First pass output
step1_file="$STEP1_DIR/${base}.step1.fastq.gz"
step1_log="$STATS_DIR/${base}.step1.log"

# Second pass output
final_file="$FINAL_DIR/${base}.final.fastq.gz"
final_log="$STATS_DIR/${base}.final.log"

# Skip both steps if final output exists
if [[ -s "$final_file" ]]; then
    echo "Skipping $base (already fully trimmed)"
    exit 0
fi

# If step1 doesn't exist, run it
if [[ ! -s "$step1_file" ]]; then
    echo "First Cutadapt pass on $base..."

    cutadapt \
      -j "$CPUS_PER_SAMPLE" \
      -u 12 \
      --rename '{id}_{cut_prefix}' \
      -a "$POLYA_ADAPTER" \
      -a "$POLYG_ADAPTER" \
      -e 0.2 \
      -O 8 \
      -q 20 \
      -m 15 \
      -o "$step1_file" \
      "$FASTQ_FILE" > "$step1_log"

    echo "First Cutadapt pass completed on $base."
else
    echo "Step 1 output found, skipping to step 2 for $base."
fi

# Always run step 2 if final file doesn't exist
echo "Second Cutadapt pass on $base..."

cutadapt \
  -j "$CPUS_PER_SAMPLE" \
  -u 3 \
  -a "$ILLUMINA_ADAPTER" \
  -e 0.1 \
  -O 12 \
  -q 20 \
  -m 15 \
  -o "$final_file" \
  "$step1_file" > "$final_log"

echo "Second Cutadapt pass completed on $base."

