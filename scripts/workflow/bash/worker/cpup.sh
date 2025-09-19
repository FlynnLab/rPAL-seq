#!/bin/bash
# Usage: bash cpup.sh input.bam output.txt
set -euo pipefail

# Input arguments
INPUT_BAM=$1
OUTPUT_TXT=$2

# References
REF="/path/to/transcriptome.fa"

echo "Running mpileup on $INPUT_BAM..."

# Run
samtools mpileup \
  -d 10000 \
  -Q 10 \
  --reverse-del \
  -f "$REF" \
  "$INPUT_BAM" | cpup > "$OUTPUT_TXT" 2>&1

  echo "Finished mpileup for $INPUT_BAM at $(date)"
