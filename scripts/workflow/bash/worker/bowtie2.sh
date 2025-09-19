#!/bin/bash
# Usage: bash bowtie2.sh input.fastq.gz
set -euo pipefail

# Define WORKDIR fallback if not passed
WORKDIR="${WORKDIR:-/path/to/workdir}"

OUTPUT_DIR="$WORKDIR/bowtie_alignment"
ALIGNED_DIR="$OUTPUT_DIR/aligned"
UNALIGNED_DIR="$OUTPUT_DIR/unaligned"
STATS_DIR="$OUTPUT_DIR/stats"

# Reference transcriptome index
REF="/path/to/bowtie2/index/prefix"

# Create necessary directories
mkdir -p "$ALIGNED_DIR" "$UNALIGNED_DIR" "$STATS_DIR"

# Resources per sample
CPUS_PER_SAMPLE=8

# Input file passed as argument
FASTQ_FILE="$1"

# Extract sample name
sample_name=$(basename "$FASTQ_FILE" .fastq.gz)
bam_file="$ALIGNED_DIR/${sample_name}_aligned.unsorted.bam"

# Skip if already aligned
if [[ -f "$bam_file" ]]; then
    echo "Skipping $sample_name (already aligned)"
    exit 0
fi

echo "Aligning and preparing $sample_name..."

# Run bowtie2 -> write compressed UNSORTED BAM
bowtie2 \
    --very-sensitive-local \
    -k 10 \
    --norc \
    -p "$CPUS_PER_SAMPLE" \
    -x "$REF" \
    -U "$FASTQ_FILE" \
    --un-gz "$UNALIGNED_DIR/${sample_name}_unaligned.fastq.gz" \
    --no-unal \
    2> "$STATS_DIR/${sample_name}_stats.txt" | \
samtools view -b -@ "$CPUS_PER_SAMPLE" -o "$bam_file" -

echo "Finished $sample_name at $(date)"

