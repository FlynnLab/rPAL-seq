#!/bin/bash
# Usage: bash umi_em_rarefy.sh input.bam
set -euo pipefail

# Activate environment
source /path/to/venv/bin/activate

# Set WORKDIR if not externally defined
WORKDIR="${WORKDIR:-/path/to/workdir}"

SCRIPTS_DIR="/path/to/rPAL-seq/scripts/workflow"
OUTPUT_DIR="$WORKDIR/umi_rarefy"
STATS_DIR="$OUTPUT_DIR/stats"
STEP1_DIR="$OUTPUT_DIR/step1"
FINAL_DIR="$OUTPUT_DIR/final"

# Create output directories
mkdir -p "$STEP1_DIR" "$FINAL_DIR" "$STATS_DIR"

# Resources per sample
CPUS_PER_SAMPLE=2
MEM_PER_SAMPLE=6G

# Input file passed as argument
INPUT_BAM="$1"

# Extract sample name
BAM_BASENAME=$(basename "$INPUT_BAM" .bam)

# Paths (phase-rarefy artifacts)
QNAME_BAM="$STEP1_DIR/$BAM_BASENAME.qname.bam"                  # name-sorted input for rarefaction
RAREFY_TSV="$FINAL_DIR/$BAM_BASENAME.rarefy.tsv"                # main rarefaction output
TAGGED_TSV="$FINAL_DIR/$BAM_BASENAME.rarefy.tagged.tsv"         # optional combined-friendly table
RAREFY_LOG="$STATS_DIR/$BAM_BASENAME.rarefy.log"
SQLITE_DB="$STEP1_DIR/${BAM_BASENAME}.rarefy.sqlite"            # marker file (we'll create it)

# UMI / ingest parameters (kept parallel to EM worker)
UMI_FROM="qname"            # "qname" or "RX"
QNAME_UMI_SPLIT="_"         # last field after split
UMI_LENGTH="12"              # "" to accept all; set "12" for UMI12 runs
WEIGHT_MODE="AS"            # "AS" or "NM"
MIN_MAPQ="0"

# Rarefaction grid and determinism
RAREFY_FRACTIONS="0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0"
RAREFY_SEED="42"            # deterministic sub-sampling by QNAME

# Fanout pruning
PERREAD_TOPK="1"
DELTA_AS="2"
# DELTA_NM="1"              # if WEIGHT_MODE=NM

# Optional knobs
REDO_FROM_STEP1="${REDO_FROM_STEP1:-0}"
WRITE_TAGGED="${WRITE_TAGGED:-0}"
QUIET="${QUIET:-1}"

# If requested, purge artifacts so we definitely redo from Step 1
if [[ "$REDO_FROM_STEP1" == "1" ]]; then
    echo ">>> REDO_FROM_STEP1=1 – removing prior Step 1/rarefy artifacts for a clean rebuild"
    rm -f "$QNAME_BAM" "$SQLITE_DB" "$RAREFY_TSV" "$TAGGED_TSV"
fi

# Skip if rarefaction already exists + marker present
if [[ -f "$RAREFY_TSV" && -f "$SQLITE_DB" ]]; then
    echo ">>> Skipping $BAM_BASENAME – rarefaction already complete (TSV + marker present)"
    exit 0
fi

echo ">>> Starting rarefaction for $BAM_BASENAME..."

# Step 0: Save pre-rarefy flagstat (for context)
samtools flagstat "$INPUT_BAM" > "$STATS_DIR/$BAM_BASENAME.before_rarefy.flagstat.txt"

# Step 1: Name-sort the input BAM
if [[ -f "$QNAME_BAM" ]]; then
    echo ">>> Skipping Step 1 for $BAM_BASENAME – name-sorted BAM already exists"
else
    echo ">>> Step 1: name-sorting $INPUT_BAM"
    samtools sort -n -@ "$CPUS_PER_SAMPLE" -m "$MEM_PER_SAMPLE" \
      -T "$STEP1_DIR/${BAM_BASENAME}.name.sort" \
      -o "$QNAME_BAM" "$INPUT_BAM"
fi

# Step 2: Rarefaction only (stream BAM once; write TSV)
if [[ -f "$RAREFY_TSV" && -f "$SQLITE_DB" ]]; then
    echo ">>> Skipping Step 2 for $BAM_BASENAME – rarefaction outputs already exist"
else
    echo ">>> Step 2: running umi_em_dedup.py (phase=rarefy)"

    UMI_LEN_OPT=()
    [[ -n "$UMI_LENGTH" ]] && UMI_LEN_OPT=(--umi-length "$UMI_LENGTH")

    QUIET_OPT=()
    [[ "$QUIET" == "1" ]] && QUIET_OPT=(--quiet)

    python "$SCRIPTS_DIR/python/umi_em_dedup.py" "$QNAME_BAM" \
      --phase rarefy \
      --rarefy-out "$RAREFY_TSV" \
      --rarefy-fractions "$RAREFY_FRACTIONS" \
      --rarefy-seed "$RAREFY_SEED" \
      --umi-from "$UMI_FROM" \
      --qname-umi-split "$QNAME_UMI_SPLIT" \
      "${UMI_LEN_OPT[@]}" \
      --weight-mode "$WEIGHT_MODE" \
      --min-mapq "$MIN_MAPQ" \
      ${DELTA_AS:+--perread-delta-as "$DELTA_AS"} \
      ${DELTA_NM:+--perread-delta-nm "$DELTA_NM"} \
      ${PERREAD_TOPK:+--perread-topk "$PERREAD_TOPK"} \
      --tempdb "$SQLITE_DB" \
      "${QUIET_OPT[@]}" \
      2> "$RAREFY_LOG"

    # Create a small marker file so the skip check works next time
    echo '{"phase":"rarefy","ok":true}' > "$SQLITE_DB"
fi

# Step 3: Optional tagged TSV for downstream concatenation
if [[ "$WRITE_TAGGED" == "1" ]]; then
    echo ">>> Writing tagged TSV for $BAM_BASENAME"
    {
      echo -e "sample\tfraction\treads\tunique\tsaturation"
      awk -v FS="\t" -v OFS="\t" -v S="$BAM_BASENAME" 'NR>1{print S, $1, $2, $3, $4}' "$RAREFY_TSV"
    } > "$TAGGED_TSV"
fi

echo ">>> Finished rarefaction for $BAM_BASENAME at $(date)"
