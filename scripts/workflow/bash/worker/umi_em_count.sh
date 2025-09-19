#!/bin/bash
# Usage: bash umi_em_count.sh input.bam
set -euo pipefail

# Activate environment
source /path/to/venv/bin/activate

# Set WORKDIR if not externally defined
WORKDIR="${WORKDIR:-/path/to/workdir}"

SCRIPTS_DIR="/path/to/rPAL-seq/scripts"
OUTPUT_DIR="$WORKDIR/umi_em"
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

# Paths (phase-1 artifacts)
QNAME_BAM="$STEP1_DIR/$BAM_BASENAME.qname.bam"                 # name-sorted input for EM
SQLITE_DB="$STEP1_DIR/${BAM_BASENAME}.umi_em.sqlite"          # persisted for phase-2 reuse (marker path)
COUNTS_TSV="$FINAL_DIR/$BAM_BASENAME.counts.tsv"              # counts live in final
EM_LOG="$STATS_DIR/$BAM_BASENAME.umi_em.log"

# UMI/EM parameters
UMI_FROM="qname"            # "qname" or "RX"
QNAME_UMI_SPLIT="_"         # last field after split
UMI_LENGTH="12"             # set to "" to accept all lengths
WEIGHT_MODE="AS"            # "AS" or "NM"
BETA="0.1"
GAMMA="1.0"
MIN_MAPQ="0"

# Speed knobs for shared storage
EM_MAX_ITERS="30"
# EM_TOL is RELATIVE when < 1.0; Python multiplies by usable groups internally.
EM_TOL="1e-3"               # default; may be dynamically overridden below unless DYNAMIC_EM_TOL=0

# Fanout pruning:
# If your aligner is bowtie2 -k 10, setting PERREAD_TOPK=10 (or 12) has no effect on each read.
# Setting it lower (e.g., 8) will drop the lowest-scoring 2 after delta filtering (deterministic by AS/NM+rname).
PERREAD_TOPK="0"            # e.g., 8–12 to cap per-read fanout; 0=disable
DELTA_AS="2"                # only used when WEIGHT_MODE=AS (tight; reduce per-read spread)
# DELTA_NM="1"              # uncomment if using WEIGHT_MODE=NM instead

# Optional: force a clean rebuild from Step 1 to avoid reusing possibly corrupt artifacts
REDO_FROM_STEP1="${REDO_FROM_STEP1:-0}"   # set to 1 to delete prior Step1/EM outputs and rerun
RTOL="${RTOL:-1e-3}"                      # relative tolerance factor (used when DYNAMIC_EM_TOL=1)
DYNAMIC_EM_TOL="${DYNAMIC_EM_TOL:-1}"     # set to 0 to keep EM_TOL as-is

# If requested, purge artifacts so we definitely redo from Step 1
if [[ "$REDO_FROM_STEP1" == "1" ]]; then
    echo ">>> REDO_FROM_STEP1=1 – removing prior Step 1/EM artifacts for a clean rebuild"
    rm -f "$QNAME_BAM" "$SQLITE_DB" "$COUNTS_TSV"
fi

# Skip if counts already exist + DB present (phase-1 complete)
if [[ -f "$COUNTS_TSV" && -f "$SQLITE_DB" ]]; then
    echo ">>> Skipping $BAM_BASENAME – EM already complete (counts + DB present)"
    exit 0
fi

echo ">>> Starting EM (phase-1) for $BAM_BASENAME..."

# Step 0: Save pre-EM flagstat
samtools flagstat "$INPUT_BAM" > "$STATS_DIR/$BAM_BASENAME.before_em.flagstat.txt"

# Step 1: Name-sort the input BAM
if [[ -f "$QNAME_BAM" ]]; then
    echo ">>> Skipping Step 1 for $BAM_BASENAME – name-sorted BAM already exists"
else
    echo ">>> Step 1: name-sorting $INPUT_BAM"
    samtools sort -n -@ "$CPUS_PER_SAMPLE" -m "$MEM_PER_SAMPLE" \
      -T "$STEP1_DIR/${BAM_BASENAME}.name.sort" \
      -o "$QNAME_BAM" "$INPUT_BAM"
fi

# --- Dynamic EM_TOL computation (relative; Python multiplies by usable groups) ---
if [[ "$DYNAMIC_EM_TOL" == "1" ]]; then
    echo ">>> Estimating UMI group count (for logging only) and setting relative EM_TOL (RTOL=$RTOL)..."
    BAM_FOR_COUNT="$QNAME_BAM"
    [[ -f "$BAM_FOR_COUNT" ]] || BAM_FOR_COUNT="$INPUT_BAM"

    # Count distinct UMIs (last field after QNAME_UMI_SPLIT) with optional length filter
    G=$(
      LC_ALL=C samtools view -@ "$CPUS_PER_SAMPLE" "$BAM_FOR_COUNT" \
      | awk -v FS="\t" -v sc="$QNAME_UMI_SPLIT" -v ul="$UMI_LENGTH" '
          {
            n=split($1,a,sc); u=a[n];
            if (ul=="" || length(u)==ul) print u;
          }' \
      | sort -u ${TMPDIR:+--temporary-directory="$TMPDIR"} \
      | wc -l
    )
    [[ -z "$G" || "$G" -lt 1 ]] && G=1
    EM_TOL="$RTOL"  # pass as relative; the Python script converts to absolute
    echo ">>> Dynamic EM_TOL (relative): groups~$G, RTOL=$RTOL -> tol_rel=$EM_TOL" | tee -a "$EM_LOG"
fi
# ---------------------------------------------------------------------------------

# Step 2: EM only (build state + EM -> counts). No hard BAM in phase-1.
if [[ -f "$COUNTS_TSV" && -f "$SQLITE_DB" ]]; then
    echo ">>> Skipping Step 2 for $BAM_BASENAME – EM outputs already exist"
else
    echo ">>> Step 2: running umi_em_dedup.py (phase=em)"
    UMI_LEN_OPT=""
    if [[ -n "${UMI_LENGTH}" ]]; then
        UMI_LEN_OPT="--umi-length ${UMI_LENGTH}"
    fi

    python "$SCRIPTS_DIR/umi_em_dedup.py" "$QNAME_BAM" \
      --phase em \
      --out-counts "$COUNTS_TSV" \
      --umi-from "$UMI_FROM" \
      --qname-umi-split "$QNAME_UMI_SPLIT" \
      $UMI_LEN_OPT \
      --weight-mode "$WEIGHT_MODE" \
      --beta "$BETA" --gamma "$GAMMA" \
      --max-iters "$EM_MAX_ITERS" --tol "$EM_TOL" \
      --tempdb "$SQLITE_DB" \
      --min-mapq "$MIN_MAPQ" \
      ${DELTA_AS:+--perread-delta-as "$DELTA_AS"} \
      ${DELTA_NM:+--perread-delta-nm "$DELTA_NM"} \
      ${PERREAD_TOPK:+--perread-topk "$PERREAD_TOPK"} \
      2> "$EM_LOG"
fi

echo ">>> Finished EM (phase-1) for $BAM_BASENAME at $(date)"

