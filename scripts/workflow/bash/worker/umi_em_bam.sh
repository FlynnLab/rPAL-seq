#!/bin/bash
# Usage: bash umi_em_bam.sh input.bam
set -euo pipefail

# Activate environment
source /path/to/venv/bin/activate

# Set WORKDIR if not externally defined
WORKDIR="${WORKDIR:-/path/to/workdir}"

SCRIPTS_DIR="/path/to/rPAL-seq/scripts/workflow"
OUTPUT_DIR="$WORKDIR/umi_em"
STATS_DIR="$OUTPUT_DIR/stats"
STEP1_DIR="$OUTPUT_DIR/step1"
FINAL_DIR="$OUTPUT_DIR/final"
BAMS_DIR="$OUTPUT_DIR/bam"

# Create output directories
mkdir -p "$STEP1_DIR" "$FINAL_DIR" "$STATS_DIR" "$BAMS_DIR"

# Resources per sample
CPUS_PER_SAMPLE=4
MEM_PER_SAMPLE=3G

# Phase-2 knobs (override via env if desired)
ASSIGN_WORKERS=3             # parallel shard DBM build
IO_THREADS=2                 # pysam/htslib IO threads
STREAM_SORT="1"              # stream hard BAM directly into samtools sort

# Input file passed as argument
INPUT_BAM="$1"
BAM_BASENAME=$(basename "$INPUT_BAM" .bam)

# Reuse artifacts from phase-1
QNAME_BAM="$STEP1_DIR/$BAM_BASENAME.qname.bam"           # must exist from phase-1
SQLITE_DB="$STEP1_DIR/${BAM_BASENAME}.umi_em.sqlite"     # marker path from phase-1
COUNTS_TSV="$FINAL_DIR/$BAM_BASENAME.counts.tsv"         # must exist from phase-1

# Phase-2 outputs
EM_QNAME_BAM="$STEP1_DIR/$BAM_BASENAME.emhard.qname.bam" # name-sorted hard pick (optional if STREAM_SORT=1)
FINAL_BAM="$BAMS_DIR/$BAM_BASENAME.emhard.coord.bam"    # coordinate-sorted final
FINAL_INDEX="$FINAL_BAM.bai"
DEDUP_LOG="$STATS_DIR/$BAM_BASENAME.umi_dedup.log"

# Hard-pick parameters (do not affect EM)
ASSIGNED_MAPQ="${ASSIGNED_MAPQ:-255}"
MIN_MAPQ="0"                 # will be overridden by meta.json if available
KEEP_UNASSIGNED_BAM=""       # set to a path for debugging, else leave empty
OUT_HARD_TSV=""              # optional TSV of hard assignments; set path or leave empty
NO_TAGS="1"                  # set to 1 to skip HP/HC/HB tags

# ---- Load EM meta so phase-2 uses exactly the same settings as phase-1 ----
# Defaults (used only if meta.json is missing)
UMI_FROM="qname"
QNAME_UMI_SPLIT="_"
UMI_LENGTH=""                # empty -> accept all lengths
WEIGHT_MODE="AS"
BETA="0.1"
GAMMA="1.0"

META_JSON="$STEP1_DIR/${BAM_BASENAME}.umi_em.sqlite.state/meta.json"
if [[ -f "$META_JSON" ]]; then
    echo ">>> Loading EM meta: $META_JSON"
    # Emit shell assignments safely quoted and eval them
    eval "$(
      python - "$META_JSON" <<'PY'
import json, sys, shlex
p = sys.argv[1]
with open(p, 'r', encoding='utf-8') as fh:
    m = json.load(fh)

def setvar(k, v):
    print(f"{k}={shlex.quote(str(v))}")

setvar("UMI_FROM",        m.get("umi_from", "qname"))
setvar("QNAME_UMI_SPLIT", m.get("qname_split", "_"))
setvar("UMI_LENGTH",      m.get("umi_length", ""))          # "" means 'no filter'
setvar("WEIGHT_MODE",     m.get("weight_mode", "AS"))
setvar("BETA",            m.get("beta", "0.1"))
setvar("GAMMA",           m.get("gamma", "1.0"))
setvar("MIN_MAPQ",        m.get("min_mapq", "0"))
PY
    )"
else
    echo ">>> WARNING: meta.json not found; using default UMI/weight parameters for phase-2"
fi
# ---------------------------------------------------------------------------

# Sanity: required inputs from phase-1
[[ -f "$QNAME_BAM" ]] || { echo "ERROR: Missing name-sorted BAM: $QNAME_BAM" >&2; exit 2; }
[[ -f "$SQLITE_DB" ]] || { echo "ERROR: Missing SQLite DB marker: $SQLITE_DB" >&2; exit 2; }
[[ -f "$COUNTS_TSV" ]] || { echo "ERROR: Missing counts TSV:    $COUNTS_TSV" >&2; exit 2; }

# Skip if final outputs already present
if [[ -f "$FINAL_BAM" && -f "$FINAL_INDEX" ]]; then
    echo ">>> Skipping $BAM_BASENAME – hard-pick BAM already complete"
    exit 0
fi

echo ">>> Starting hard-pick (phase-2) for $BAM_BASENAME..."

# Step A: Hard-pick on name-sorted BAM using existing state + counts
if [[ -f "$EM_QNAME_BAM" && "$STREAM_SORT" != "1" ]]; then
    echo ">>> Skipping Step A – hard-assigned qname BAM already exists"
else
    echo ">>> Step A: running umi_em_dedup.py (phase=dedup)"
    OUT_HARD_OPT=()
    [[ -n "$OUT_HARD_TSV" ]] && OUT_HARD_OPT=(--out-hard "$OUT_HARD_TSV")

    KEEP_UNASSIGNED_OPT=()
    [[ -n "$KEEP_UNASSIGNED_BAM" ]] && KEEP_UNASSIGNED_OPT=(--keep-unassigned-bam "$KEEP_UNASSIGNED_BAM")

    UMI_LEN_OPT=()
    [[ -n "$UMI_LENGTH" ]] && UMI_LEN_OPT=(--umi-length "$UMI_LENGTH")

    NO_TAGS_OPT=()
    [[ "$NO_TAGS" == "1" ]] && NO_TAGS_OPT=(--no-tags)

    if [[ "$STREAM_SORT" == "1" ]]; then
      echo ">>> Streaming hard-pick into samtools sort (no intermediate qname BAM)"
      set -o pipefail
      python "$SCRIPTS_DIR/python/umi_em_dedup.py" "$QNAME_BAM" \
        --phase dedup \
        --out-counts "$COUNTS_TSV" \
        --out-bam - \
        --tempdb "$SQLITE_DB" \
        --assigned-mapq "$ASSIGNED_MAPQ" \
        --min-mapq "$MIN_MAPQ" \
        --umi-from "$UMI_FROM" \
        --qname-umi-split "$QNAME_UMI_SPLIT" \
        "${UMI_LEN_OPT[@]}" \
        --weight-mode "$WEIGHT_MODE" \
        --beta "$BETA" --gamma "$GAMMA" \
        --assign-workers "$ASSIGN_WORKERS" \
        --io-threads "$IO_THREADS" \
        "${NO_TAGS_OPT[@]}" \
        "${OUT_HARD_OPT[@]}" \
        "${KEEP_UNASSIGNED_OPT[@]}" \
        2> "$DEDUP_LOG" \
      | samtools sort -@ "$CPUS_PER_SAMPLE" -m "$MEM_PER_SAMPLE" \
          -T "$BAMS_DIR/${BAM_BASENAME}.coord.sort" \
          -o "$FINAL_BAM"
      samtools index "$FINAL_BAM"
      # mark that EM_QNAME_BAM is 'produced' via stream to keep idempotence logic happy
      : > "$EM_QNAME_BAM"
    else
      python "$SCRIPTS_DIR/python/umi_em_dedup.py" "$QNAME_BAM" \
        --phase dedup \
        --out-counts "$COUNTS_TSV" \
        --out-bam "$EM_QNAME_BAM" \
        --tempdb "$SQLITE_DB" \
        --assigned-mapq "$ASSIGNED_MAPQ" \
        --min-mapq "$MIN_MAPQ" \
        --umi-from "$UMI_FROM" \
        --qname-umi-split "$QNAME_UMI_SPLIT" \
        "${UMI_LEN_OPT[@]}" \
        --weight-mode "$WEIGHT_MODE" \
        --beta "$BETA" --gamma "$GAMMA" \
        --assign-workers "$ASSIGN_WORKERS" \
        --io-threads "$IO_THREADS" \
        "${NO_TAGS_OPT[@]}" \
        "${OUT_HARD_OPT[@]}" \
        "${KEEP_UNASSIGNED_OPT[@]}" \
        2> "$DEDUP_LOG"
    fi
fi

# Step B: Coordinate-sort & index the collapsed BAM (skipped if we already streamed)
if [[ "$STREAM_SORT" == "1" ]]; then
  : # already sorted & indexed above
elif [[ -f "$FINAL_BAM" && -f "$FINAL_INDEX" ]]; then
    echo ">>> Skipping Step B – final BAM already exists"
else
    echo ">>> Step B: coordinate-sorting EM-collapsed BAM"
    samtools sort -@ "$CPUS_PER_SAMPLE" -m "$MEM_PER_SAMPLE" \
      -T "$BAMS_DIR/${BAM_BASENAME}.coord.sort" \
      -o "$FINAL_BAM" "$EM_QNAME_BAM"
    samtools index "$FINAL_BAM"
fi

# Step C: Post-EM flagstat
samtools flagstat "$FINAL_BAM" > "$STATS_DIR/$BAM_BASENAME.after_em.flagstat.txt"

echo ">>> Finished hard-pick (phase-2) for $BAM_BASENAME at $(date)"

