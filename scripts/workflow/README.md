# Preprocessing workflow

FASTQ to the two matrices the analysis half consumes: per-transcript UMI counts
(`count_matrix.csv`, via EM assignment) and per-position pileups (`counts_long.csv`, via
`cpup`). Written for a cluster with a SLURM scheduler, but every step also runs standalone
(see below).

For the general convention on how scripts take their inputs, see "How the scripts are
parameterized" in the [top-level README](../../README.md).

## Layout

```
workflow/
  bash/
    slurm/       # sbatch wrappers: find the inputs, fan out over samples with GNU parallel
    worker/      # the actual work, one sample per invocation; plain bash, no scheduler needed
  python/
    umi_em_dedup.py   # UMI grouping + EM posterior assignment over multi-mapped alignments
```

## Order of processing

| # | worker | invocation | produces |
|---|---|---|---|
| 1 | `cutadapt.sh` | `bash cutadapt.sh input.fastq.gz` | adapter/UMI-trimmed FASTQ |
| 2 | `bowtie2.sh` | `bash bowtie2.sh input.fastq.gz` | BAM, multi-mappings retained |
| 3 | `umi_em_count.sh` | `bash umi_em_count.sh input.bam` | per-transcript UMI counts (`.tsv`) |
| 4 | `umi_em_bam.sh` | `bash umi_em_bam.sh input.bam` | posterior-deduplicated BAM |
| 5 | `cpup.sh` | `bash cpup.sh input.bam output.txt` | per-position base counts |

`umi_em_rarefy.sh` is an optional branch off step 3 for saturation/rarefaction curves.

Each worker reads `WORKDIR` from the environment and writes into subdirectories beneath it,
so a standalone run looks like:

```bash
export WORKDIR=/your/workdir          # must contain raw_fastq/ for step 1
bash worker/cutadapt.sh "$WORKDIR"/raw_fastq/sample.fastq.gz
```

With a scheduler, submit the matching wrapper instead, which globs the input directory and
runs several samples at once:

```bash
sbatch slurm/cutadapt_batch.sh
```

## Resource requests

The `#SBATCH` directives in `slurm/` ask for, per job across a whole cohort:

| step | tasks x cpus | memory | walltime |
| --- | --- | --- | --- |
| `cutadapt_batch.sh` | 6 x 8 | 64 GB | 4 h |
| `bowtie2_batch.sh` | 4 x 8 | 32 GB | 8 h |
| `umi_em_count_batch.sh` | 4 x 2 | 64 GB | 12 h |
| `umi_em_bam_batch.sh` | 4 x 4 | 64 GB | 16 h |
| `cpup_batch.sh` | 4 x 1 | 32 GB | 8 h |
| `umi_em_rarefy_batch.sh` | 4 x 2 | 64 GB | 6 h |

These are the allocations the jobs **request**, not measured peak usage, and each covers
several samples running concurrently. A single sample needs roughly the memory divided by the
task count, so a workstation with 8 or more cores and 32 GB of RAM can process samples
serially through the workers.

## Placeholders to set

Beyond `WORKDIR`:

- **`SCRIPTS_DIR`** resolves to this directory, i.e. `/path/to/rPAL-seq/scripts/workflow` on a plain
  clone. Each call site appends its own subpath from there: the `slurm/` wrappers invoke
  `$SCRIPTS_DIR/bash/worker/<worker>.sh`, and the EM workers invoke
  `$SCRIPTS_DIR/python/umi_em_dedup.py`. One value serves both.
- **`REF`** is set separately in two workers: `bowtie2.sh` wants a bowtie2 index prefix
  (`/path/to/bowtie2/index/prefix`, built from `data/transcriptome/ncrna_*.fa`), while
  `cpup.sh` wants the FASTA itself (`/path/to/transcriptome.fa`).

The EM workers also `source /path/to/venv/bin/activate`; either point that at your Python
environment or delete the line if the interpreter is already on `PATH`.

## Requirements

`cutadapt`, `bowtie2`, `samtools`, [`cpup`](https://github.com/y9c/cpup), `python3` with
`pysam`, and [GNU `parallel`](https://www.gnu.org/software/parallel/), which the `slurm/`
wrappers use to fan out over samples. Versions in the [top-level README](../../README.md).
