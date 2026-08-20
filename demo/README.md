# Worked example: HeLa glycoRNA calls from public processed data

This reproduces the transcript-level glycoRNA calls for HeLa using only a 0.2 MB download. It
exercises the analysis half of the pipeline (`scripts/analysis/`) and needs no cluster, no alignment
index and no FASTQ. Total run time is about 10 to 12 seconds on a laptop.

It is also the shortest way to confirm your installation works before pointing the scripts at your own
data. For the general convention on how scripts take their inputs, see "How the scripts are
parameterized" in the [top-level README](../README.md).

## Contents

| file | purpose |
| --- | --- |
| `metadata_demo.csv` | the four HeLa biological replicates, in the shape the scripts expect |

The count matrix itself is not vendored here; it is downloaded from GEO in the first step below.

## 1. Get the demo data

The processed count matrix is a supplementary file of GEO accession
[GSE308686](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308686) and is 0.2 MB compressed:

Run these from the **repository root**, so that `../data` and `../scripts` resolve:

```bash
mkdir -p demo_run && cd demo_run
curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE308nnn/GSE308686/suppl/GSE308686_count_matrix_homo.csv.gz
gunzip GSE308686_count_matrix_homo.csv.gz
cp ../demo/metadata_demo.csv .
```

`GSE308686_count_matrix_homo.csv` holds 2,317 transcripts across 196 libraries.
`demo/metadata_demo.csv` selects the four HeLa biological replicates (`VH1` to `VH4`); the scripts
append the `i`/`p`/`h` suffixes to reach the Input, IP and Control library of each replicate, so the
example runs on 12 of those 196 columns.

## 2. Run it

Both scripts read their inputs from paths declared at the top of the file (see Instructions for use).
Point them at the demo files, then run from inside `demo_run`:

```bash
sed -e 's|/path/to/count_matrix.csv|GSE308686_count_matrix_homo.csv|' \
    -e 's|/path/to/metadata.csv|metadata_demo.csv|' \
    -e 's|/path/to/transcriptome/annotation.csv|../data/annotation/transcript2family_homo.csv|' \
    ../scripts/analysis/R/counts_analysis/deseq2_pairwise_lrrna_norm.R > step1.R

sed -e 's|/path/to/metadata.csv|metadata_demo.csv|' \
    ../scripts/analysis/R/counts_analysis/tp_wald_sizeadjust.R > step2.R

Rscript step1.R    # DESeq2 IP-vs-Input and Control-vs-Input, long-rRNA normalized
Rscript step2.R    # post hoc Wald test against the mock-release control
```

## 3. Expected output

```
deseq2_results_pI_HeLa.csv    IP vs Input contrast
deseq2_results_cI_HeLa.csv    Control vs Input contrast
wald_all_HeLa.csv             every transcript with both contrasts and the delta statistic
wald_dir_HeLa.csv             the directional subset (log2FC IP/Input > 0.5)
wald_TP_HeLa.csv              216 true-positive glycoRNA calls
```

`wald_TP_HeLa.csv` carries 216 rows and 21 columns. The strongest calls are small nuclear and
nucleolar RNAs and a pre-miRNA, for example `snRNA_LOC124906135_29211` at log2FC 2.96, matching the
biotype composition of the example volcano in `doc/volcano_tp_example.png`.

## 4. Expected run time

About 10 to 12 seconds total: 7 to 9 seconds for `step1.R` and about 3 seconds for `step2.R`, over two
runs on an 8-core Apple M3 with 8 GB of RAM under macOS 14.2. That measurement used R 4.6.0 with DESeq2 1.52.0
and limma 3.68.3 rather than the pinned versions above; the timing is not sensitive to the difference,
but exact numeric output can shift slightly between DESeq2 releases.

To exercise the per-base variant arm as well, `GSE308686_counts_long.csv.gz` from the same accession
is the corresponding input for `scripts/analysis/R/variant_analysis/`. It is 117 MB compressed and so
is not part of the quick demo.
