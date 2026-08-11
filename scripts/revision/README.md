# Revision analyses

Scripts for the analyses added during peer review. They are kept separate from `scripts/analysis/` and
`scripts/workflow/` because they follow a different convention: these read parameters from a `config.yaml`
where present, whereas the original scripts carry their parameters inline as `/path/to/...` placeholders. The
original scripts are unchanged — they are the code that produced the submitted results.

## What produces what

| panel | script |
|---|---|
| Fig. 1j — recovered glycan % vs nominal, three read-assignment methods | `em_benchmark/fig1_real_5gly.R` |
| Fig. 1k — recovered vs true fraction, 1-nt divergence | `em_benchmark/fig8_synth_snp.R` |
| Fig. 1l — IP enrichment vs nominal glycan occupancy, inverse-variance-weighted meta-regression | `spike_titration/denovo_doseresponse_metareg.R` |
| Fig. 2e, Fig. 4b — ROC with leave-one-out cross-validation | `classification/counts_analysis_homo_2/revision_roc_var5.r` |
| Fig. 3f,g — spike coverage profile and 5' coverage retention | `coverage/spike_coverage_uncapped.R` |
| Extended Data Fig. 3a — synthetic paralogs | `em_benchmark/fig3_sim_synth.R` |
| Extended Data Fig. 3b — simulated tRNA pair, 1 SNP | `em_benchmark/fig2_sim_trna.R` |
| Extended Data Fig. 3c — molecule-count recovery under PCR duplication | `em_benchmark/fig4_umi_count.R` |
| Extended Data Fig. 3d — modification-induced mismatch | `em_benchmark/fig5_mod_mismatch.R` |
| Extended Data Fig. 3e — mismatch converting the base to the paralog's | `em_benchmark/fig6_mod_decoy_bounded.R` |

`spike_titration/denovo_enrichment_doseresponse.R` computes the per-dose enrichment values that the
meta-regression fits, and `spike_titration/detection_correction.R` is sourced by both.

`fig3_sim_synth.R` plots three divergence facets (1, 4 and 12 nt); Extended Data Fig. 3a shows the 4 and 12 nt
facets, because the 1-nt facet is Fig. 1k.

## Inputs

Sequencing-derived inputs are **not** duplicated here — they come from GEO accession `GSE313898`, which
includes sequencing runs first deposited under `GSE308686` for the preprint version of
this work, or from the upstream steps in `scripts/workflow/` and `scripts/analysis/`:

| script | input | where it comes from |
|---|---|---|
| `coverage/*` | per-position pileups of the two spike constructs, at unlimited depth | deposited in GEO as `spike_pileup_uncapped.csv`. `spike_coverage_uncapped.R` expects one `uncapped_pileups/<library>.d0.txt` per library — the same values in per-library form, with `chr`, `pos`, `ref_base` and a packed `depth,A,C,G,T,N,Skip,Gap,Insert,Delete` field |
| `spike_titration/*` | count matrix, sample table, and the per-group DESeq2 `pI` results | the count matrix is in GEO; the DESeq2 results come from `scripts/analysis/R/counts_analysis/deseq2_pairwise_lrrna_norm.R` |
| `classification/*` | MDS coordinates and distance objects for the cell-line and EV panels | derived from the GEO count matrices by the UMAP/z-delta scripts in `scripts/analysis/R/umap_clustering/` |

The tables under `em_benchmark/` **are** included, because they cannot be obtained from GEO: they are the
output of simulations (synthetic paralog pairs, modification-induced mismatch grids) rather than of any
sequencing run. `real_spike_assignment_n3.tsv` and `denovo_input_glycan_titration_summary.csv` are kept with
them as the observed-data counterparts the same figures compare against.

The two spike-in reference sequences used throughout are in
[`data/transcriptome/spike_references.fa`](../../data/transcriptome/spike_references.fa): `Fluc_spike`, the
unconjugated 89-nt reference, and `glyco_Fluc_spike`, the 51-nt glycan-conjugated construct.

## Running

```sh
cd <directory> && Rscript <script>
```

`spike_titration/config.yaml` uses the same placeholders as the rest of the repository —
`/path/to/metadata.csv`, `/path/to/transcriptome/annotation.csv`, and `count_matrix.csv` relative to the
working directory. Point them at your own files before running. The `em_benchmark/` scripts are
self-contained and run as they are.

`coverage/recompute_5f_uncapped.py` recomputes the 5'-retention metric directly from the pileups as an
independent check of the R implementation. It takes the pileup directory as its argument:

```sh
python3 recompute_5f_uncapped.py uncapped_pileups
```

The depth cap used elsewhere in the pipeline (`mpileup_max_depth`) saturates on these constructs and flattens
the coverage profile, so the 5'-retention metric must be computed from uncapped pileups
(`samtools mpileup -d 0`).

## Dependencies

Beyond the packages listed in the top-level README: `dplyr`, `ggplot2`, `paletteer`, `pROC`, `e1071`,
`vegan`, `ragg`.
