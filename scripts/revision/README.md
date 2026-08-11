# Revision analyses

Scripts and inputs for the analyses added during peer review. They are kept separate from
`scripts/analysis/` and `scripts/workflow/` because they follow a different convention: these read
parameters from a `config.yaml` where present, whereas the original scripts carry their parameters inline as
`/path/to/...` placeholders. The original scripts are unchanged — they are the code that produced the
submitted results.

Each directory ships the input tables its scripts read, so every panel below can be regenerated from a clone
without access to the sequencing data. Raw and processed sequencing data are in GEO (see the paper's Data
availability statement).

## What produces what

| panel | script |
|---|---|
| Fig. 1j — recovered glycan % vs nominal, three read-assignment methods | `em_benchmark/fig1_real_5gly.R` |
| Fig. 1k — recovered vs true fraction, 1-nt divergence | `em_benchmark/fig8_synth_snp.R` |
| Fig. 1l — IP enrichment vs nominal glycan occupancy, inverse-variance-weighted meta-regression | `spike_titration/denovo_doseresponse_metareg.R` |
| Fig. 2e, Fig. 4b — ROC with leave-one-out cross-validation | `classification/counts_analysis_homo_2/revision_roc_var5.r` |
| Fig. 3f,g — spike coverage profile and 5' coverage retention | `coverage/spike_coverage_uncapped.R` |
| Extended Data Fig. 3a — synthetic paralogs, 4 and 12 nt divergence | `em_benchmark/fig3_sim_synth.R` |
| Extended Data Fig. 3b — simulated tRNA pair, 1 SNP | `em_benchmark/fig2_sim_trna.R` |
| Extended Data Fig. 3c — molecule-count recovery under PCR duplication | `em_benchmark/fig4_umi_count.R` |
| Extended Data Fig. 3d — modification-induced mismatch | `em_benchmark/fig5_mod_mismatch.R` |
| Extended Data Fig. 3e — mismatch converting the base to the paralog's | `em_benchmark/fig6_mod_decoy_bounded.R` |

`spike_titration/denovo_enrichment_doseresponse.R` computes the per-dose enrichment values that the
meta-regression fits, and `spike_titration/detection_correction.R` is sourced by both.

## Directories

```text
revision/
├── theme_nature.R      # shared ggplot2 theme, sourced by every plotting script
├── em_benchmark/       # read-assignment benchmark; self-contained, ships its result tables
├── spike_titration/    # synthetic sialoglycoRNA standard, dose-response (needs config.yaml)
├── coverage/           # 5' coverage retention, from uncapped pileups
└── classification/     # ROC / LOOCV lineage and EV classification
```

## Running

```sh
cd <directory> && Rscript <script>
```

Two notes on inputs:

- `coverage/` reads `uncapped_pileups/*.d0.txt`, generated with `samtools mpileup -d 0` (unlimited depth)
  restricted to the two spike references. The depth cap used elsewhere in the pipeline
  (`mpileup_max_depth`) saturates on these constructs and flattens the coverage profile, so the 5'-retention
  metric must be computed from uncapped pileups. `recompute_5f_uncapped.py` recomputes the metric directly
  from the same files as an independent check of the R implementation.
- `spike_titration/` reads a count matrix and sample table through `config.yaml`; the paths there are
  relative to the analysis directory as run. Point them at your own files.

## Dependencies

Beyond the packages listed in the top-level README: `dplyr`, `ggplot2`, `paletteer`, `pROC`, `e1071`,
`vegan`, `ragg`.
