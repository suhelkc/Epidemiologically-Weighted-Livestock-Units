# EWLU (FMD) replication project

Open the folder as an RStudio Project; all paths are relative to the project root.

## Inputs (data/raw/)
| File | Source |
|---|---|
| WAHIS_data.xlsx (sheet `infur_20260702`) | WOAH WAHIS early-warning extract. Required citation text is in the workbook's "Copyright and disclaimer" sheet. |
| GLW4-2020_D-DA_{CTL,BFL,PGS,SHP,GTS}.tif | FAO Gridded Livestock of the World v4, 2020, dasymetric density (animals/km2), 5 arcmin |
| GLW4-2020_D-DA_GLEAM3-ALL-LU.tif | GLEAM3 all-species livestock units (TLU-type comparator) |
| world-administrative-boundaries.* | Country polygons with field `iso3` |

## Order of scripts
1. `R/01_theory_checks.R` Propositions 1-3, assortativity experiment, induced weights, variance shares (Tables 2, 3, S-variance)
2. `R/02_wahis_attack_rates.R` pooled attack rates with event-clustered bootstrap (Table 4)
3. `R/03_build_blocks.R` rasters -> counts -> 0.5-degree block x country units; outbreaks -> units
4. `R/04_validation.R <start> <end>` pre-specified validation (Table 5, decision rule, diagnostics)
5. `R/05_figures.R` figures
6. `R/06_secondary_pig_weight.R` secondary: pig weight estimated in training regions, tested in the held-out region
7. `R/07_event_size_exploratory.R` EXPLORATORY (post hoc): event size vs neighbourhood indices (H2)
8. `R/08_within_outbreak_receiver.R` EXPLORATORY (post hoc): within-outbreak species attack-odds ratios (H3)

9. `R/09_review_sensitivity.R <tag> <start> <end> <offset 0/1> <all_countries 0/1>` POST-REVIEW sensitivity analyses of the occurrence test (narrow window, area offset, all countries); run once per row of Table 8

Scripts 06-09 and the Section 5.3-5.4 analyses were formulated after the pre-specified test (script 04) was null; see the manuscript for the full record.

`R/00_setup.R` holds the frozen priors (see `output/frozen_priors.json`, SHA-256 recorded at freeze).

## Notes
- Reading the WAHIS workbook with readxl needs several GB of RAM. Script 02 exports the 14 columns it uses to
  `data/processed/wahis_subset.csv` on first run (or reads that CSV if present).
- All random processes use seed 2026.
- Required packages: data.table, terra, sf, readxl, lme4, pROC, ggplot2, jsonlite.

- Script 01 (section E) numerically checks Proposition 4 (individual-level vs unit-level spatial matrices).
- Background runs of script 04/09 take about 10 minutes each (leave-one-country-out); the 2015-2025 weight-draw sensitivity output did not complete in the authors' environment.
