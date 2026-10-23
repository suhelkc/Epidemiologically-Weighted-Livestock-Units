# 00_setup.R -------------------------------------------------------------
# EWLU replication project. Run from the project root (open the .Rproj in RStudio).
# Raw inputs expected in data/raw/ :
#   WAHIS_data.xlsx                         (sheet "infur_20260702", WOAH early-warning extract)
#   GLW4-2020_D-DA_{CTL,BFL,PGS,SHP,GTS,GLEAM3-ALL-LU}.tif   (animals per km2, 5 arcmin)
#   world-administrative-boundaries.shp (+ .shx .dbf .prj .cpg)  (country polygons, field iso3)
suppressPackageStartupMessages({
  library(data.table); library(terra); library(sf); library(readxl)
  library(lme4); library(pROC); library(ggplot2); library(jsonlite)
})
SEED <- 2026L
set.seed(SEED)
terraOptions(progress = 0)

P <- list(
  raw  = "data/raw",
  proc = "data/processed",
  out  = "output"
)
invisible(lapply(c(P$proc, P$out), dir.create, recursive = TRUE, showWarnings = FALSE))

SPECIES <- c("cattle", "buffalo", "pigs", "sheep", "goats")
GLW_FILES <- c(
  cattle  = "GLW4-2020_D-DA_CTL.tif",
  buffalo = "GLW4-2020_D-DA_BFL.tif",
  pigs    = "GLW4-2020_D-DA_PGS.tif",
  sheep   = "GLW4-2020_D-DA_SHP.tif",
  goats   = "GLW4-2020_D-DA_GTS.tif",
  LU      = "GLW4-2020_D-DA_GLEAM3-ALL-LU.tif"   # TLU-type comparator
)

# --- Frozen FMD priors (log-uniform, relative to cattle = 1) -----------------
# Mirrors output/frozen_priors.json; do not edit after the freeze timestamp.
PRIORS <- list(
  buffalo = list(Lambda = c(0.5, 2),  sigma = c(0.5, 1.5), k = c(0.5, 2)),
  pigs    = list(Lambda = c(3, 100),  sigma = c(0.1, 1),   k = c(0.1, 2)),
  sheep   = list(Lambda = c(0.05, 1), sigma = c(0.5, 1.5), k = c(0.25, 2)),
  goats   = list(Lambda = c(0.05, 1), sigma = c(0.3, 1.5), k = c(0.25, 2))
)
# Median of a product of independent log-uniforms = exp(sum of log-midpoints)
prior_median <- function(p) exp(sum(vapply(p, function(r) mean(log(r)), numeric(1))))
W_MEDIAN <- c(cattle = 1, vapply(PRIORS, prior_median, numeric(1)))

# Draw from a log-uniform(lo, hi)
rlunif <- function(n, lo, hi) exp(runif(n, log(lo), log(hi)))
draw_weights <- function(n) {
  out <- sapply(names(PRIORS), function(s) {
    p <- PRIORS[[s]]
    rlunif(n, p$Lambda[1], p$Lambda[2]) * rlunif(n, p$sigma[1], p$sigma[2]) * rlunif(n, p$k[1], p$k[2])
  })
  cbind(cattle = 1, out)
}

# One row per Outbreak_id x Species: WAHIS holds immediate-notification (IN) and follow-up (FUR) reports for the same
# outbreak. Keep the first usable row, preferring IN over FUR, so cases and susceptibles are not double counted.
dedupe_species_rows <- function(dt) {
  dt <- copy(dt)[order(Outbreak_id, Species, IN_FUR != "IN")]
  dt[, .SD[1], by = .(Outbreak_id, Species)]
}
