# 07_event_size_exploratory.R ----------------------------------------------
# EXPLORATORY (post hoc) analysis, formulated AFTER the pre-specified occurrence test (script 04) was null.
# Rationale (Section 2.8): EWLU is proportional to R0, so it should relate to how far an introduced
# epidemic spreads, not to whether a notification occurs (which also depends on introduction pressure
# and reporting). Hypothesis: among notified FMD events, event size increases with the kernel-weighted
# EWLU around the index outbreak. Weights are the frozen prior medians; nothing is tuned.
source("R/00_setup.R"); set.seed(SEED)
B_BOOT <- 2000
units <- fread(file.path(P$proc, "units.csv")); obu <- fread(file.path(P$proc, "fmd_outbreaks_units.csv"))
wah   <- fread(file.path(P$proc, "wahis_subset.csv"))
fmd   <- wah[grepl("Foot and mouth", disease_eng, ignore.case = TRUE) & (as.logical(is_wild) %in% FALSE)]
fmd   <- unique(fmd, by = "Outbreak_id")
fmd[, wt := fifelse(as.logical(Is_cluster) %in% TRUE & !is.na(N_outbreaks_cluster), as.numeric(N_outbreaks_cluster), 1)]
ev <- fmd[, .(size_adj = sum(wt), size_raw = .N, start = min(Outbreak_start_date)), by = epi_event_id]
idx <- fmd[order(Outbreak_start_date, Outbreak_id)][, .SD[1], by = epi_event_id][, .(epi_event_id, Outbreak_id)]
idx <- merge(idx, obu[, .(Outbreak_id, iso3, unit_id)], by = "Outbreak_id")[!is.na(unit_id)]
idx <- merge(idx, units[, .(unit_id, br, bc)], by = "unit_id")
ev  <- merge(ev, idx, by = "epi_event_id")
cat("Events with an index outbreak matched to a unit:", nrow(ev), "in", uniqueN(ev$iso3), "countries\n")

## Block-level indices (all countries pooled), then neighbourhood sums -------------------------
blk <- units[, lapply(.SD, sum), by = .(br, bc), .SDcols = c(SPECIES, "LU")]
blk[, `:=`(EWLU = as.numeric(as.matrix(blk[, SPECIES, with = FALSE]) %*% W_MEDIAN[SPECIES]),
           UNW = cattle + buffalo + pigs + sheep + goats, CAT = cattle)]
nb <- function(r) {
  off <- CJ(dr = -r:r, dc = -r:r)
  j <- off[, .(br = ev$br + dr, bc = ev$bc + dc), by = .(dr, dc)][, ev_row := rep(seq_len(nrow(ev)), times = nrow(off))]
  j <- merge(j, blk[, .(br, bc, EWLU, LU, UNW, CAT)], by = c("br", "bc"))
  j[, lapply(.SD, sum), by = ev_row, .SDcols = c("EWLU", "LU", "UNW", "CAT")][order(ev_row)]
}
run_radius <- function(r) {
  s <- nb(r); d <- copy(ev)[s$ev_row]
  for (p in c("EWLU", "LU", "UNW", "CAT")) set(d, j = paste0("x_", p), value = as.numeric(scale(log1p(s[[p]]))))
  d[, y := log10(size_adj)]
  fitp <- function(tr, te, x) { f <- as.formula(paste("y ~", x, "+ (1 | iso3)"))
    m <- suppressWarnings(lmer(f, data = tr, REML = TRUE)); list(p = predict(m, newdata = te, re.form = NA), m = m) }
  cts <- unique(d$iso3)
  cvres <- rbindlist(lapply(cts, function(cn) {
    tr <- d[iso3 != cn]; te <- d[iso3 == cn]
    out <- data.table(iso3 = cn, y = te$y, p_null = mean(tr$y))
    for (p in c("EWLU", "LU", "UNW", "CAT")) out[[paste0("p_", p)]] <- fitp(tr, te, paste0("x_", p))$p
    out }))
  mse <- sapply(c("null", "EWLU", "LU", "UNW", "CAT"), function(k) mean((cvres$y - cvres[[paste0("p_", k)]])^2))
  boot <- function(k) { pc <- cvres[, .(e = sum((y - get(paste0("p_", k)))^2), n = .N), by = iso3]
    pe <- cvres[, .(e = sum((y - p_EWLU)^2)), by = iso3]$e
    dd <- pe - pc$e; bs <- replicate(B_BOOT, { i <- sample.int(nrow(pc), replace = TRUE); sum(dd[i]) / sum(pc$n[i]) })
    data.table(radius_blocks = r, comparator = k, delta_mse = sum(dd) / sum(pc$n), lo = quantile(bs, .025), hi = quantile(bs, .975)) }
  bd <- rbindlist(lapply(c("null", "LU", "UNW", "CAT"), boot))
  m <- fitp(d, d, "x_EWLU")$m; b <- summary(m)$coefficients["x_EWLU", ]
  list(mse = data.table(radius_blocks = r, as.data.table(as.list(mse))),
       diff = bd,
       slope = data.table(radius_blocks = r, n_events = nrow(d), beta_log10size_per_SD = b[1], se = b[2],
                          lo = b[1] - 1.96 * b[2], hi = b[1] + 1.96 * b[2],
                          spearman_size_EWLU = cor(d$size_adj, d$x_EWLU, method = "spearman"),
                          spearman_size_LU = cor(d$size_adj, d$x_LU, method = "spearman"),
                          spearman_size_UNW = cor(d$size_adj, d$x_UNW, method = "spearman")))
}
R <- lapply(c(0, 1, 2), run_radius)       # 1x1, 3x3 (PRIMARY), 5x5 blocks
fwrite(rbindlist(lapply(R, `[[`, "mse")),   file.path(P$out, "E1_event_size_loco_mse.csv"))
fwrite(rbindlist(lapply(R, `[[`, "diff")),  file.path(P$out, "E1_event_size_delta_mse.csv"))
fwrite(rbindlist(lapply(R, `[[`, "slope")), file.path(P$out, "E1_event_size_slopes.csv"))
for (nm in c("mse", "diff", "slope")) { cat("\n==", nm, "\n"); print(rbindlist(lapply(R, `[[`, nm))[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 4) else z)]) }
