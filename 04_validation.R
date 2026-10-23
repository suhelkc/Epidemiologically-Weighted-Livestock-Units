# 04_validation.R ---------------------------------------------------------
# Pre-specified validation of FMD EWLU against WAHIS notifications (Section 5 of the manuscript).
# Weights are fixed at the analytic prior medians (output/frozen_priors.json).
source("R/00_setup.R")
args   <- commandArgs(trailingOnly = TRUE)
WINDOW <- if (length(args) >= 2) as.integer(args[1:2]) else c(2005L, 2026L)
TAG    <- sprintf("%d_%d", WINDOW[1], WINDOW[2])
B_BOOT <- 2000
set.seed(SEED)

units <- fread(file.path(P$proc, "units.csv"))
ob    <- fread(file.path(P$proc, "fmd_outbreaks_units.csv"))
ob    <- ob[!is.na(unit_id) & year >= WINDOW[1] & year <= WINDOW[2]]

## At-risk set: units with susceptible animals in countries notifying >=1 FMD outbreak ----------
units[, n5 := cattle + buffalo + pigs + sheep + goats]
country_region <- ob[, .(region = names(which.max(table(region)))), by = iso3]
risk <- units[iso3 %in% ob$iso3 & n5 >= 1]
risk[, y := as.integer(unit_id %in% ob$unit_id)]
risk <- merge(risk, country_region, by = "iso3", all.x = TRUE)
cat(sprintf("Window %s | countries: %d | units: %d | units with outbreak: %d (%.1f%%)\n",
            TAG, uniqueN(risk$iso3), nrow(risk), sum(risk$y), 100 * mean(risk$y)))

## Predictors (log1p, standardised on the at-risk set) --------------------------
risk[, EWLU := as.matrix(risk[, SPECIES, with = FALSE]) %*% W_MEDIAN[SPECIES]]
risk[, UNW  := n5]
risk[, CAT  := cattle]
preds <- c(EWLU = "EWLU", LU = "LU", UNW = "UNW", CAT = "CAT")
for (p in preds) set(risk, j = paste0("x_", p), value = as.numeric(scale(log1p(risk[[p]]))))

## Helpers -----------------------------------------------------------------------------------
logscore <- function(y, p) { p <- pmin(pmax(p, 1e-12), 1 - 1e-12); y * log(p) + (1 - y) * log(1 - p) }
avg_prec <- function(y, p) { o <- order(p, decreasing = TRUE); y <- y[o]; sum(cumsum(y) / seq_along(y) * y) / sum(y) }
fit_pred <- function(train, test, xname) {
  f <- as.formula(paste("y ~", xname, "+ (1 | iso3)"))
  m <- suppressWarnings(glmer(f, data = train, family = binomial, nAGQ = 0,
                              control = glmerControl(optimizer = "bobyqa")))
  list(p = as.numeric(predict(m, newdata = test, re.form = NA, type = "response")), beta = fixef(m)[2])
}
cv <- function(folds) {                       # folds: named list of iso3 vectors held out together
  res <- rbindlist(lapply(names(folds), function(fn) {
    te <- risk[iso3 %in% folds[[fn]]]; tr <- risk[!iso3 %in% folds[[fn]]]
    rbindlist(lapply(names(preds), function(pn) {
      fp <- fit_pred(tr, te, paste0("x_", preds[[pn]]))
      data.table(fold = fn, model = pn, iso3 = te$iso3, region = te$region, y = te$y, p = fp$p)
    }))
  }))
  res[, ls := logscore(y, p)]
  res
}

## Leave-one-country-out CV ------------------------------------------------------------------
cat("Running leave-one-country-out CV ...\n")
loco <- cv(setNames(as.list(unique(risk$iso3)), unique(risk$iso3)))
## Leave-one-region-out CV (spatial-block CV) -----------------------------------------------
regs <- na.omit(unique(risk$region))
lor  <- cv(setNames(lapply(regs, function(r) unique(risk[region == r, iso3])), regs))

summ <- function(res) {
  pc <- res[, .(L = sum(ls), n = .N, k = sum(y)), by = .(model, iso3, region)]
  tot <- res[, .(log_score = mean(ls),
                 AUC_pooled = as.numeric(auc(roc(y, p, quiet = TRUE))),
                 PR_AUC = avg_prec(y, p),
                 calib_slope = unname(coef(glm(y ~ qlogis(pmin(pmax(p, 1e-12), 1 - 1e-12)), family = binomial))[2])),
             by = model]
  wc <- res[, .(auc = if (sum(y) >= 3 && sum(1 - y) >= 3) as.numeric(auc(roc(y, p, quiet = TRUE))) else NA_real_,
                n = .N), by = .(model, iso3)]
  tot <- merge(tot, wc[, .(AUC_within_country = mean(auc, na.rm = TRUE)), by = model], by = "model")
  list(per_country = pc, total = tot[match(names(preds), model)])
}
S_loco <- summ(loco); S_lor <- summ(lor)
fwrite(S_loco$total, file.path(P$out, sprintf("T5_loco_%s.csv", TAG)))
fwrite(S_lor$total,  file.path(P$out, sprintf("T5_lor_%s.csv",  TAG)))

## Paired country-clustered bootstrap of differences in mean log score --------------------------
boot_diff <- function(pc, comp) {
  w  <- dcast(pc, iso3 + region ~ model, value.var = "L"); nn <- unique(pc[, .(iso3, n)])
  w  <- merge(w, nn, by = "iso3")
  d  <- w$EWLU - w[[comp]]
  idx <- replicate(B_BOOT, sample.int(nrow(w), replace = TRUE))
  bs <- apply(idx, 2, function(i) sum(d[i]) / sum(w$n[i]))
  data.table(comparator = comp, delta_log_score = sum(d) / sum(w$n),
             lo = quantile(bs, 0.025), hi = quantile(bs, 0.975))
}
pc_n <- loco[, .(L = sum(ls), n = .N), by = .(model, iso3, region)]
diffs <- rbindlist(lapply(c("LU", "UNW", "CAT"), boot_diff, pc = pc_n))
fwrite(diffs, file.path(P$out, sprintf("T5_bootstrap_diffs_%s.csv", TAG)))

## Region rule -----------------------------------------------------------------------------------
reg_n <- risk[, .(n_countries = uniqueN(iso3)), by = region][n_countries >= 5]
reg_d <- dcast(pc_n[region %in% reg_n$region, .(L = sum(L), n = sum(n)), by = .(model, region)],
               region ~ model, value.var = "L")
reg_nn <- pc_n[region %in% reg_n$region, .(n = sum(n[model == "EWLU"])), by = region]
reg_d <- merge(reg_d, reg_nn, by = "region")
reg_d[, `:=`(d_vs_LU = (EWLU - LU) / n, d_vs_UNW = (EWLU - UNW) / n, d_vs_CAT = (EWLU - CAT) / n)]
fwrite(reg_d, file.path(P$out, sprintf("T5_region_diffs_%s.csv", TAG)))
maj <- function(v) mean(v >= 0) > 0.5
decision <- list(
  window = TAG,
  ci_excludes_zero_vs_LU  = diffs[comparator == "LU",  lo] > 0,
  ci_excludes_zero_vs_UNW = diffs[comparator == "UNW", lo] > 0,
  majority_regions_not_worse_vs_LU  = maj(reg_d$d_vs_LU),
  majority_regions_not_worse_vs_UNW = maj(reg_d$d_vs_UNW),
  n_regions = nrow(reg_d))
decision$EWLU_supported <- with(decision, ci_excludes_zero_vs_LU && ci_excludes_zero_vs_UNW &&
                                  majority_regions_not_worse_vs_LU && majority_regions_not_worse_vs_UNW)
write_json(decision, file.path(P$out, sprintf("decision_%s.json", TAG)), auto_unbox = TRUE, pretty = TRUE)

## Diagnostics: do scalar EWLU and the (U, V) pair carry the same information? ----------------------
med <- function(s, comp) exp(mean(log(PRIORS[[s]][[comp]])))
wS <- c(cattle = 1, sapply(names(PRIORS), function(s) med(s, "sigma") * sqrt(med(s, "k"))))
wI <- c(cattle = 1, sapply(names(PRIORS), function(s) sqrt(med(s, "k")) * med(s, "Lambda")))
risk[, U := as.matrix(risk[, SPECIES, with = FALSE]) %*% wS[SPECIES]]
risk[, V := as.matrix(risk[, SPECIES, with = FALSE]) %*% wI[SPECIES]]
uv <- data.table(window = TAG,
                 spearman_U_V   = cor(risk$U, risk$V, method = "spearman"),
                 spearman_EWLU_LU = cor(risk$EWLU, risk$LU, method = "spearman"),
                 spearman_EWLU_UNW = cor(risk$EWLU, risk$UNW, method = "spearman"))
fwrite(uv, file.path(P$out, sprintf("T_UV_diagnostic_%s.csv", TAG)))

## In-sample sensitivity to the weight draws (country fixed effects, AIC) -------------------------------
set.seed(SEED)
Wd <- draw_weights(200)
cc <- model.matrix(~ 0 + iso3, risk)
aic_of <- function(x) { m <- glm.fit(cbind(cc, x), risk$y, family = binomial()); m$aic }
aic_lu  <- aic_of(as.numeric(scale(log1p(risk$LU))))
aic_unw <- aic_of(as.numeric(scale(log1p(risk$UNW))))
Nm <- as.matrix(risk[, SPECIES, with = FALSE])
aics <- apply(Wd, 1, function(w) aic_of(as.numeric(scale(log1p(Nm %*% w[SPECIES])))))
sens <- data.table(window = TAG, aic_LU = aic_lu, aic_UNW = aic_unw,
                   aic_median_EWLU = aic_of(risk$x_EWLU),
                   share_draws_better_than_LU  = mean(aics < aic_lu),
                   share_draws_better_than_UNW = mean(aics < aic_unw),
                   delta_aic_vs_LU_median = median(aics - aic_lu),
                   delta_aic_vs_LU_q025 = quantile(aics - aic_lu, 0.025),
                   delta_aic_vs_LU_q975 = quantile(aics - aic_lu, 0.975))
for (s in setdiff(SPECIES, "cattle")) sens[[paste0("spearman_dAIC_", s)]] <- cor(aics, Wd[, s], method = "spearman")
fwrite(sens, file.path(P$out, sprintf("T_weight_draw_sensitivity_%s.csv", TAG)))

cat("\n== LOCO ==\n");  print(S_loco$total[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 4) else z)])
cat("\n== LORO ==\n");  print(S_lor$total[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 4) else z)])
cat("\n== Bootstrap differences (EWLU minus comparator, mean log score) ==\n"); print(diffs[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 5) else z)])
cat("\n== Region differences ==\n"); print(reg_d[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 5) else z)])
cat("\n== Decision ==\n"); print(decision)
cat("\n== U/V diagnostic ==\n"); print(uv)
cat("\n== Weight-draw sensitivity ==\n"); print(sens)
saveRDS(list(loco = loco, lor = lor, risk_n = nrow(risk)), file.path(P$proc, sprintf("cv_%s.rds", TAG)))
