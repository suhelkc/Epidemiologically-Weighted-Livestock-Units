# 06_secondary_pig_weight.R -------------------------------------------------
# SECONDARY analysis (declared as secondary in the protocol): estimate the pig weight from training
# regions, evaluate on the held-out region. Other weights stay at their frozen prior medians.
# This never alters the primary weights.
source("R/00_setup.R"); set.seed(SEED)
units <- fread(file.path(P$proc, "units.csv")); ob <- fread(file.path(P$proc, "fmd_outbreaks_units.csv"))[!is.na(unit_id)]
units[, n5 := cattle + buffalo + pigs + sheep + goats]
creg <- ob[, .(region = names(which.max(table(region)))), by = iso3]
risk <- merge(units[iso3 %in% ob$iso3 & n5 >= 1], creg, by = "iso3", all.x = TRUE)
risk[, y := as.integer(unit_id %in% ob$unit_id)]
N <- as.matrix(risk[, SPECIES, with = FALSE])
xfun <- function(wpig) { w <- W_MEDIAN[SPECIES]; w["pigs"] <- wpig; as.numeric(scale(log1p(N %*% w))) }
grid <- sort(unique(c(exp(seq(log(0.05), log(30), length.out = 12)), W_MEDIAN["pigs"])))
risk[, x_LU := as.numeric(scale(log1p(LU)))]; risk[, x_UNW := as.numeric(scale(log1p(n5)))]
risk[, x_P := xfun(W_MEDIAN["pigs"])]
fitg <- function(d) suppressWarnings(glmer(y ~ x + (1 | iso3), data = d, family = binomial, nAGQ = 0, control = glmerControl(optimizer = "bobyqa")))
ls_ <- function(y, p) { p <- pmin(pmax(p, 1e-12), 1 - 1e-12); y * log(p) + (1 - y) * log(1 - p) }
res <- rbindlist(lapply(na.omit(unique(risk$region)), function(r) {
  tr <- risk$region != r; te <- risk$region == r
  ll <- sapply(grid, function(g) { d <- data.table(y = risk$y, x = xfun(g), iso3 = risk$iso3)[tr]; as.numeric(logLik(fitg(d))) })
  gs <- grid[which.max(ll)]
  pred <- function(xv) { d <- data.table(y = risk$y, x = xv, iso3 = risk$iso3); predict(fitg(d[tr]), newdata = d[te], re.form = NA, type = "response") }
  data.table(region = r, n = sum(te), pig_weight_hat = gs,
             ls_estimated = mean(ls_(risk$y[te], pred(xfun(gs)))),
             ls_prior_median = mean(ls_(risk$y[te], pred(risk$x_P))),
             ls_LU = mean(ls_(risk$y[te], pred(risk$x_LU))), ls_UNW = mean(ls_(risk$y[te], pred(risk$x_UNW))))
}))
pooled <- res[, lapply(.SD, weighted.mean, w = n), .SDcols = c("ls_estimated", "ls_prior_median", "ls_LU", "ls_UNW")]
fwrite(res, file.path(P$out, "S1_pig_weight_by_region.csv")); fwrite(pooled, file.path(P$out, "S1_pooled_log_scores.csv"))
print(res); print(pooled)
