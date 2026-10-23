# 09_review_sensitivity.R -----------------------------------------------------
# POST-REVIEW sensitivity analyses of the occurrence test (added after peer review; not pre-specified).
# Usage: Rscript R/09_review_sensitivity.R <tag> <start_year> <end_year> <offset 0/1> <all_countries 0/1>
#   offset=1      : add log(unit area km2) as an offset (blocks differ in area with latitude)
#   all_countries : include countries with no FMD notification (10-fold country-grouped CV instead of LOCO)
source("R/00_setup.R")
args <- commandArgs(trailingOnly = TRUE)
TAG <- args[1]; WIN <- as.integer(args[2:3]); USE_OFF <- args[4] == "1"; ALLC <- args[5] == "1"
B_BOOT <- 1000; set.seed(SEED)
units <- fread(file.path(P$proc, "units.csv")); ob <- fread(file.path(P$proc, "fmd_outbreaks_units.csv"))[!is.na(unit_id)]
ob <- ob[year >= WIN[1] & year <= WIN[2]]
units[, n5 := cattle + buffalo + pigs + sheep + goats]
Rkm <- 6371.0088; dd <- (5 / 60) * pi / 180
units[, lat := 90 - (br * 6 + 3) * (5 / 60)]
units[, area_km2 := n_cells * Rkm^2 * dd * dd * cos(lat * pi / 180)]
if (ALLC) {                                          # descriptives of excluded countries (written once)
  inc <- unique(ob$iso3)
  desc <- units[n5 >= 1, .(units = .N, cattle = sum(cattle), buffalo = sum(buffalo), pigs = sum(pigs),
                           sheep = sum(sheep), goats = sum(goats)), by = .(included = iso3 %in% inc)]
  desc[, pig_share := pigs / (cattle + buffalo + pigs + sheep + goats)]
  ncty <- units[n5 >= 1, .(n_countries = uniqueN(iso3)), by = .(included = iso3 %in% inc)]
  desc <- merge(desc, ncty, by = "included")
  fwrite(desc, file.path(P$out, "R1_excluded_country_descriptives.csv")); print(desc)
}
risk <- if (ALLC) units[n5 >= 1] else units[iso3 %in% ob$iso3 & n5 >= 1]
risk[, y := as.integer(unit_id %in% ob$unit_id)]
risk[, log_area := log(area_km2) - mean(log(area_km2))]   # centred: only the intercept shifts
cat(sprintf("%s | window %d-%d | offset %s | countries %d | units %d | with outbreak %d\n", TAG, WIN[1], WIN[2], USE_OFF,
            uniqueN(risk$iso3), nrow(risk), sum(risk$y)))
risk[, EWLU := as.matrix(risk[, SPECIES, with = FALSE]) %*% W_MEDIAN[SPECIES]]
risk[, UNW := n5]; risk[, CAT := cattle]
preds <- c("EWLU", "LU", "UNW", "CAT")
for (p in preds) set(risk, j = paste0("x_", p), value = as.numeric(scale(log1p(risk[[p]]))))
ls_ <- function(y, p) { p <- pmin(pmax(p, 1e-12), 1 - 1e-12); y * log(p) + (1 - y) * log(1 - p) }
fit_pred <- function(tr, te, x) {
  f <- as.formula(paste("y ~", x, if (USE_OFF) "+ offset(log_area)" else "", "+ (1 | iso3)"))
  m <- tryCatch(suppressWarnings(glmer(f, data = tr, family = binomial, nAGQ = 0,
                                       control = glmerControl(optimizer = "bobyqa"))),
                error = function(e) tryCatch(suppressWarnings(glmer(f, data = tr, family = binomial, nAGQ = 0,
                                       control = glmerControl(optimizer = "Nelder_Mead"))), error = function(e2) NULL))
  if (!is.null(m)) return(as.numeric(predict(m, newdata = te, re.form = NA, type = "response")))
  # last resort: fixed-effects logistic (no country intercept); recorded in the log
  cat("fallback glm for a fold\n")
  g <- glm(as.formula(paste("y ~", x, if (USE_OFF) "+ offset(log_area)" else "")), data = tr, family = binomial)
  as.numeric(predict(g, newdata = te, type = "response"))
}
ctry <- unique(risk$iso3)
folds <- if (ALLC) split(sample(ctry), rep_len(1:10, length(ctry))) else setNames(as.list(ctry), ctry)
res <- rbindlist(lapply(folds, function(cs) {
  te <- risk[iso3 %in% cs]; tr <- risk[!iso3 %in% cs]
  rbindlist(lapply(preds, function(p) data.table(model = p, iso3 = te$iso3, y = te$y, ls = ls_(te$y, fit_pred(tr, te, paste0("x_", p))))))
}))
pc <- res[, .(L = sum(ls), n = .N), by = .(model, iso3)]
tot <- res[, .(mean_log_score = mean(ls)), by = model][match(preds, model)]
w <- dcast(pc, iso3 ~ model, value.var = "L"); nn <- unique(pc[, .(iso3, n)])[, .(n = n[1]), by = iso3]; w <- merge(w, nn, by = "iso3")
diffs <- rbindlist(lapply(c("LU", "UNW", "CAT"), function(cp) {
  d <- w$EWLU - w[[cp]]
  bs <- replicate(B_BOOT, { i <- sample.int(nrow(w), replace = TRUE); sum(d[i]) / sum(w$n[i]) })
  data.table(comparator = cp, delta_log_score = sum(d) / sum(w$n), lo = quantile(bs, .025), hi = quantile(bs, .975), se = sd(bs))
}))
fwrite(tot,   file.path(P$out, sprintf("R1_%s_log_scores.csv", TAG)))
fwrite(diffs, file.path(P$out, sprintf("R1_%s_diffs.csv", TAG)))
print(tot); print(diffs)
