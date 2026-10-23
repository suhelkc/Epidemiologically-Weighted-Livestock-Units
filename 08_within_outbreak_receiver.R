# 08_within_outbreak_receiver.R ----------------------------------------------
# EXPLORATORY (post hoc). Within a single outbreak, different species share the same exposure and the same
# reporting practice, so the ratio of per-animal attack rates estimates the receiver weight ratio
# w^S = sigma * kappa_in (Proposition 4). Binomial mixed model with an outbreak random intercept.
source("R/00_setup.R"); set.seed(SEED)
wah <- fread(file.path(P$proc, "wahis_subset.csv"))
fmd <- wah[grepl("Foot and mouth", disease_eng, ignore.case = TRUE) & (as.logical(is_wild) %in% FALSE)]
grp <- c(Cattle = "cattle", Buffaloes = "buffalo", Swine = "pigs", Sheep = "sheep", Goats = "goats")
d <- fmd[, sp := grp[Species]][!is.na(sp) & susceptible > 0 & !is.na(cases) & cases <= susceptible]
d <- dedupe_species_rows(d)   # one row per outbreak x species (IN preferred over FUR)
d <- d[, .(cases = sum(cases), susc = sum(susceptible)), by = .(Outbreak_id, epi_event_id, sp)]
multi <- d[, .(k = uniqueN(sp)), by = Outbreak_id][k >= 2, Outbreak_id]
d <- d[Outbreak_id %in% multi & Outbreak_id %in% d[sp == "cattle", Outbreak_id]]   # cattle present as reference
d[, sp := factor(sp, levels = SPECIES)]
cat("Multi-species outbreaks with cattle:", uniqueN(d$Outbreak_id), "| events:", uniqueN(d$epi_event_id), "\n")
print(d[, .(outbreaks = uniqueN(Outbreak_id), animals = sum(susc), cases = sum(cases)), by = sp][order(sp)])
d[, fail := susc - cases]
m <- suppressWarnings(glmer(cbind(cases, fail) ~ sp + (1 | Outbreak_id), data = d, family = binomial,
                            control = glmerControl(optimizer = "bobyqa")))
co <- summary(m)$coefficients; co <- co[grepl("^sp", rownames(co)), , drop = FALSE]
est <- data.table(species = sub("^sp", "", rownames(co)), OR = exp(co[, 1]),
                  lo = exp(co[, 1] - 1.96 * co[, 2]), hi = exp(co[, 1] + 1.96 * co[, 2]), p = co[, 4])
med <- function(s, comp) exp(mean(log(PRIORS[[s]][[comp]])))
est[, prior_wS := sapply(species, function(s) med(s, "sigma") * sqrt(med(s, "k")))]
est[, prior_wS_lo := sapply(species, function(s) PRIORS[[s]]$sigma[1] * sqrt(PRIORS[[s]]$k[1]))]
est[, prior_wS_hi := sapply(species, function(s) PRIORS[[s]]$sigma[2] * sqrt(PRIORS[[s]]$k[2]))]
est[, inside_prior_range := OR >= prior_wS_lo & OR <= prior_wS_hi]
est[, ci_overlaps_prior_range := hi >= prior_wS_lo & lo <= prior_wS_hi]
# Wald intervals ignore extra-binomial variation and event clustering, so they are far too narrow here.
# Report event-clustered bootstrap intervals instead (resample events; refit).
evs <- split(d, d$epi_event_id); B <- 1000
bs <- rbindlist(lapply(seq_len(B), function(b) {
  dd <- rbindlist(evs[sample(names(evs), length(evs), replace = TRUE)], idcol = "rep")
  dd[, Outbreak_id := paste(rep, Outbreak_id)]
  mm <- suppressWarnings(try(glmer(cbind(cases, fail) ~ sp + (1 | Outbreak_id), data = dd, family = binomial,
                                   control = glmerControl(optimizer = "bobyqa")), silent = TRUE))
  if (inherits(mm, "try-error")) return(NULL)
  cc <- fixef(mm); cc <- cc[grepl("^sp", names(cc))]
  data.table(species = sub("^sp", "", names(cc)), OR = exp(cc))
}))
cib <- bs[, .(boot_lo = quantile(OR, .025), boot_hi = quantile(OR, .975), n_ok = .N), by = species]
est <- merge(est, cib, by = "species")
est[, `:=`(boot_overlaps_prior_range = boot_hi >= prior_wS_lo & boot_lo <= prior_wS_hi,
           boot_inside_prior_range = boot_lo >= prior_wS_lo & boot_hi <= prior_wS_hi)]
fwrite(est, file.path(P$out, "E2_within_outbreak_receiver.csv"))
print(est[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])
