# run_all.R ---------------------------------------------------------------
# Runs the full pipeline from the project root. Step 04 takes roughly 20-30 min per window on one core.
steps <- c("R/01_theory_checks.R", "R/02_wahis_attack_rates.R", "R/03_build_blocks.R")
for (s in steps) { message("== ", s); source(s, echo = FALSE) }
system2("Rscript", c("R/04_validation.R", "2005", "2026"))   # primary window
system2("Rscript", c("R/04_validation.R", "2015", "2025"))   # sensitivity window
source("R/05_figures.R")
for (s in c("R/06_secondary_pig_weight.R", "R/07_event_size_exploratory.R", "R/08_within_outbreak_receiver.R")) { message("== ", s); source(s, echo = FALSE) }
sink("output/sessionInfo.txt"); print(sessionInfo()); sink()

# Post-review sensitivity analyses (Table 8)
system2("Rscript", c("R/09_review_sensitivity.R", "w2018_2022", "2018", "2022", "0", "0"))
system2("Rscript", c("R/09_review_sensitivity.R", "w2018_2022_offset", "2018", "2022", "1", "0"))
system2("Rscript", c("R/09_review_sensitivity.R", "primary_offset", "2005", "2026", "1", "0"))
system2("Rscript", c("R/09_review_sensitivity.R", "allcountries", "2005", "2026", "0", "1"))
