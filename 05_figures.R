# 05_figures.R ------------------------------------------------------------
# Manuscript figures (vector PDF + 300 dpi PNG). Colour-blind-safe palette (Okabe-Ito).
source("R/00_setup.R")
pal <- c(cattle = "#0072B2", buffalo = "#009E73", pigs = "#D55E00", sheep = "#E69F00", goats = "#CC79A7")
th  <- theme_classic(base_size = 10) + theme(legend.position = "bottom")
save_fig <- function(p, name, w = 7, h = 4) {
  ggsave(file.path(P$out, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(P$out, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

## Fig 1. Induced FMD weights (log scale) ---------------------------------
t3 <- fread(file.path(P$out, "T3_induced_weights.csv"))[species != "cattle"]
f1 <- ggplot(t3, aes(x = median, y = reorder(species, median), colour = species)) +
  geom_pointrange(aes(xmin = q025, xmax = q975), linewidth = 0.7) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  scale_x_log10() + scale_colour_manual(values = pal, guide = "none") +
  labs(x = "Weight relative to cattle (log scale; median and 95% interval)", y = NULL) + th
save_fig(f1, "Fig1_induced_weights", 6, 3)

## Fig 2. Departure from separable mixing -----------------------------------
t2 <- melt(fread(file.path(P$out, "T2_assortativity_spearman.csv")), id.vars = "alpha",
           variable.name = "index", value.name = "rho")
f2 <- ggplot(t2, aes(alpha, rho, colour = index, linetype = index)) + geom_line(linewidth = 0.8) + geom_point() +
  scale_colour_manual(values = c("#0072B2", "#56B4E9", "#D55E00", "grey40")) +
  labs(x = "Within-species assortativity (alpha)", y = "Spearman rho with exact R0", colour = NULL, linetype = NULL) + th
save_fig(f2, "Fig2_assortativity", 6, 4)

## Fig 3. First-order variance share ----------------------------------------
t6 <- melt(fread(file.path(P$out, "T_variance_shares.csv"))[species != "cattle"],
           id.vars = "species", measure.vars = c("first_order_analytic", "first_order_binned_MC"))
f3 <- ggplot(t6, aes(species, value, fill = variable)) + geom_col(position = "dodge") +
  scale_fill_manual(values = c("#0072B2", "#E69F00"), labels = c("Analytic (pi^2 V)", "Binned Monte Carlo")) +
  labs(x = NULL, y = "Share of Var(ln EWLU)", fill = NULL) + th
save_fig(f3, "Fig3_variance_shares", 6, 4)

## Fig 4. Attack rates -------------------------------------------------------
t4 <- fread(file.path(P$out, "T4_attack_rates.csv"))
f4 <- ggplot(t4, aes(100 * ar, reorder(sp, ar))) + geom_pointrange(aes(xmin = 100 * lo, xmax = 100 * hi)) +
  labs(x = "Pooled within-outbreak attack rate (%, event-clustered bootstrap 95% CI)", y = NULL) + th
save_fig(f4, "Fig4_attack_rates", 6, 3.5)

## Fig 5. Held-out log-score differences by region ---------------------------
rf <- list.files(P$out, "^T5_region_diffs_2005_2026.csv$", full.names = TRUE)
if (length(rf)) {
  r <- melt(fread(rf), id.vars = "region", measure.vars = c("d_vs_LU", "d_vs_UNW", "d_vs_CAT"))
  f5 <- ggplot(r, aes(value, region, colour = variable)) + geom_point(size = 2.5) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_colour_manual(values = c("#D55E00", "#0072B2", "#009E73"), labels = c("vs LU", "vs unweighted", "vs cattle")) +
    labs(x = "Mean held-out log score, EWLU minus comparator", y = NULL, colour = NULL) + th
  save_fig(f5, "Fig5_region_log_score_diff", 6.5, 3.5)
}
