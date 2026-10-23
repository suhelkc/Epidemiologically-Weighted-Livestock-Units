# 01_theory_checks.R -----------------------------------------------------
# Numerical verification of Propositions 1-3 and the weight-uncertainty results
source("R/00_setup.R")
set.seed(SEED)

## A. Exactness: spectral radius of K vs c * sum(w N) ---------------------
err_d <- err_f <- numeric(2000)
for (i in seq_len(2000)) {
  S    <- sample(2:6, 1)
  sig  <- rlnorm(S); kin <- rlnorm(S); kout <- rlnorm(S); Lam <- rlnorm(S)
  N    <- rlnorm(S, 8, 1.5); cc <- rlnorm(1, -8, 1)
  w    <- sig * kin * kout * Lam
  K    <- cc * outer(sig * kin * N, kout * Lam)
  r    <- max(Mod(eigen(K, only.values = TRUE)$values))
  err_d[i] <- abs(r - cc * sum(w * N)) / r
  Kf   <- K / sum(N)
  rf   <- max(Mod(eigen(Kf, only.values = TRUE)$values))
  err_f[i] <- abs(rf - cc * sum(w * N) / sum(N)) / rf
}
exactness <- data.table(max_rel_err_density = max(err_d), max_rel_err_frequency = max(err_f))
fwrite(exactness, file.path(P$out, "T_exactness.csv"))

## B. Departure from separable mixing (assortativity alpha) ---------------
sp   <- c("buffalo", "cattle", "pigs", "small_rum")
Lam  <- c(1, 1, 5, 0.3); sig <- c(1, 1, 0.4, 0.8); kin <- kout <- c(1, 1, 0.5, 0.8)
w    <- sig * kin * kout * Lam
tlu  <- c(1, 1, 0.2, 0.1); eq <- rep(1, 4)
M    <- 5000; S <- 4
Nmat <- exp(sapply(log(c(2e4, 4e4, 1e4, 6e4)), function(m) rnorm(M, m, 1)))
ref  <- colMeans(Nmat)
exactR0 <- function(Nrow, al) {
  C <- 1 + al * (S * diag(S) - 1)
  K <- (sig * kin * Nrow) * C * rep(kout * Lam, each = S)   # row s scaled, column j scaled
  max(Mod(eigen(K, only.values = TRUE)$values))
}
grad_w <- function(al, eps = 1e-4) {
  g <- numeric(S)
  for (s in seq_len(S)) {
    a <- b <- ref; a[s] <- a[s] * (1 + eps); b[s] <- b[s] * (1 - eps)
    g[s] <- (exactR0(a, al) - exactR0(b, al)) / (ref[s] * 2 * eps)
  }
  g
}
tabB <- rbindlist(lapply(c(0, 0.1, 0.25, 0.5, 0.75), function(al) {
  R  <- apply(Nmat, 1, exactR0, al = al)
  gw <- grad_w(al)
  data.table(alpha = al,
             EWLU_separable = cor(R, as.numeric(Nmat %*% w),  method = "spearman"),
             EWLU_gradient  = cor(R, as.numeric(Nmat %*% gw), method = "spearman"),
             TLU_type       = cor(R, as.numeric(Nmat %*% tlu), method = "spearman"),
             Equal_weights  = cor(R, as.numeric(Nmat %*% eq),  method = "spearman"))
}))
fwrite(tabB, file.path(P$out, "T2_assortativity_spearman.csv"))

## C. Induced FMD weights from the frozen priors --------------------------
set.seed(SEED)
Wd <- draw_weights(1e5)
tabC <- data.table(species = colnames(Wd),
                   median = apply(Wd, 2, median),
                   q025 = apply(Wd, 2, quantile, 0.025),
                   q975 = apply(Wd, 2, quantile, 0.975),
                   analytic_median = W_MEDIAN[colnames(Wd)])
fwrite(tabC, file.path(P$out, "T3_induced_weights.csv"))

## D. First-order variance shares at global composition -------------------
# Global totals are computed in 03_build_blocks.R (GLW4 counts); fall back to 2020 totals (millions)
Ng <- c(cattle = 1528.1, buffalo = 203.8, pigs = 952.5, sheep = 1237.0, goats = 1131.0)
if (file.exists(file.path(P$proc, "global_totals.csv"))) {
  g <- fread(file.path(P$proc, "global_totals.csv")); Ng <- setNames(g$total_millions, g$species)[SPECIES]
}
E    <- Wd[, SPECIES] %*% Ng[SPECIES]
logE <- log(E[, 1])
V    <- apply(log(Wd[, SPECIES]), 2, var)
pi_s <- (W_MEDIAN[SPECIES] * Ng[SPECIES]) / sum(W_MEDIAN[SPECIES] * Ng[SPECIES])
analytic <- pi_s^2 * V / sum(pi_s^2 * V)
binned <- sapply(SPECIES, function(s) {
  x <- Wd[, s]
  if (s == "cattle") return(0)
  br  <- unique(quantile(x, seq(0, 1, length.out = 51)))
  idx <- cut(x, br, include.lowest = TRUE, labels = FALSE)
  cm  <- tapply(logE, idx, mean)
  var(cm[idx]) / var(logE)
})
tabD <- data.table(species = SPECIES, share_of_EWLU_at_median = as.numeric(pi_s),
                   first_order_analytic = as.numeric(analytic), first_order_binned_MC = as.numeric(binned))
fwrite(tabD, file.path(P$out, "T_variance_shares.csv"))
cat("\nVar(log EWLU) =", round(var(logE), 3), "; geometric SD factor =", round(exp(sd(logE)), 2), "\n")
print(exactness); print(tabB[, lapply(.SD, round, 4)]); print(tabC); print(tabD[, c(list(species = species), lapply(.SD, round, 4)), .SDcols = !"species"])

## E. Spatial extension: individual-level NGM vs unit-level coarse-graining (Section 2.5) ----------
# Individual-level NGM over (unit, species) states, with one infected ANIMAL as the unit of infection:
#   K[(i,s),(j,t)] = c * k_ij * sigma_s kin_s N_is * kout_t Lambda_t
# Proposition 4 (individual level): rho(K) = c * rho(diag(E) %*% k), E_j = sum_s w_s N_js (scalar unit EWLU).
set.seed(SEED)
errE <- errRev <- errU <- numeric(1000)
for (r in seq_len(1000)) {
  n <- sample(3:7, 1); S <- sample(2:5, 1)
  sig <- rlnorm(S); kin <- rlnorm(S); kout <- rlnorm(S); Lam <- rlnorm(S); cc <- rlnorm(1, -2, 0.5)
  N   <- matrix(rlnorm(n * S, 5, 1), n, S)
  xy  <- matrix(runif(2 * n), n); dm <- as.matrix(dist(xy)); kk <- exp(-dm / 0.3)
  w   <- sig * kin * kout * Lam
  a   <- sig * kin; b <- kout * Lam
  idx <- expand.grid(s = seq_len(S), i = seq_len(n))          # state (i,s)
  K <- matrix(0, n * S, n * S)
  for (p in seq_len(n * S)) for (q in seq_len(n * S)) {
    i <- idx$i[p]; s <- idx$s[p]; j <- idx$i[q]; t <- idx$s[q]
    K[p, q] <- cc * kk[i, j] * a[s] * N[i, s] * b[t]
  }
  rhoK <- max(Mod(eigen(K, only.values = TRUE)$values))
  E    <- as.numeric(N %*% w)
  rhoE <- cc * max(Mod(eigen(diag(E) %*% kk, only.values = TRUE)$values))
  rhoRev <- cc * max(Mod(eigen(kk, only.values = TRUE)$values)) * sum(E)   # claimed alternative
  U <- as.numeric(N %*% a); V <- as.numeric(N %*% b)
  rhoUV <- cc * max(Mod(eigen(diag(U) %*% kk %*% diag(V), only.values = TRUE)$values))
  errE[r] <- abs(rhoK - rhoE) / rhoK; errRev[r] <- abs(rhoK - rhoRev) / rhoK; errU[r] <- abs(rhoK - rhoUV) / rhoK
}
spatial <- data.table(model = c("individual-level K vs c*rho(diag(E) k)", "individual-level K vs c*rho(k)*sum(E)", "individual-level K vs unit-level c*rho(diag(U) k diag(V))"),
                      median_rel_err = c(median(errE), median(errRev), median(errU)),
                      max_rel_err = c(max(errE), max(errRev), max(errU)))
fwrite(spatial, file.path(P$out, "T_spatial_ngm_check.csv")); print(spatial)
