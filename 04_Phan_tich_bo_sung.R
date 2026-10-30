# =============== 04. SUPPLEMENTARY ANALYSES FOR THE MANUSCRIPT ===============
# Omega reliability, ordinal (WLSMV) re-estimation, ESEM check of the three facets,
# latent mean differences under scalar invariance, and a MIMIC model testing whether
# group differences in AI career anxiety remain after controlling threat appraisal.
# Writes ket_qua/so_lieu_bai_bao.json so the manuscript numbers can be checked.
# Run after 00_CONFIG.R (CHAY_TAT_CA.R sources it in order).

options(repos = c(CRAN = "https://cloud.r-project.org"))
for (p in c("lavaan", "semTools", "jsonlite")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressMessages({ library(lavaan); library(semTools); library(jsonlite) })

dat <- read.csv(DATA_FILE, fileEncoding = "UTF-8-BOM")
ANX <- paste0("ANX", 1:9); THR <- paste0("THR", 1:3); ASE <- paste0("ASE", 1:3); CSE <- paste0("CSE", 1:3)
ALL <- c(ANX, THR, ASE, CSE)
dat$grp_disc   <- ifelse(dat$khoi_nganh == "Công nghệ thông tin", "IT", "Non-IT")
dat$grp_year   <- ifelse(dat$nam_hoc %in% c("Năm 1", "Năm 2"), "Junior", "Senior")
dat$grp_gender <- ifelse(dat$gioi_tinh == "Nam", "Male", ifelse(dat$gioi_tinh == "Nữ", "Female", NA))

R <- list(N = nrow(dat))
out <- file.path(OUT_TXT, "Ket_qua_bo_sung.txt")
sink(out, split = TRUE)
cat("SUPPLEMENTARY ANALYSES | N =", nrow(dat), "|", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")

m_uni <- paste("ANX =~", paste(ANX, collapse = " + "))
m_3f  <- "OBS =~ ANX1 + ANX2 + ANX3\nOPP =~ ANX4 + ANX5 + ANX6\nCMP =~ ANX7 + ANX8 + ANX9"
meas  <- paste(m_uni, paste("THR =~", paste(THR, collapse = " + ")),
               paste("ASE =~", paste(ASE, collapse = " + ")), paste("CSE =~", paste(CSE, collapse = " + ")), sep = "\n")
one   <- paste("F =~", paste(c(ANX, THR), collapse = " + "))
two   <- paste(m_uni, paste("THR =~", paste(THR, collapse = " + ")), sep = "\n")

## ---- A. Omega (from the MLR measurement model) ----
f_meas <- cfa(meas, dat, estimator = "MLR")
om <- sapply(compRelSEM(f_meas), as.numeric)  # McDonald's omega per factor (semTools >= 0.5-9 returns a list)
f_uni  <- cfa(m_uni, dat, estimator = "MLR")
om_uni <- as.numeric(unlist(compRelSEM(f_uni)))
cat("\n== A. Omega ==\n"); print(round(om, 3)); cat("ANX (one-factor model):", round(om_uni, 3), "\n")
R$omega <- as.list(round(unlist(om), 3)); R$omega_ANX_onefactor <- round(unname(unlist(om_uni)), 3)

## ---- B. Ordinal re-estimation (WLSMV, items as ordered) ----
ix <- c("chisq.scaled", "df.scaled", "pvalue.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr")
wfit <- function(model, vars) cfa(model, dat, ordered = vars, estimator = "WLSMV")
w_uni  <- wfit(m_uni, ANX)
w_3f   <- wfit(m_3f, ANX)
w_meas <- wfit(meas, ALL)
w_one  <- wfit(one, c(ANX, THR)); w_two <- wfit(two, c(ANX, THR))
W <- rbind(OneFactor = fitMeasures(w_uni, ix), ThreeFactor = fitMeasures(w_3f, ix),
           Measurement = fitMeasures(w_meas, ix), ANX_THR_one = fitMeasures(w_one, ix),
           ANX_THR_two = fitMeasures(w_two, ix))
cat("\n== B. WLSMV fit ==\n"); print(round(W, 3))
cl3w <- lavInspect(w_3f, "cor.lv"); cat("Facet correlations (WLSMV):\n"); print(round(cl3w, 3))
dt <- lavTestLRT(w_one, w_two)                 # scaled difference test for WLSMV
cat("Difference test, one vs two factors:\n"); print(dt)
sw <- standardizedSolution(w_meas)
lw <- sw[sw$op == "=~" & sw$lhs == "ANX", "est.std"]
rw <- sw[sw$op == "~~" & sw$lhs == "ANX" & sw$rhs == "THR", c("est.std", "ci.lower", "ci.upper")]
cat("ANX loadings (WLSMV):", round(range(lw), 2), "| ANX-THR r =", round(unlist(rw), 3), "\n")
R$wlsmv <- list(
  fit = lapply(split(round(W, 3), rownames(W)), as.list),
  fit_rows = rownames(W), fit_table = round(W, 3),
  facet_r = round(cl3w[lower.tri(cl3w)], 3),
  anx_loadings_range = round(range(lw), 2),
  anx_thr_r = round(unlist(rw), 3),
  diff_chisq = round(dt[2, "Chisq diff"], 2), diff_df = dt[2, "Df diff"], diff_p = signif(dt[2, "Pr(>Chisq)"], 3))

## ---- C. ESEM: three correlated facets with cross-loadings allowed ----
esem <- paste("efa('b')*F1 + efa('b')*F2 + efa('b')*F3 =~", paste(ANX, collapse = " + "))
R$esem <- list()
for (rot in c("geomin", "oblimin")) {
  f_esem <- try(suppressWarnings(sem(esem, dat, estimator = "MLR", rotation = rot)), silent = TRUE)
  ok <- !inherits(f_esem, "try-error") && lavInspect(f_esem, "converged")
  cat("\n== C. ESEM (3 factors,", rot, "rotation): converged =", ok, "==\n")
  if (ok) {
    fe <- fitMeasures(f_esem, c("chisq.scaled", "df.scaled", "pvalue.scaled", "cfi.robust", "tli.robust", "rmsea.robust", "srmr"))
    le <- lavInspect(f_esem, "std")$lambda; ce <- lavInspect(f_esem, "cor.lv")
    primary <- apply(abs(le), 1, which.max)
    print(round(fe, 3)); print(round(le, 2)); print(round(ce, 2)); cat("Primary factor per item:", primary, "\n")
    R$esem[[rot]] <- list(converged = TRUE, fit = round(fe, 3), loadings = round(unclass(le), 2),
                          factor_r = round(ce[lower.tri(ce)], 2), primary = unname(primary))
  } else R$esem[[rot]] <- list(converged = FALSE)
}

## ---- D. Latent mean differences under scalar invariance ----
lat_diff <- function(g, d) {
  f <- cfa(m_uni, d, estimator = "MLR", group = g, group.equal = c("loadings", "intercepts"))
  pe <- parameterEstimates(f)
  k  <- pe[pe$op == "~1" & pe$lhs == "ANX" & pe$group == 2, ]
  v  <- pe[pe$op == "~~" & pe$lhs == "ANX" & pe$rhs == "ANX", "est"]
  list(groups = lavInspect(f, "group.label"), ns = unname(lavInspect(f, "nobs")),
       diff = round(k$est, 3), se = round(k$se, 3), p = signif(k$pvalue, 3),
       d = round(k$est / sqrt(mean(v)), 2))
}
LD <- list(discipline = lat_diff("grp_disc", dat), year = lat_diff("grp_year", dat),
           gender = lat_diff("grp_gender", subset(dat, !is.na(grp_gender))))
cat("\n== D. Latent mean differences (group 2 minus group 1; d = difference / pooled latent SD) ==\n")
for (n in names(LD)) cat(n, ":", paste(LD[[n]]$groups, collapse = " vs "), "| diff =", LD[[n]]$diff,
                         "SE =", LD[[n]]$se, "p =", LD[[n]]$p, "d =", LD[[n]]$d, "\n")
R$latent_means <- LD

## ---- E. MIMIC: do group differences in ANX survive control for THR? ----
dat$female <- ifelse(dat$gioi_tinh == "Nữ", 1, ifelse(dat$gioi_tinh == "Nam", 0, NA))
dat$year1  <- as.integer(dat$nam_hoc == "Năm 1")
dat$other_disc <- as.integer(dat$khoi_nganh == "Khác")
dat$no_ai_training <- as.integer(dat$da_hoc_AI == "Chưa")
mm <- subset(dat, !is.na(female))
cov <- "female + year1 + other_disc + no_ai_training"
mimic_thr <- paste(two, paste("THR ~", cov), paste("ANX ~ THR +", cov), sep = "\n")
mimic_nothr <- paste(two, paste("THR ~", cov), paste("ANX ~", cov), "ANX ~~ THR", sep = "\n")
f_m1 <- sem(mimic_nothr, mm, estimator = "MLR")
f_m2 <- sem(mimic_thr, mm, estimator = "MLR")
take <- function(f, lhs) { s <- standardizedSolution(f); s <- s[s$op == "~" & s$lhs == lhs, c("rhs", "est.std", "pvalue")]
  s$est.std <- round(s$est.std, 3); s$pvalue <- signif(s$pvalue, 3); s }
cat("\n== E. MIMIC (n =", nrow(mm), ") ==\n")
cat("ANX on covariates, THR not controlled:\n"); print(take(f_m1, "ANX"), row.names = FALSE)
cat("THR on covariates:\n"); print(take(f_m2, "THR"), row.names = FALSE)
cat("ANX on covariates, THR controlled:\n"); print(take(f_m2, "ANX"), row.names = FALSE)
r2 <- lavInspect(f_m2, "r2")
cat("R2:", round(r2, 3), "\n"); print(round(fitMeasures(f_m2, c("cfi.robust", "tli.robust", "rmsea.robust", "srmr")), 3))
R$mimic <- list(n = nrow(mm), anx_no_thr = take(f_m1, "ANX"), thr = take(f_m2, "THR"),
                anx_with_thr = take(f_m2, "ANX"), r2 = round(r2, 3),
                fit = round(fitMeasures(f_m2, c("cfi.robust", "tli.robust", "rmsea.robust", "srmr")), 3))

sink()
write_json(R, file.path(OUT_TXT, "so_lieu_bai_bao.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message(">>> Supplementary analyses saved: ", out, " and so_lieu_bai_bao.json")
