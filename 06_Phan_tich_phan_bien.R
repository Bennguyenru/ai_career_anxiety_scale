# =============== 06. ANALYSES REQUESTED IN PRE-SUBMISSION REVIEW ===============
# A. Polychoric EFA and parallel analysis (half A) as a check on the Pearson-based EFA
# B. Two-factor ESEM of the 12 ANX and THR items: item-level overlap (cross-loadings)
# C. Local dependence in the ordinal one-factor model (residual correlations)
# D. Ordinal invariance (WLSMV, theta parameterization): thresholds, loadings, intercepts
# E. MIMIC coefficients with 95% CIs, before and after controlling THR
suppressMessages({ library(lavaan); library(semTools); library(psych); library(jsonlite) })
dat <- read.csv(DATA_FILE, fileEncoding = "UTF-8-BOM")
ANX <- paste0("ANX", 1:9); THR <- paste0("THR", 1:3)
Z <- list()
sink(file.path(OUT_TXT, "Ket_qua_phan_bien.txt"), split = TRUE)

## ---- A. Polychoric EFA in half A (same seed as the main split) ----
set.seed(2026); idx <- sample(nrow(dat)); hA <- dat[idx[1:(nrow(dat) / 2)], ]
pc <- polychoric(hA[ANX])$rho
ev <- eigen(pc)$values
set.seed(1); pa <- suppressMessages(fa.parallel(hA[ANX], cor = "poly", fa = "fa", fm = "ml", plot = FALSE, n.iter = 200))
e1 <- fa(hA[ANX], nfactors = 1, fm = "ml", cor = "poly")
cat("A. Polychoric eigenvalues:", round(ev[1:3], 2), "| PA factors:", pa$nfact, "| loadings", round(range(e1$loadings), 2), "| var", round(e1$Vaccounted[2, 1], 3), "\n")
Z$poly <- list(ev = round(ev[1:3], 2), pa = pa$nfact, lmin = round(min(e1$loadings), 2), lmax = round(max(e1$loadings), 2), var = round(e1$Vaccounted[2, 1], 3))

## ---- B. ESEM, two factors, ANX + THR items (geomin), full sample ----
es <- paste0('efa("b")*F1 + efa("b")*F2 =~ ', paste(c(ANX, THR), collapse = " + "))
f_es <- cfa(es, dat, estimator = "MLR", rotation = "geomin")
L <- lavInspect(f_es, "std")$lambda[c(ANX, THR), ]
# name the factors: the one where THR items load highest is THR
thrcol <- which.max(colSums(L[THR, , drop = FALSE]))
anxcol <- setdiff(1:2, thrcol)
Lt <- data.frame(item = rownames(L), ANX_f = round(L[, anxcol], 2), THR_f = round(L[, thrcol], 2))
print(Lt, row.names = FALSE)
phi <- lavInspect(f_es, "std")$psi[1, 2]
fm <- fitMeasures(f_es, c("cfi.robust", "rmsea.robust", "srmr"))
cat("ESEM factor r =", round(phi, 2), "| fit", round(fm, 3), "\n")
Z$esem <- list(loadings = Lt, r = round(phi, 2), fit = as.list(round(fm, 3)),
               anx_cross_max = max(abs(Lt$THR_f[1:9])), anx_cross_items = Lt$item[1:9][abs(Lt$THR_f[1:9]) >= .20])
write.csv(Lt, file.path(DIR_BANG, "TableS10_ESEM_ANX_THR.csv"), row.names = FALSE)

## ---- C. Local dependence, ordinal one-factor ANX model ----
f_ord <- cfa(paste("ANX =~", paste(ANX, collapse = " + ")), dat, ordered = ANX, estimator = "WLSMV", std.lv = TRUE)
rc <- residuals(f_ord, type = "cor")$cov; rv <- rc[lower.tri(rc)]
cat("C. Residual correlations: max |r| =", round(max(abs(rv)), 3), "| n > .10 =", sum(abs(rv) > .10), "of", length(rv), "\n")
Z$localdep <- list(max = round(max(abs(rv)), 3), n_gt10 = sum(abs(rv) > .10), n = length(rv),
                   fit = as.list(round(fitMeasures(f_ord, c("cfi.scaled", "rmsea.scaled", "srmr")), 3)))

## ---- D. Ordinal invariance (Wu & Estabrook, 2016) ----
dat$grp_disc <- ifelse(dat$khoi_nganh == "Công nghệ thông tin", "IT", "Non-IT")
dat$grp_year <- ifelse(dat$nam_hoc %in% c("Năm 1", "Năm 2"), "Junior", "Senior")
dat$grp_gender <- ifelse(dat$gioi_tinh == "Nam", "Male", ifelse(dat$gioi_tinh == "Nữ", "Female", NA))
m1 <- paste("ANX =~", paste(ANX, collapse = " + "))
ordinv <- function(g) {
  d <- dat[!is.na(dat[[g]]), ]
  steps <- list(configural = NULL, thresholds = "thresholds", loadings = c("thresholds", "loadings"),
                intercepts = c("thresholds", "loadings", "intercepts"))
  out <- t(sapply(names(steps), function(s) {
    syn <- as.character(measEq.syntax(m1, data = d, ordered = ANX, parameterization = "theta", ID.fac = "std.lv",
                                      ID.cat = "Wu.Estabrook.2016", group = g, group.equal = steps[[s]]))
    f <- cfa(syn, d, ordered = ANX, estimator = "WLSMV", parameterization = "theta", group = g)
    fitMeasures(f, c("chisq.scaled", "df.scaled", "cfi.scaled", "rmsea.scaled"))
  }))
  out <- cbind(out, dCFI = c(NA, diff(out[, "cfi.scaled"])))
  round(out, 3)
}
OI <- lapply(c(discipline = "grp_disc", year = "grp_year", gender = "grp_gender"), ordinv)
print(OI)
Z$ordinv <- lapply(OI, function(x) list(cfi = as.list(x[, "cfi.scaled"]), dcfi_max = round(max(abs(x[-1, "dCFI"])), 3)))

## ---- E. MIMIC with CIs ----
two <- paste(m1, paste("THR =~", paste(THR, collapse = " + ")), sep = "\n")
dat$female <- ifelse(dat$gioi_tinh == "Nữ", 1, ifelse(dat$gioi_tinh == "Nam", 0, NA))
dat$year1 <- as.integer(dat$nam_hoc == "Năm 1"); dat$other_disc <- as.integer(dat$khoi_nganh == "Khác")
dat$no_ai_training <- as.integer(dat$da_hoc_AI == "Chưa")
mm <- subset(dat, !is.na(female)); cv <- "female + year1 + other_disc + no_ai_training"
f_a <- sem(paste(two, paste("THR ~", cv), paste("ANX ~", cv), "ANX ~~ THR", sep = "\n"), mm, estimator = "MLR")
f_b <- sem(paste(two, paste("THR ~", cv), paste("ANX ~ THR +", cv), sep = "\n"), mm, estimator = "MLR")
ci <- function(f, lhs) { s <- standardizedSolution(f); s <- s[s$op == "~" & s$lhs == lhs, c("rhs", "est.std", "ci.lower", "ci.upper", "pvalue")]
  s[, 2:4] <- round(s[, 2:4], 2); s$pvalue <- round(s$pvalue, 3); s }
M <- list(anx_no = ci(f_a, "ANX"), thr = ci(f_b, "THR"), anx_thr = ci(f_b, "ANX"))
print(M)
Z$mimic_ci <- M
sink()
js <- file.path(OUT_TXT, "so_lieu_bai_bao.json")
J <- if (file.exists(js)) fromJSON(js, simplifyVector = FALSE) else list()
J$phan_bien <- Z
write(toJSON(J, auto_unbox = TRUE, pretty = TRUE, digits = NA, dataframe = "rows"), js)
cat("06 done\n")
