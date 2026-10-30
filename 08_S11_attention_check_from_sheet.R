# 08_S11_attention_check_from_sheet.R
# Sensitivity of the main results to the attention-check exclusion, using the
# response sheet without Timestamp/E-mail ("Raw Data", 24 columns, 589 rows).
#
# The analytic sample is taken as the 464 rows of the public cleaned file
# (data/ai_career_anxiety_oct2026_n464.csv); every one of them is found in the
# sheet. Failures (ATTN != 2) and straightliners are identified from the answers.
# Because the sheet has no e-mail column, duplicate submissions cannot be
# identified, so the failure group may contain one duplicate (83 here vs. 82 in
# Table S1).

suppressMessages({library(lavaan); library(psych); library(semTools)})
set.seed(2026)
dir.create("output", showWarnings = FALSE)

SHEET_FILE <- "data/Tong_hop_form.xlsx"
PUB_FILE   <- "data/ai_career_anxiety_oct2026_n464.csv"

pub <- read.csv(PUB_FILE, stringsAsFactors = FALSE, encoding = "UTF-8")
cols <- names(pub)[-1]
sh <- as.data.frame(readxl::read_excel(SHEET_FILE, sheet = "Raw Data"))
sh <- sh[rowSums(!is.na(sh)) > 0, 1:24]
names(sh) <- cols

anx <- paste0("ANX", 1:9); thr <- paste0("THR", 1:3)
ase <- paste0("ASE", 1:3); cse <- paste0("CSE", 1:3)
items <- c(thr, anx, ase, cse)
grp <- c("khoi_nganh", "nam_hoc", "gioi_tinh", "tan_suat_AI", "da_hoc_AI")
for (v in c(items, "ATTN")) sh[[v]] <- as.numeric(sh[[v]])

sh$fail <- sh$ATTN != 2
sh$sl   <- apply(sh[, items], 1, function(r) length(unique(r)) == 1)

retained <- pub[, cols]
failed   <- sh[sh$fail, cols]
sl_pass  <- sh[!sh$fail & sh$sl, cols]
failed_sl <- sh[sh$fail & sh$sl, cols]
cat("Sheet rows:", nrow(sh), "| failed attention check:", nrow(failed),
    "(of which straightliners:", nrow(failed_sl), ") | straightliners who passed:", nrow(sl_pass), "\n")
stopifnot(nrow(retained) == 464, nrow(sl_pass) == 31)

# ---------- 1. Who failed the check? ----------
dd <- rbind(cbind(failed, status = "failed"), cbind(retained, status = "retained"))
dd$ANX_mean <- rowMeans(dd[, anx]); dd$THR_mean <- rowMeans(dd[, thr])
comp <- list()
for (g in grp) {
  tab <- table(dd$status, dd[[g]])
  ct <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 10000))
  V <- sqrt(unname(suppressWarnings(chisq.test(tab))$statistic) / (sum(tab) * (min(dim(tab)) - 1)))
  cat("\n==", g, "== (row %)\n"); print(round(100 * prop.table(tab, 1), 1))
  comp[[g]] <- data.frame(variable = g, chisq = unname(ct$statistic), p_sim = ct$p.value, cramers_V = V)
}
comp <- do.call(rbind, comp)
for (v in c("ANX_mean", "THR_mean")) {
  a <- dd[dd$status == "failed", v]; b <- dd[dd$status == "retained", v]
  tt <- t.test(a, b)
  d_ <- (mean(a) - mean(b)) / sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  comp <- rbind(comp, data.frame(variable = paste0(v, sprintf(" (failed %.2f [%.2f] vs retained %.2f [%.2f])",
                                  mean(a), sd(a), mean(b), sd(b))),
                                 chisq = unname(tt$statistic), p_sim = tt$p.value, cramers_V = d_))
}
names(comp)[2:4] <- c("statistic", "p", "effect (V or d)")
print(comp, digits = 3, row.names = FALSE)
write.csv(comp, "output/S11a_failed_vs_retained.csv", row.names = FALSE)

# ---------- 2. Main results with the exclusions undone ----------
m1 <- paste("ANX =~", paste(anx, collapse = " + "))
m3 <- paste("S =~ ANX1 + ANX2 + ANX3", "O =~ ANX4 + ANX5 + ANX6", "C =~ ANX7 + ANX8 + ANX9", sep = "\n")
mT <- paste("THR =~", paste(thr, collapse = " + "))
m4 <- paste(m1, mT, paste("ASE =~", paste(ase, collapse = " + ")),
            paste("CSE =~", paste(cse, collapse = " + ")), sep = "\n")
m2 <- paste(m1, mT, sep = "\n")
mM <- paste("F =~", paste(c(anx, thr), collapse = " + "))

fitc <- function(model, dat) {
  w <- character(0)
  f <- withCallingHandlers(cfa(model, dat, estimator = "MLR", std.lv = TRUE),
         warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  st <- c(if (!lavInspect(f, "converged")) "not converged",
          if (!isTRUE(suppressWarnings(lavInspect(f, "post.check")))) "inadmissible", w)
  attr(f, "status") <- if (length(st)) paste(unique(st), collapse = "; ") else "ok"
  f
}
fi <- function(f, m) unname(fitMeasures(f, m))
lr <- function(f, a, b) { s <- standardizedSolution(f); s[s$lhs == a & s$op == "~~" & s$rhs == b, ] }
omega <- function(f, fac) { s <- standardizedSolution(f); l <- s$est.std[s$lhs == fac & s$op == "=~"]
  sum(l)^2 / (sum(l)^2 + sum(1 - l^2)) }
ave <- function(f, fac) { s <- standardizedSolution(f); mean(s$est.std[s$lhs == fac & s$op == "=~"]^2) }

run <- function(dat, label) {
  f1 <- fitc(m1, dat); f3 <- fitc(m3, dat); f4 <- fitc(m4, dat); f2 <- fitc(m2, dat); fM <- fitc(mM, dat)
  s3 <- standardizedSolution(f3); fr <- s3$est.std[s3$op == "~~" & s3$lhs != s3$rhs]
  r  <- lr(f4, "ANX", "THR")
  dt <- lavTestLRT(f2, fM, method = "satorra.bentler.2001")
  ht <- htmt(m4, dat)
  data.frame(
    sample = label, N = nrow(dat),
    one_factor_CFI_RMSEA_SRMR = sprintf("%.3f / %.3f / %.3f", fi(f1, "cfi.robust"), fi(f1, "rmsea.robust"), fi(f1, "srmr")),
    three_factor_max_facet_r = sprintf("%.2f", max(fr)),
    four_factor_CFI_RMSEA_SRMR = sprintf("%.3f / %.3f / %.3f", fi(f4, "cfi.robust"), fi(f4, "rmsea.robust"), fi(f4, "srmr")),
    alpha_ANX = sprintf("%.3f", psych::alpha(dat[, anx])$total$raw_alpha),
    omega_ANX = sprintf("%.3f", omega(f4, "ANX")), AVE_ANX = sprintf("%.3f", ave(f4, "ANX")),
    r_ANX_THR = sprintf("%.2f [%.2f, %.2f]", r$est.std, r$ci.lower, r$ci.upper),
    HTMT_ANX_THR = sprintf("%.3f", ht["THR", "ANX"]),
    dchisq_1v2 = sprintf("%.2f, p %s", dt[2, "Chisq diff"], format.pval(dt[2, "Pr(>Chisq)"], digits = 2, eps = .001)),
    r_ANX_ASE_CSE = sprintf("%.2f / %.2f", lr(f4, "ANX", "ASE")$est.std, lr(f4, "ANX", "CSE")$est.std),
    ANX_M_SD = sprintf("%.2f (%.2f)", mean(rowMeans(dat[, anx])), sd(rowMeans(dat[, anx]))),
    status = paste(attr(f1, "status"), attr(f3, "status"), attr(f4, "status"), attr(f2, "status"), attr(fM, "status"), sep = " | "),
    check.names = FALSE)
}

res <- rbind(
  run(retained, "Analytic sample"),
  run(rbind(retained, failed[!apply(failed[, items], 1, function(r) length(unique(r)) == 1), ]),
      "Attention-check step undone (straightliners still removed)"),
  run(rbind(retained, failed, sl_pass), "Both steps undone"))
print(t(res))
write.csv(res, "output/S11b_main_results_by_sample.csv", row.names = FALSE)
