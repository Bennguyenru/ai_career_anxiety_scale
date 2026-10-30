# =============== 05. ADDITIONAL VALIDITY AND PRECISION EVIDENCE ===============
# A. Differential relations: do ANX and THR relate differently to AI self-efficacy and career
#    self-efficacy? (latent correlations constrained equal, scaled chi-square difference tests)
# B. Item response (graded, normal-ogive) parameters from the ordinal one-factor model, test
#    information and conditional reliability across the latent continuum.
# C. Percentile norms for the ANX mean score.
# Appends to ket_qua/so_lieu_bai_bao.json (key "chat_luong") and writes supplementary tables/figure.

suppressMessages({ library(lavaan); library(jsonlite) })
dat <- read.csv(DATA_FILE, fileEncoding = "UTF-8-BOM")
ANX <- paste0("ANX", 1:9); THR <- paste0("THR", 1:3); ASE <- paste0("ASE", 1:3); CSE <- paste0("CSE", 1:3)
Q <- list(N = nrow(dat))
sink(file.path(OUT_TXT, "Ket_qua_chat_luong.txt"), split = TRUE)
cat("ADDITIONAL VALIDITY / PRECISION | N =", nrow(dat), "\n")

## ---- A. Differential relations ----
base <- paste(paste("ANX =~", paste(ANX, collapse = " + ")), paste("THR =~", paste(THR, collapse = " + ")),
              paste("ASE =~", paste(ASE, collapse = " + ")), paste("CSE =~", paste(CSE, collapse = " + ")), sep = "\n")
m_free <- paste(base, "ANX ~~ a*ASE\nTHR ~~ b*ASE\nANX ~~ c*CSE\nTHR ~~ d*CSE", sep = "\n")
f_free <- cfa(m_free, dat, estimator = "MLR", std.lv = TRUE)
pe <- parameterEstimates(f_free)
getr <- function(l, r) { x <- pe[pe$lhs == l & pe$rhs == r & pe$op == "~~", ]; round(c(r = x$est, lo = x$ci.lower, hi = x$ci.upper), 3) }
rr <- list(ANX_ASE = getr("ANX", "ASE"), THR_ASE = getr("THR", "ASE"), ANX_CSE = getr("ANX", "CSE"), THR_CSE = getr("THR", "CSE"))
print(rr)
f_eqASE <- cfa(paste(m_free, "a == b", sep = "\n"), dat, estimator = "MLR", std.lv = TRUE)
f_eqCSE <- cfa(paste(m_free, "c == d", sep = "\n"), dat, estimator = "MLR", std.lv = TRUE)
tA <- lavTestLRT(f_free, f_eqASE); tC <- lavTestLRT(f_free, f_eqCSE)
print(tA); print(tC)
pick <- function(t) list(dchi = round(t[2, "Chisq diff"], 2), df = t[2, "Df diff"], p = signif(t[2, "Pr(>Chisq)"], 3))
Q$diff_rel <- list(cors = rr, eq_ASE = pick(tA), eq_CSE = pick(tC))

## ---- B. Graded (normal-ogive) item parameters from ordinal CFA ----
f_ord <- cfa(paste("ANX =~", paste(ANX, collapse = " + ")), dat, ordered = ANX, estimator = "WLSMV", std.lv = TRUE)
ps <- parameterEstimates(f_ord, standardized = TRUE)
lam <- ps$std.all[ps$op == "=~"]; names(lam) <- ANX
th <- lapply(ANX, function(v) ps$est[ps$op == "|" & ps$lhs == v]); names(th) <- ANX
a_no <- lam / sqrt(1 - lam^2)                      # normal-ogive slope
b <- lapply(ANX, function(v) th[[v]] / sqrt(1 - lam[v]^2)); names(b) <- ANX   # category boundaries on theta
a_logit <- 1.702 * a_no
theta <- seq(-3, 3, by = 0.05)
item_info <- function(a, bb, t) {                 # graded model information (normal ogive)
  P <- cbind(1, sapply(bb, function(bk) pnorm(a * (t - bk))), 0)      # P(X >= k)
  D <- cbind(0, sapply(bb, function(bk) a * dnorm(a * (t - bk))), 0)  # derivatives
  pc <- P[, -ncol(P)] - P[, -1]; dc <- D[, -ncol(D)] - D[, -1]
  rowSums(dc^2 / pmax(pc, 1e-12))
}
I_items <- sapply(ANX, function(v) item_info(a_no[v], b[[v]], theta))
TI <- rowSums(I_items); relc <- TI / (TI + 1)
rng <- range(theta[relc >= .80])
cat("\nSlopes (logistic metric):\n"); print(round(a_logit, 2))
cat("Boundary range:", round(min(unlist(b)), 2), "to", round(max(unlist(b)), 2), "\n")
cat("Conditional reliability >= .80 for theta in", rng, "; peak info", round(max(TI), 2), "at", theta[which.max(TI)], "\n")
cat("Share of sample with theta in that range (from factor scores):\n")
fs <- lavPredict(f_ord); share <- mean(fs >= rng[1] & fs <= rng[2]); print(round(share, 3))
Q$irt <- list(a_logit = as.list(round(a_logit, 2)), a_min = round(min(a_logit), 2), a_max = round(max(a_logit), 2),
              b_min = round(min(unlist(b)), 2), b_max = round(max(unlist(b)), 2),
              rel80_lo = rng[1], rel80_hi = rng[2], peak_info = round(max(TI), 2), peak_theta = theta[which.max(TI)],
              share_in_range = round(share, 3))
tabS7 <- data.frame(Item = ANX, Loading = round(lam, 2), a = round(a_logit, 2),
                    t(sapply(b, function(x) round(x, 2))))
names(tabS7)[4:7] <- paste0("b", 1:4)
write.csv(tabS7, file.path(DIR_BANG, "TableS7_IRT.csv"), row.names = FALSE)
png(file.path(DIR_HINH, "FigureS4_information.png"), width = 1600, height = 1000, res = 220)
par(mar = c(4.2, 4.2, 1, 4.2))
plot(theta, TI, type = "l", lwd = 2, xlab = expression(paste("AI career anxiety (", theta, ")")), ylab = "Test information", las = 1)
par(new = TRUE); plot(theta, relc, type = "l", lty = 2, axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
axis(4, las = 1); mtext("Conditional reliability", side = 4, line = 3); abline(h = .80, col = "grey60", lty = 3)
legend("bottom", c("Information", "Reliability"), lty = c(1, 2), lwd = c(2, 1), bty = "n", horiz = TRUE)
dev.off()

## ---- C. Percentile norms ----
dat$ANXm <- rowMeans(dat[ANX])
sc <- sort(unique(round(dat$ANXm, 2)))
pr <- sapply(sc, function(s) round(100 * mean(dat$ANXm <= s + 1e-9)))
norms <- data.frame(Score = sc, Percentile = pr)
cuts <- quantile(dat$ANXm, c(.10, .25, .50, .75, .90))
print(round(cuts, 2))
Q$norms <- list(p10 = round(unname(cuts[1]), 2), p25 = round(unname(cuts[2]), 2), p50 = round(unname(cuts[3]), 2),
                p75 = round(unname(cuts[4]), 2), p90 = round(unname(cuts[5]), 2))
write.csv(norms, file.path(DIR_BANG, "TableS8_norms.csv"), row.names = FALSE)
sink()

js <- file.path(OUT_TXT, "so_lieu_bai_bao.json")
J <- if (file.exists(js)) fromJSON(js, simplifyVector = FALSE) else list()
J$chat_luong <- Q
write(toJSON(J, auto_unbox = TRUE, pretty = TRUE, digits = NA), js)
cat("05 done\n")
