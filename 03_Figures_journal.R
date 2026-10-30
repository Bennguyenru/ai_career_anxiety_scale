# ============ 03. JOURNAL FIGURE SET (drawn at final print size, APA 7) ============
# Requires objects from 01_Phan_tich.R (dat, f_uni, f_meas, f1, f2, H, B4, B6, tr, map_*).
# Every figure is drawn at its real printed width, so no scaling is needed:
#   full page width = JOURNAL_WIDTH_MM (default 165 mm, fits US Letter and A4), 1.5 column = 140 mm. Smallest text = 8 pt.
# Titles, subtitles and notes are NOT inside the images (APA 7 puts them in the manuscript);
# they are written to ket_qua/hinh_journal/Figure_captions_APA.md with the current numbers.
suppressMessages({library(lavaan); library(semTools); library(ggplot2); library(patchwork)})
has_ragg <- requireNamespace("ragg", quietly=TRUE)
dir.create(DIR_HINH_J, recursive=TRUE, showWarnings=FALSE)
invisible(file.remove(list.files(DIR_HINH_J, full.names=TRUE)))

MM <- 1/25.4                      # mm -> inch
FULLW <- if (exists("JOURNAL_WIDTH_MM")) JOURNAL_WIDTH_MM else 165   # full text width (165 mm fits US Letter and A4 journals)
PT <- function(pt) pt/ggplot2::.pt  # point size -> ggplot text size
C_BLUE <- "#2a78d6"; C_ORANGE <- "#eb6834"; C_AQUA <- "#1baf7a"
INK <- "#1d1d1b"; INK2 <- "#4a4945"; INK3 <- "#8a8984"; GRID <- "#e6e5e1"
BLUE_LT <- "#cde2fb"; BLUE_MD <- "#86b6ef"; BLUE_DK <- "#1c5cab"
RED_DK <- "#c94141"; RED_LT <- "#f2aaa6"; NEUTRAL <- "#d9d8d4"

d2 <- function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.2f", x))
d3 <- function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.3f", x))
fp <- function(p) if(p < .001) "p < .001" else paste0("p = ", d3(p))
fpi <- function(p) sub("^p", "*p*", fp(p))

theme_j <- function(){
  theme_minimal(base_size=9) %+replace% theme(
    axis.text=element_text(colour=INK2, size=8), axis.title=element_text(colour=INK2, size=8.5),
    panel.grid.major=element_line(colour=GRID, linewidth=.3), panel.grid.minor=element_blank(),
    legend.position="bottom", legend.justification="left",
    legend.title=element_text(size=8, colour=INK2), legend.text=element_text(size=8, colour=INK2),
    legend.key.size=unit(3.2,"mm"), legend.margin=margin(0,0,0,0),
    plot.subtitle=element_text(face="bold", size=9, colour=INK, hjust=0, margin=margin(b=4)),
    plot.title.position="plot",
    plot.background=element_rect(fill="white", colour=NA), panel.background=element_rect(fill="white", colour=NA),
    plot.margin=margin(4,6,4,4))
}
save_j <- function(p, name, w_mm, h_mm){
  w <- w_mm*MM; h <- h_mm*MM; base <- file.path(DIR_HINH_J, name)
  if(has_ragg){
    ggsave(paste0(base,".png"), p, width=w, height=h, dpi=600, device=ragg::agg_png, bg="white")
    try(ggsave(paste0(base,".tiff"), p, width=w, height=h, dpi=600, device=ragg::agg_tiff, bg="white", compression="lzw"), silent=TRUE)
  } else {
    ggsave(paste0(base,".png"), p, width=w, height=h, dpi=600, bg="white")
    try(ggsave(paste0(base,".tiff"), p, width=w, height=h, dpi=600, bg="white", device="tiff", compression="lzw",
               type=if(isTRUE(capabilities("cairo"))) "cairo" else getOption("bitmapType")), silent=TRUE)
  }
  try(ggsave(paste0(base,".pdf"), p, width=w, height=h, device=cairo_pdf, bg="white"), silent=TRUE)
}

fac_name <- c(ANX="AI career anxiety", THR="Threat appraisal", ASE="AI self-efficacy", CSE="Career self-efficacy")
item_lab <- c(   # short English paraphrases; replace with the official English wording of the scale
  ANX1="Skills may become obsolete before I use them",
  ANX2="My competencies lose value as AI does similar tasks",
  ANX3="My degree and skills may be undervalued",
  ANX4="Job openings in my field will shrink",
  ANX5="Employers prefer AI over new graduates",
  ANX6="Job hunting after graduation makes me anxious",
  ANX7="I cannot add enough value compared with AI",
  ANX8="I may be less competitive in an automated market",
  ANX9="I may not keep up with employers' AI demands")
facet_of  <- setNames(rep(c("Skill obsolescence","Opportunity contraction","Competitive inadequacy"), each=3), ANX)
facet_col <- c("Skill obsolescence"=C_BLUE, "Opportunity contraction"=C_ORANGE, "Competitive inadequacy"=C_AQUA)
N <- nrow(dat)
fm <- fitMeasures(f_meas, c("chisq.scaled","df.scaled","cfi.robust","tli.robust","rmsea.robust","srmr"))
fu <- fitMeasures(f_uni,  c("chisq.scaled","df.scaled","cfi.robust","rmsea.robust","srmr"))
ss <- standardizedSolution(f_meas)
ps_all <- ss[ss$op=="~~" & ss$lhs!=ss$rhs, c("lhs","rhs","est.std","ci.lower","ci.upper","pvalue")]
cap <- list()

# =====================================================================================
# FIGURE 1. Measurement model (140 mm x 165 mm)
# =====================================================================================
blocks <- list(ANX=ANX, THR=THR, ASE=ASE, CSE=CSE)
yy <- c(); cur <- 0
for(f in names(blocks)){ for(it in blocks[[f]]){ yy[it] <- cur; cur <- cur - 1 }; cur <- cur - 1.0 }
it_df <- data.frame(item=names(yy), y=as.numeric(yy))
fy <- tapply(it_df$y, rep(names(blocks), lengths(blocks)), mean)[names(blocks)]
XF <- 0; FA <- 1.2; FB <- .95; XI <- 4.4; HW <- .45
ld <- ss[ss$op=="=~", c("lhs","rhs","est.std","pvalue")]
ld$yf <- fy[ld$lhs]; ld$yi <- yy[ld$rhs]
th <- atan2((ld$yi-ld$yf)/FB, (XI-HW-XF)/FA)
ld$x0 <- XF + FA*cos(th); ld$y0 <- ld$yf + FB*sin(th)
ld$x1 <- XI - HW - .04;   ld$y1 <- ld$yi
ld$lx <- ld$x0 + .72*(ld$x1-ld$x0); ld$ly <- ld$y0 + .72*(ld$y1-ld$y0)
ld$lw <- .25 + 1.3*pmax(ld$est.std - .40, 0)
ell <- do.call(rbind, lapply(names(blocks), function(f){
  t <- seq(0, 2*pi, length.out=120); data.frame(factor=f, x=XF+FA*cos(t), y=fy[[f]]+FB*sin(t)) }))
ell$fill <- ifelse(ell$factor=="ANX", BLUE_LT, "#f1f0ec"); ell$col <- ifelse(ell$factor=="ANX", C_BLUE, INK3)
pr <- ps_all[(ps_all$lhs=="ANX" & ps_all$rhs=="THR") | (ps_all$lhs=="THR" & ps_all$rhs=="ANX"), ][1, ]
ya <- fy[["ANX"]]; yb <- fy[["THR"]]; dd <- .3 + .16*abs(ya-yb); tt <- seq(0,1,length.out=80); x0 <- XF-FA-.03
arc <- data.frame(x=x0 - dd*sin(pi*tt), y=ya + (yb-ya)*tt)
cn <- names(blocks); CW <- .78; CH <- 1.0; MX <- -3.7; MY <- fy[["ASE"]] + .8
cell <- do.call(rbind, lapply(2:4, function(i) do.call(rbind, lapply(1:(i-1), function(j){
  a <- cn[i]; b <- cn[j]; r <- ps_all[(ps_all$lhs==a & ps_all$rhs==b) | (ps_all$lhs==b & ps_all$rhs==a), ][1,]
  data.frame(x=MX+(j-.5)*CW, y=MY-(i-1.5)*CH, r=r$est.std, hl=(a=="THR" & b=="ANX")) }))))
p1 <- ggplot() +
  geom_path(data=arc, aes(x, y), colour=C_ORANGE, linewidth=1, arrow=arrow(length=unit(1.4,"mm"), ends="both", type="closed")) +
  geom_label(aes(x=x0-dd, y=(ya+yb)/2, label=d2(pr$est.std)), colour=C_ORANGE, fontface="bold", size=PT(8.5),
             label.size=0, fill="white", label.padding=unit(.1,"lines")) +
  geom_segment(data=ld, aes(x=x0, y=y0, xend=x1, yend=y1), linewidth=ld$lw, colour=INK2,
               arrow=arrow(length=unit(1.3,"mm"), type="closed")) +
  geom_label(data=ld, aes(lx, ly, label=d2(est.std)), size=PT(8), colour=INK, label.size=0, fill="white",
             label.padding=unit(.06,"lines")) +
  geom_polygon(data=ell, aes(x, y, group=factor), fill=ell$fill, colour=ell$col, linewidth=.45) +
  geom_text(data=data.frame(f=names(fy), y=as.numeric(fy)), aes(XF, y+.22, label=f), fontface="bold", size=PT(9.5), colour=INK) +
  geom_text(data=data.frame(f=names(fy), y=as.numeric(fy)), aes(XF, y-.3, label=fac_name[f]), size=PT(8), colour=INK2) +
  geom_tile(data=it_df, aes(XI, y), width=2*HW, height=.74, fill="#f7f6f3", colour=INK3, linewidth=.3) +
  geom_text(data=it_df, aes(XI, y, label=item), size=PT(8), fontface="bold", colour=INK) +
  geom_tile(data=cell, aes(x, y, fill=abs(r)), width=CW-.05, height=CH-.08) +
  geom_tile(data=cell[cell$hl,], aes(x, y), width=CW-.05, height=CH-.08, fill=NA, colour=C_ORANGE, linewidth=.8) +
  scale_fill_gradient(low="#f2f6fc", high=BLUE_MD, limits=c(0,1), guide="none") +
  geom_text(data=cell, aes(x, y, label=d2(r)), size=PT(8.5), fontface="bold", colour=INK) +
  geom_text(data=data.frame(x=MX+(1:3-.5)*CW, lab=cn[1:3]), aes(x, MY+.12, label=lab), size=PT(8), fontface="bold", colour=INK2, vjust=0) +
  geom_text(data=data.frame(y=MY-(2:4-1.5)*CH, lab=cn[2:4]), aes(MX-.1, y, label=lab), size=PT(8), fontface="bold", colour=INK2, hjust=1) +
  annotate("text", x=MX-.6, y=MY+.85, label="Latent correlations", hjust=0, size=PT(8.5), fontface="bold", colour=INK) +
  coord_cartesian(xlim=c(-4.35, 5.0), ylim=c(min(it_df$y)-.6, .6), clip="off", expand=FALSE) +
  theme_void() + theme(plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(3,3,3,3))
save_j(p1, "Figure1_measurement_model", 140, 165)
cap$f1 <- sprintf(paste0(
  "**Figure 1**\n\n*Measurement Model of AI Career Anxiety and Related Constructs*\n\n",
  "*Note.* Standardized estimates from a four-factor confirmatory factor analysis with robust maximum likelihood (MLR), *N* = %d; ",
  "χ²(%d) = %.2f, CFI = %s, TLI = %s, RMSEA = %s, SRMR = %s. Straight arrows are standardized loadings%s. ",
  "All factors were free to correlate; the inset gives the latent correlations, and the curved arrow marks the ",
  "correlation between ANX and THR examined in the discriminant validity analysis. ANX = AI career anxiety; ",
  "THR = threat appraisal; ASE = AI self-efficacy; CSE = career self-efficacy."),
  N, as.integer(fm["df.scaled"]), fm["chisq.scaled"], d3(fm["cfi.robust"]), d3(fm["tli.robust"]),
  d3(fm["rmsea.robust"]), d3(fm["srmr"]), if(all(ld$pvalue < .001)) " (all *p* < .001)" else "")

# =====================================================================================
# FIGURE 2. Item responses + loadings (174 mm x 120 mm)
# =====================================================================================
lu <- ss[ss$op=="=~" & ss$lhs=="ANX", c("rhs","est.std","ci.lower","ci.upper")]
agree <- sapply(ANX, function(it) mean(dat[[it]] >= 4))
ord <- names(sort(agree)); ypos <- setNames(seq_along(ord), ord)
wrapl <- function(s, w=30) paste(strwrap(s, w), collapse="\n")
ylabs <- sapply(ord, function(it) paste0(it, "  ", wrapl(item_lab[[it]])))
lik <- do.call(rbind, lapply(ANX, function(it){
  p <- as.numeric(table(factor(dat[[it]], levels=1:5)))/N; left <- -(p[1]+p[2]+p[3]/2)
  xmin <- left + c(0, cumsum(p)[1:4]); data.frame(resp=factor(1:5), xmin=xmin, xmax=xmin+p, y=ypos[[it]]) }))
side <- do.call(rbind, lapply(ANX, function(it){
  p <- as.numeric(table(factor(dat[[it]], levels=1:5)))/N
  data.frame(y=ypos[[it]], dis=p[1]+p[2], agr=p[4]+p[5], xl=-(p[1]+p[2]+p[3]/2), xr=p[3]/2+p[4]+p[5]) }))
p2a <- ggplot(lik) +
  geom_rect(aes(xmin=xmin, xmax=xmax, ymin=y-.36, ymax=y+.36, fill=resp), colour="white", linewidth=.35) +
  geom_vline(xintercept=0, colour=INK3, linewidth=.35) +
  geom_text(data=side, aes(x=xl-.015, y=y, label=paste0(round(100*dis),"%")), hjust=1, size=PT(8), colour=INK2) +
  geom_text(data=side, aes(x=xr+.015, y=y, label=paste0(round(100*agr),"%")), hjust=0, size=PT(8), colour=INK, fontface="bold") +
  scale_fill_manual(values=c(RED_DK, RED_LT, NEUTRAL, BLUE_MD, C_BLUE), name=NULL,
                    labels=c("Strongly disagree","Disagree","Neutral","Agree","Strongly agree")) +
  scale_x_continuous(limits=c(-.5,.9), breaks=c(-.25,0,.25,.5,.75), labels=function(b) paste0(abs(round(100*b)),"%")) +
  scale_y_continuous(breaks=seq_along(ord), labels=ylabs, limits=c(.4, length(ord)+.6), expand=c(0,0)) +
  labs(x="Share of respondents", y=NULL, subtitle="(a) Response distribution") +
  theme_j() + theme(panel.grid.major.y=element_blank(), axis.text.y=element_text(hjust=0, colour=INK, size=8, lineheight=.95)) +
  guides(fill=guide_legend(nrow=1))
lu$y <- ypos[lu$rhs]; lu$facet <- factor(facet_of[lu$rhs], levels=names(facet_col))
p2b <- ggplot(lu, aes(y=y)) +
  geom_vline(xintercept=.50, linetype="dashed", colour=INK3, linewidth=.35) +
  geom_vline(xintercept=.70, linetype="dotted", colour=INK3, linewidth=.45) +
  geom_linerange(aes(xmin=ci.lower, xmax=ci.upper, colour=facet), linewidth=.9) +
  geom_point(aes(x=est.std, colour=facet), size=2.2) +
  geom_text(aes(x=.95, label=d2(est.std)), hjust=0, size=PT(8), colour=INK, fontface="bold") +
  scale_colour_manual(values=facet_col, name="Intended facet") +
  scale_x_continuous(limits=c(.35,1.03), breaks=c(.4,.6,.8), labels=d2) +
  scale_y_continuous(breaks=seq_along(ord), labels=NULL, limits=c(.4, length(ord)+.6), expand=c(0,0)) +
  labs(x="Loading (95% CI)", y=NULL, subtitle="(b) Loading") +
  theme_j() + theme(panel.grid.major.y=element_blank()) +
  guides(colour=guide_legend(nrow=2, title.position="left", byrow=TRUE))
p2 <- (p2a | p2b) + plot_layout(widths=c(2.55,1), guides="collect") &
  theme(legend.position="bottom", legend.box="vertical", legend.box.just="left", legend.justification="left",
        legend.spacing.y=unit(1,"mm"))
save_j(p2, "Figure2_item_responses_loadings", FULLW, 125)
cap$f2 <- sprintf(paste0(
  "**Figure 2**\n\n*Response Distribution and Standardized Loadings of the Nine AI Career Anxiety Items*\n\n",
  "*Note.* *N* = %d. Panel (a): items are ordered by the share of agreement; numbers to the left and right of each bar are ",
  "total disagreement (1 or 2) and agreement (4 or 5), and the neutral share is split at zero. Panel (b): standardized loadings ",
  "with 95%% confidence intervals from the measurement model in Figure 1, colored by the facet each item was written for; ",
  "dashed line = .50, dotted line = .70. A one-factor model for the nine items fit well, χ²(%d) = %.2f, CFI = %s, RMSEA = %s, SRMR = %s."),
  N, as.integer(fu["df.scaled"]), fu["chisq.scaled"], d3(fu["cfi.robust"]), d3(fu["rmsea.robust"]), d3(fu["srmr"]))

# =====================================================================================
# FIGURE 3. Discriminant validity (174 mm x 72 mm)
# =====================================================================================
Hm <- as.matrix(H); cons <- rownames(Hm)
hd <- do.call(rbind, lapply(2:length(cons), function(i) do.call(rbind, lapply(1:(i-1), function(j)
  data.frame(row=cons[i], col=cons[j], v=Hm[i,j])))))
hd$row <- factor(hd$row, levels=rev(cons)); hd$col <- factor(hd$col, levels=cons); hd$flag <- hd$v >= .85
p3a <- ggplot(hd, aes(col, row)) +
  geom_tile(aes(fill=v), colour="white", linewidth=1.2) +
  geom_tile(data=hd[hd$flag,], fill=NA, colour=C_ORANGE, linewidth=.9, width=.93, height=.93) +
  geom_text(aes(label=d2(v), colour=v > .55, fontface=ifelse(flag,"bold","plain")), size=PT(9)) +
  scale_fill_gradient(low="#eef4fc", high=BLUE_DK, limits=c(0,1), guide="none") +
  scale_colour_manual(values=c(`FALSE`=INK, `TRUE`="white"), guide="none") +
  scale_x_discrete(position="top") + labs(x=NULL, y=NULL, subtitle="(a) HTMT ratio") +
  theme_j() + theme(panel.grid.major=element_blank(), axis.text=element_text(face="bold", colour=INK, size=8.5))
lc <- ps_all; lc$pair <- paste(lc$lhs, "–", lc$rhs); lc$hl <- lc$est.std >= .85
lc <- lc[order(lc$est.std), ]; lc$pair <- factor(lc$pair, levels=lc$pair)
p3b <- ggplot(lc, aes(y=pair)) +
  scale_y_discrete(expand=expansion(add=c(.6,1))) + geom_blank() +
  annotate("rect", xmin=.85, xmax=1, ymin=-Inf, ymax=Inf, fill="#fdf0ea") +
  geom_vline(xintercept=1, colour=INK2, linewidth=.45) +
  geom_vline(xintercept=.85, colour=C_ORANGE, linetype="dashed", linewidth=.4) +
  geom_linerange(aes(xmin=ci.lower, xmax=ci.upper, colour=hl), linewidth=.9) +
  geom_point(aes(x=est.std, colour=hl), size=2.2) +
  geom_text(aes(x=1.07, label=sprintf("%s [%s, %s]", d2(est.std), d2(ci.lower), d2(ci.upper)), fontface=ifelse(hl,"bold","plain")),
            hjust=0, size=PT(8), colour=INK) +
  annotate("text", x=1.07, y=nrow(lc)+.75, label="r [95% CI]", hjust=0, size=PT(8), colour=INK2) +
  scale_colour_manual(values=c(`FALSE`=C_BLUE, `TRUE`=C_ORANGE), guide="none") +
  scale_x_continuous(limits=c(-.15,1.55), breaks=c(0,.25,.5,.75,1), labels=d2) +
  coord_cartesian(clip="off") +
  labs(x="Latent correlation", y=NULL, subtitle="(b) Latent factor correlation") +
  theme_j() + theme(panel.grid.major.y=element_blank(), axis.text.y=element_text(face="bold", colour=INK, size=8))
p3 <- (p3a | p3b) + plot_layout(widths=c(1,1.55))
save_j(p3, "Figure3_discriminant_validity", FULLW, 74)
an <- as.data.frame(anova(f1, f2)); rr <- lc[lc$pair %in% c("ANX – THR","THR – ANX"), ]
cap$f3 <- sprintf(paste0(
  "**Figure 3**\n\n*Discriminant Validity Among AI Career Anxiety, Threat Appraisal, AI Self-Efficacy, and Career Self-Efficacy*\n\n",
  "*Note.* *N* = %d. Panel (a): heterotrait-monotrait ratios (Henseler et al., 2015); the outlined cell is at or above .85. ",
  "Panel (b): latent correlations with 95%% confidence intervals; the shaded band marks the .85 to 1.0 range. ANX and THR correlate ",
  "*r* = %s, 95%% CI [%s, %s]; the interval excludes 1.0, and a two-factor model fit better than a single merged factor, ",
  "Satorra-Bentler Δχ²(%d) = %.2f, %s."),
  N, d2(rr$est.std), d2(rr$ci.lower), d2(rr$ci.upper), as.integer(an[2,"Df diff"]), an[2,"Chisq diff"], fpi(an[2,"Pr(>Chisq)"]))

# =====================================================================================
# FIGURE 4. Measurement invariance (174 mm x 130 mm)
# =====================================================================================
iv <- B6; iv$Step <- factor(tools::toTitleCase(iv$Step), levels=c("Configural","Metric","Scalar"))
short <- c("Discipline (IT vs Non-IT)"="Discipline", "Year (Junior vs Senior)"="Year", "Gender (Male vs Female)"="Gender")
iv$G <- factor(short[as.character(iv$Grouping)], levels=c("Discipline","Year","Gender"))
gcol <- setNames(c(C_BLUE, C_ORANGE, C_AQUA), levels(iv$G))
pd <- position_dodge(width=.14)
mk_inv <- function(var, ttl, ylim, good, goodlab, dirn){
  ggplot(iv, aes(Step, .data[[var]], colour=G, group=G)) +
    scale_x_discrete(expand=expansion(add=c(.3,1.05))) + geom_blank() + coord_cartesian(clip="off") +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=if(dirn=="above") good else -Inf, ymax=if(dirn=="above") Inf else good, fill="#eef7f2") +
    geom_hline(yintercept=good, linetype="dashed", colour=INK3, linewidth=.35) +
    annotate("text", x=.65, y=good, label=goodlab, hjust=0, vjust=if(dirn=="above") 1.6 else -.7, size=PT(8), colour=INK2) +
    geom_line(linewidth=.8, position=pd) + geom_point(size=2.2, position=pd) +
    geom_text(data=iv[iv$Step=="Scalar",], aes(label=G), position=position_nudge(x=.17), hjust=0, size=PT(8), show.legend=FALSE) +
    scale_colour_manual(values=gcol, guide="none") + scale_y_continuous(limits=ylim, labels=d3) +
    labs(x=NULL, y=NULL, subtitle=ttl) + theme_j()
}
p4a <- mk_inv("cfi.robust", "(a) Robust CFI", c(.94, 1.004), .95, "≥ .95", "above")
p4b <- mk_inv("rmsea.robust", "(b) Robust RMSEA", c(0, .08), .06, "≤ .06", "below")
dc <- iv[iv$Step!="Configural", c("G","Step","dCFI")]
dc$lab <- paste0(dc$G, ": ", ifelse(dc$Step=="Metric", "configural → metric", "metric → scalar"))
dc$lab <- factor(dc$lab, levels=rev(dc$lab))
p4c <- ggplot(dc, aes(dCFI, lab, colour=G)) +
  scale_y_discrete() + geom_blank() +
  annotate("rect", xmin=-Inf, xmax=.010, ymin=-Inf, ymax=Inf, fill="#eef7f2") +
  geom_vline(xintercept=0, colour=INK3, linewidth=.35) +
  geom_vline(xintercept=.010, linetype="dashed", colour=C_ORANGE, linewidth=.45) +
  annotate("text", x=.0103, y=6.4, label="cutoff .010", hjust=0, size=PT(8), colour=C_ORANGE) +
  geom_segment(aes(x=0, xend=dCFI, yend=lab), linewidth=.8) + geom_point(size=2.2) +
  geom_text(aes(x=pmax(dCFI,0)+.0008, label=d3(dCFI)), hjust=0, size=PT(8), colour=INK) +
  scale_colour_manual(values=gcol, guide="none") +
  scale_x_continuous(limits=c(-.008,.016), breaks=c(-.005,0,.005,.010), labels=d3) +
  coord_cartesian(clip="off") +
  labs(x="Drop in robust CFI after adding equality constraints", y=NULL, subtitle="(c) Change in CFI between steps") +
  theme_j() + theme(panel.grid.major.y=element_blank(), axis.text.y=element_text(colour=INK, size=8))
p4 <- wrap_elements(full=(p4a | p4b)) / wrap_elements(full=p4c) + plot_layout(heights=c(1, .95))
save_j(p4, "Figure4_measurement_invariance", FULLW, 132)
cap$f4 <- paste0(
  "**Figure 4**\n\n*Measurement Invariance of the AI Career Anxiety Scale Across Discipline, Year of Study, and Gender*\n\n",
  "*Note.* Multigroup confirmatory factor analysis of the nine-item one-factor model with MLR. Discipline = information technology versus other ",
  "disciplines; Year = Years 1 and 2 versus Year 3 and above; Gender = male versus female (respondents who reported another gender or preferred not ",
  "to say were excluded from this comparison). Shaded areas mark good fit. Panel (c) shows the drop in robust CFI when loadings (metric) and then ",
  "intercepts (scalar) were constrained equal; values at or below .010 support invariance (Chen, 2007), and negative values indicate improved fit.")

# =====================================================================================
# FIGURE 5. Group comparisons (174 mm x 95 mm)
# =====================================================================================
gm <- mean(dat$ANX_mean); gsd <- sd(dat$ANX_mean)
yr <- tr(dat$nam_hoc, map_nam); yr[yr %in% c("Year 4", ">4 years")] <- "Year 4+"
gd <- tr(dat$gioi_tinh, map_gioi); gd[!gd %in% c("Male","Female")] <- NA
ds <- tr(dat$khoi_nganh, map_nganh); ds[ds=="Business/Management"] <- "Business"; ds[ds=="Design/Media"] <- "Design"
long <- rbind(data.frame(panel="Year of study", g=yr, s=dat$ANX_mean),
              data.frame(panel="Gender", g=gd, s=dat$ANX_mean),
              data.frame(panel="Discipline", g=ds, s=dat$ANX_mean))
long <- long[!is.na(long$g), ]
lev <- list(`Year of study`=c("Year 1","Year 2","Year 3","Year 4+"), Gender=c("Female","Male"),
            Discipline=c("IT","Business","Design","Other"))
tests <- list()
test_txt <- function(d, pn){
  d$g <- factor(d$g)
  if(nlevels(d$g)==2){
    t <- t.test(s ~ g, data=d); a <- split(d$s, d$g)
    sp <- sqrt(((length(a[[1]])-1)*var(a[[1]]) + (length(a[[2]])-1)*var(a[[2]]))/(length(d$s)-2))
    dd <- abs(diff(sapply(a,mean)))/sp
    tests[[pn]] <<- sprintf("*t*(%.1f) = %.2f, *p* %s, *d* = %s", t$parameter, abs(t$statistic), sub("^p ","",fp(t$p.value)), d2(dd))
    sprintf("t(%.0f) = %.2f, %s, d = %s", t$parameter, abs(t$statistic), fp(t$p.value), d2(dd))
  } else {
    w <- oneway.test(s ~ g, data=d); av <- summary(aov(s ~ g, data=d))[[1]]; eta <- av[1,"Sum Sq"]/sum(av[,"Sum Sq"])
    tests[[pn]] <<- sprintf("*F*(%d, %.1f) = %.2f, *p* %s, η² = %s", as.integer(w$parameter[1]), w$parameter[2], w$statistic, sub("^p ","",fp(w$p.value)), d2(eta))
    sprintf("F(%d, %.0f) = %.2f, %s, η² = %s", as.integer(w$parameter[1]), w$parameter[2], w$statistic, fp(w$p.value), d2(eta))
  }
}
tag <- c(`Year of study`="(a)", Gender="(b)", Discipline="(c)")
mk_grp <- function(pn){
  d <- long[long$panel==pn, ]; d$g <- factor(d$g, levels=lev[[pn]])
  sm <- do.call(rbind, lapply(split(d$s, d$g), function(x){ se <- sd(x)/sqrt(length(x)); q <- qt(.975, length(x)-1)
    data.frame(m=mean(x), lo=mean(x)-q*se, hi=mean(x)+q*se, n=length(x)) }))
  sm$g <- factor(rownames(sm), levels=lev[[pn]])
  xl <- setNames(paste0(gsub("/", "/\n", levels(d$g)), "\nn = ", sm$n[match(levels(d$g), sm$g)]), levels(d$g))
  set.seed(7)
  ggplot(d, aes(g, s)) +
    geom_hline(yintercept=gm, linetype="dashed", colour=INK3, linewidth=.35) +
    geom_violin(fill=BLUE_LT, colour=NA, width=.85, scale="width") +
    geom_jitter(width=.12, height=0, size=.45, alpha=.35, colour=BLUE_DK) +
    geom_boxplot(width=.15, outlier.shape=NA, fill="white", colour=INK2, linewidth=.3, alpha=.85) +
    geom_linerange(data=sm, aes(x=g, ymin=lo, ymax=hi), inherit.aes=FALSE, colour=C_ORANGE, linewidth=.9, position=position_nudge(x=.23)) +
    geom_point(data=sm, aes(x=g, y=m), inherit.aes=FALSE, colour=C_ORANGE, size=2, shape=18, position=position_nudge(x=.23)) +
    scale_x_discrete(labels=xl, expand=expansion(add=c(.5,.5))) +
    scale_y_continuous(limits=c(1,5), breaks=1:5) +
    labs(x=NULL, y=if(pn=="Year of study") "AI career anxiety (mean score, 1–5)" else NULL,
         subtitle=paste0(tag[[pn]], " ", pn, "\n", test_txt(d, pn))) +
    theme_j() + theme(panel.grid.major.x=element_blank(), axis.text.x=element_text(colour=INK, size=8, lineheight=.95),
                      plot.subtitle=element_text(face="plain", size=8, colour=INK, lineheight=1.15, hjust=0, margin=margin(b=4)))
}
p5 <- (mk_grp("Year of study") | mk_grp("Gender") | mk_grp("Discipline")) + plot_layout(widths=c(4,2.9,4))
save_j(p5, "Figure5_group_comparisons", FULLW, 98)
cap$f5 <- sprintf(paste0(
  "**Figure 5**\n\n*AI Career Anxiety by Year of Study, Gender, and Discipline*\n\n",
  "*Note.* Composite score = mean of the nine items (*M* = %.2f, *SD* = %.2f; dashed line). Violins show the distribution, boxes the median and ",
  "interquartile range, and dots individual students; orange diamonds and bars are group means with 95%% confidence intervals. Tests are Welch's ",
  "one-way ANOVA or Welch's t test: year of study, %s; gender, %s; discipline, %s. Year 4 and students beyond Year 4 were combined (*n* = %d in the latter). ",
  "Business = business or management; Design = design or media. The gender panel excludes %d respondents who reported another gender or preferred not to say. Scalar invariance (Figure 4) supports these comparisons."),
  gm, gsd, tests[["Year of study"]], tests[["Gender"]], tests[["Discipline"]], sum(dat$nam_hoc=="Trên 4 năm"), sum(is.na(gd)))

writeLines(c("# Figure captions (APA 7), generated from the current data", "",
             paste0("Data: ", DATA_FILE, ", N = ", N, ". Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), "."), "",
             "Item labels in Figure 2 are short English paraphrases. Replace them with the official English wording before submission (item_lab in 03_Figures_journal.R).", "",
             cap$f1, "", cap$f2, "", cap$f3, "", cap$f4, "", cap$f5),
           file.path(DIR_HINH_J, "Figure_captions_APA.md"), useBytes=FALSE)
message(">>> Journal figures saved to ", DIR_HINH_J)
print(list.files(DIR_HINH_J))
