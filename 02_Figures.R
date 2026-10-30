# =================== 02. PUBLICATION FIGURES (English, APA style) ===================
# Requires objects created by 01_Phan_tich.R (dat, f_uni, f_meas, f1, f2, H, B4, B6).
suppressMessages({library(ggplot2); library(patchwork)})
has_ragg <- requireNamespace("ragg", quietly=TRUE)
invisible(file.remove(list.files(DIR_HINH, full.names=TRUE)))

# ---------- palette (validated categorical order; blue <-> red diverging) ----------
C_BLUE <- "#2a78d6"; C_ORANGE <- "#eb6834"; C_AQUA <- "#1baf7a"
INK <- "#1d1d1b"; INK2 <- "#52514e"; INK3 <- "#8a8984"; GRID <- "#ebeae6"
BLUE_LT <- "#cde2fb"; BLUE_MD <- "#86b6ef"; BLUE_DK <- "#1c5cab"
RED_DK <- "#c94141"; RED_LT <- "#f2aaa6"; NEUTRAL <- "#d9d8d4"

theme_pub <- function(base=11){
  theme_minimal(base_size=base) %+replace% theme(
    plot.title=element_text(face="bold", size=base+2, colour=INK, hjust=0, margin=margin(b=4)),
    plot.subtitle=element_text(size=base-0.5, colour=INK2, hjust=0, lineheight=1.15, margin=margin(b=10)),
    plot.caption=element_text(size=base-2, colour=INK2, hjust=0, lineheight=1.2, margin=margin(t=10)),
    plot.title.position="plot", plot.caption.position="plot",
    axis.text=element_text(colour=INK2, size=base-1.5), axis.title=element_text(colour=INK2, size=base-1),
    panel.grid.major=element_line(colour=GRID, linewidth=.35), panel.grid.minor=element_blank(),
    legend.position="top", legend.justification="left",
    legend.title=element_text(size=base-1.5, colour=INK2), legend.text=element_text(size=base-1.5, colour=INK2),
    strip.text=element_text(face="bold", colour=INK, hjust=0, size=base),
    plot.background=element_rect(fill="white", colour=NA), panel.background=element_rect(fill="white", colour=NA),
    plot.margin=margin(14,16,12,14))
}
d2  <- function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.2f", x))          # APA: no leading zero
d3  <- function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.3f", x))
fp  <- function(p) if(p < .001) "p < .001" else paste0("p = ", d3(p))
save_fig <- function(p, name, w, h){
  f <- file.path(DIR_HINH, paste0(name, ".png"))
  if(has_ragg) ggsave(f, p, width=w, height=h, dpi=300, device=ragg::agg_png, bg="white")
  else ggsave(f, p, width=w, height=h, dpi=300, bg="white")
  try(ggsave(file.path(DIR_HINH, paste0(name, ".pdf")), p, width=w, height=h, device=cairo_pdf, bg="white"), silent=TRUE)
}

# ---------- labels (short English paraphrases; replace with the paper's official wording) ----------
fac_name <- c(ANX="AI career anxiety", THR="Threat appraisal", ASE="AI self-efficacy", CSE="Career self-efficacy")
item_lab <- c(
  ANX1="Skills may become obsolete before I use them",
  ANX2="My competencies lose value as AI does similar tasks",
  ANX3="My degree and skills may be undervalued",
  ANX4="Job openings in my field will shrink",
  ANX5="Employers prefer AI over new graduates",
  ANX6="Job hunting after graduation makes me anxious",
  ANX7="I cannot add enough value compared with AI",
  ANX8="I may be less competitive in an automated market",
  ANX9="I may not keep up with employers' AI demands",
  THR1="AI could harm my career prospects", THR2="AI threatens my career goals",
  THR3="AI may take key opportunities in my field",
  ASE1="I can use AI tools for study and work", ASE2="I can learn a new AI tool when needed",
  ASE3="I can master common AI tools",
  CSE1="I can build a successful career", CSE2="I will find a suitable job",
  CSE3="I can meet the demands of my field")
facet_of <- setNames(rep(c("Skill obsolescence","Opportunity contraction","Competitive inadequacy"), each=3), ANX)
facet_col <- c("Skill obsolescence"=C_BLUE, "Opportunity contraction"=C_ORANGE, "Competitive inadequacy"=C_AQUA)
N <- nrow(dat)

fm <- fitMeasures(f_meas, c("chisq.scaled","df.scaled","cfi.robust","tli.robust","rmsea.robust","srmr"))
fit_txt <- sprintf("χ²(%d) = %.2f, CFI = %s, TLI = %s, RMSEA = %s, SRMR = %s, N = %d",
  as.integer(fm["df.scaled"]), fm["chisq.scaled"], d3(fm["cfi.robust"]), d3(fm["tli.robust"]),
  d3(fm["rmsea.robust"]), d3(fm["srmr"]), N)
fu <- fitMeasures(f_uni, c("cfi.robust","rmsea.robust"))

# =====================================================================================
# FIGURE 1. Measurement model (custom path diagram)
# =====================================================================================
blocks <- list(ANX=ANX, THR=THR, ASE=ASE, CSE=CSE)
yy <- c(); cur <- 0
for(f in names(blocks)){ for(it in blocks[[f]]){ yy[it] <- cur; cur <- cur - 1 }; cur <- cur - 1.0 }
it_df <- data.frame(item=names(yy), y=as.numeric(yy), factor=rep(names(blocks), lengths(blocks)))
fy <- tapply(it_df$y, it_df$factor, mean)[names(blocks)]
XF <- 0; FA <- 1.05; FB <- 0.9; XI <- 4.3; HW <- 0.43

ss <- standardizedSolution(f_meas)
ld <- ss[ss$op=="=~", c("lhs","rhs","est.std","pvalue")]
ld$yf <- fy[ld$lhs]; ld$yi <- yy[ld$rhs]
th <- atan2((ld$yi-ld$yf)/FB, (XI-HW-XF)/FA)
ld$x0 <- XF + FA*cos(th); ld$y0 <- ld$yf + FB*sin(th)
ld$x1 <- XI - HW - 0.04;  ld$y1 <- ld$yi
ld$lx <- ld$x0 + .70*(ld$x1-ld$x0); ld$ly <- ld$y0 + .70*(ld$y1-ld$y0)
ld$lw <- .25 + 1.6*pmax(ld$est.std - .40, 0)

ell <- do.call(rbind, lapply(names(blocks), function(f){
  t <- seq(0, 2*pi, length.out=120)
  data.frame(factor=f, x=XF+FA*cos(t), y=fy[[f]]+FB*sin(t))
}))
ell$fill <- ifelse(ell$factor=="ANX", BLUE_LT, "#f1f0ec")
ell$col  <- ifelse(ell$factor=="ANX", C_BLUE, INK3)

ps <- ss[ss$op=="~~" & ss$lhs!=ss$rhs, c("lhs","rhs","est.std","ci.lower","ci.upper")]
ps_all <- ps
ps <- ps[paste(pmin(ps$lhs,ps$rhs), pmax(ps$lhs,ps$rhs))=="ANX THR", ]
arcs <- do.call(rbind, lapply(seq_len(nrow(ps)), function(i){
  ya <- fy[[ps$lhs[i]]]; yb <- fy[[ps$rhs[i]]]; span <- abs(ya-yb); d <- .25 + .16*span
  t <- seq(0, 1, length.out=80); x0 <- XF - FA - .03
  data.frame(id=i, x=x0 - d*sin(pi*t), y=ya + (yb-ya)*t, r=ps$est.std[i],
             key=paste(sort(c(ps$lhs[i], ps$rhs[i])), collapse="-"), lx=x0-d, ly=(ya+yb)/2)
}))
arc_lab <- arcs[!duplicated(arcs$id), ]
arcs$hl <- arcs$key=="ANX-THR"; arc_lab$hl <- arc_lab$key=="ANX-THR"

rel <- setNames(lapply(seq_len(nrow(B4)), function(i) B4[i,]), B4$Scale)
stat_lab <- data.frame(factor=names(blocks), y=fy - FB - .38,
  lab=sapply(names(blocks), function(f) sprintf("α = %s · CR = %s · AVE = %s", d2(rel[[f]]$Alpha), d2(rel[[f]]$CR), d2(rel[[f]]$AVE))))
all_sig <- all(ld$pvalue < .001)
# inset matrix of latent correlations (lower triangle) with 95% CI
cn <- names(blocks); CW <- .68; CH <- 1.05; MX <- -3.4; MY <- fy[["ASE"]] + .9
cell <- do.call(rbind, lapply(2:4, function(i) do.call(rbind, lapply(1:(i-1), function(j){
  a <- cn[i]; b <- cn[j]
  r <- ps_all[(ps_all$lhs==a & ps_all$rhs==b) | (ps_all$lhs==b & ps_all$rhs==a), ][1,]
  data.frame(x=MX + (j-.5)*CW, y=MY - (i-1.5)*CH, r=r$est.std, lo=r$ci.lower, hi=r$ci.upper, hl=(a=="THR" & b=="ANX"))
}))))
col_hd <- data.frame(x=MX + (1:3-.5)*CW, y=MY + .15, lab=cn[1:3])
row_hd <- data.frame(x=MX - .12, y=MY - (2:4-1.5)*CH, lab=cn[2:4])

p1 <- ggplot() +
  geom_path(data=arcs, aes(x, y, group=id, linewidth=hl, colour=hl),
            arrow=arrow(length=unit(1.6,"mm"), ends="both", type="closed"), lineend="round") +
  geom_label(data=arc_lab, aes(lx, ly, label=d2(r), fontface=ifelse(hl,"bold","plain"), colour=hl),
             size=3.1, label.size=0, fill="white", label.padding=unit(.12,"lines")) +
  scale_linewidth_manual(values=c(`FALSE`=.35, `TRUE`=1.3), guide="none") +
  scale_colour_manual(values=c(`FALSE`=INK3, `TRUE`=C_ORANGE), guide="none") +
  geom_segment(data=ld, aes(x=x0, y=y0, xend=x1, yend=y1), linewidth=ld$lw, colour=INK2,
               arrow=arrow(length=unit(1.7,"mm"), type="closed")) +
  geom_label(data=ld, aes(lx, ly, label=d2(est.std)), size=2.8, colour=INK, label.size=0,
             fill="white", label.padding=unit(.08,"lines")) +
  geom_tile(data=cell, aes(x, y, fill=abs(r)), width=CW-.05, height=CH-.08, colour=NA) +
  geom_tile(data=cell[cell$hl,], aes(x, y), width=CW-.05, height=CH-.08, fill=NA, colour=C_ORANGE, linewidth=1) +
  scale_fill_gradient(low="#f2f6fc", high="#86b6ef", limits=c(0,1), guide="none") +
  geom_text(data=cell, aes(x, y+.17, label=d2(r)), size=3.1, fontface="bold", colour=INK) +
  geom_text(data=cell, aes(x, y-.25, label=sprintf("[%s, %s]", d2(lo), d2(hi))), size=2.05, colour=INK2) +
  geom_text(data=col_hd, aes(x, y, label=lab), size=2.9, fontface="bold", colour=INK2, vjust=0) +
  geom_text(data=row_hd, aes(x, y, label=lab), size=2.9, fontface="bold", colour=INK2, hjust=1) +
  annotate("text", x=MX - .55, y=MY + 1.05, label="Latent correlations [95% CI]", hjust=0, size=3.1, fontface="bold", colour=INK) +
  geom_polygon(data=ell, aes(x, y, group=factor), fill=ell$fill, colour=ell$col, linewidth=.6) +
  geom_text(data=data.frame(f=names(fy), y=as.numeric(fy)), aes(XF, y+.2, label=f),
            fontface="bold", size=4.4, colour=INK) +
  geom_text(data=data.frame(f=names(fy), y=as.numeric(fy)), aes(XF, y-.28, label=fac_name[f]),
            size=2.55, colour=INK2) +
  geom_text(data=stat_lab, aes(XF, y, label=lab), size=2.45, colour=INK2) +
  geom_tile(data=it_df, aes(XI, y), width=2*HW, height=.74, fill="#f7f6f3", colour=INK3, linewidth=.35) +
  geom_text(data=it_df, aes(XI, y, label=item), size=3, fontface="bold", colour=INK) +
  geom_text(data=it_df, aes(XI+HW+.18, y, label=item_lab[item]), size=2.85, hjust=0, colour=INK2) +
  coord_cartesian(xlim=c(-3.95, 10.6), ylim=c(min(it_df$y)-1.5, .6), clip="off", expand=FALSE) +
  labs(title="Measurement model with standardized estimates",
       subtitle=paste0("Four-factor CFA (MLR). ", fit_txt),
       caption=paste0("Note. Straight arrows show standardized loadings",
                      if(all_sig) " (all p < .001)" else "",
                      ". All factors were free to correlate (inset, lower left). The orange arc and outlined\n",
                      "cell mark the ANX–THR overlap examined in the discriminant validity analysis. Reliability is shown under each factor.")) +
  theme_void(base_size=11) + theme(
    plot.title=element_text(face="bold", size=13, colour=INK, margin=margin(b=4)),
    plot.subtitle=element_text(size=10.5, colour=INK2, margin=margin(b=12)),
    plot.caption=element_text(size=9, colour=INK2, hjust=0, lineheight=1.2, margin=margin(t=10)),
    plot.title.position="plot", plot.caption.position="plot",
    plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(16,18,12,18))
save_fig(p1, "Fig1_measurement_model", 11, 10)

# =====================================================================================
# FIGURE 2. Item response distribution (diverging) + loadings with 95% CI
# =====================================================================================
lu <- standardizedSolution(f_meas); lu <- lu[lu$op=="=~" & lu$lhs=="ANX", c("rhs","est.std","ci.lower","ci.upper")]
agree <- sapply(ANX, function(it) mean(dat[[it]] >= 4))
ord <- names(sort(agree))                      # lowest agreement at bottom
ypos <- setNames(seq_along(ord), ord)
wrapl <- function(s, w=34) paste(strwrap(s, w), collapse="\n")
ylabs <- sapply(ord, function(it) paste0(it, "  ", wrapl(item_lab[[it]])))

lik <- do.call(rbind, lapply(ANX, function(it){
  p <- as.numeric(table(factor(dat[[it]], levels=1:5)))/N
  left <- -(p[1]+p[2]+p[3]/2)
  xmin <- left + c(0, p[1], p[1]+p[2], p[1]+p[2]+p[3], p[1]+p[2]+p[3]+p[4])
  data.frame(item=it, resp=factor(1:5), xmin=xmin, xmax=xmin+p, y=ypos[[it]])
}))
resp_lab <- c("1 Strongly disagree","2 Disagree","3 Neutral","4 Agree","5 Strongly agree")
resp_col <- c(RED_DK, RED_LT, NEUTRAL, BLUE_MD, C_BLUE)
side <- do.call(rbind, lapply(ANX, function(it){
  p <- as.numeric(table(factor(dat[[it]], levels=1:5)))/N
  data.frame(y=ypos[[it]], dis=p[1]+p[2], agr=p[4]+p[5], xl=-(p[1]+p[2]+p[3]/2), xr=p[3]/2+p[4]+p[5])
}))

p2a <- ggplot(lik) +
  geom_rect(aes(xmin=xmin, xmax=xmax, ymin=y-.36, ymax=y+.36, fill=resp), colour="white", linewidth=.5) +
  geom_vline(xintercept=0, colour=INK3, linewidth=.4) +
  geom_text(data=side, aes(x=xl-.02, y=y, label=paste0(round(100*dis),"%")), hjust=1, size=3, colour=INK2) +
  geom_text(data=side, aes(x=xr+.02, y=y, label=paste0(round(100*agr),"%")), hjust=0, size=3, colour=INK, fontface="bold") +
  scale_fill_manual(values=resp_col, labels=resp_lab, name=NULL) +
  scale_x_continuous(limits=c(-.44,.86), breaks=seq(-.25,.75,.25),
                     labels=function(b) paste0(abs(round(100*b)),"%")) +
  scale_y_continuous(breaks=seq_along(ord), labels=ylabs, expand=expansion(add=.6)) +
  labs(x="Share of respondents (disagree ← | → agree)", y=NULL, subtitle="Response distribution") +
  theme_pub() + theme(panel.grid.major.y=element_blank(), legend.position="bottom",
                      plot.subtitle=element_text(face="bold", colour=INK, size=11),
                      axis.text.y=element_text(hjust=0, lineheight=.95, colour=INK, size=9.2)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))

lu$y <- ypos[lu$rhs]; lu$facet <- factor(facet_of[lu$rhs], levels=names(facet_col))
p2b <- ggplot(lu, aes(y=y)) +
  annotate("rect", xmin=.70, xmax=1, ymin=-Inf, ymax=Inf, fill="#f3f7fd") +
  geom_vline(xintercept=.50, linetype="dashed", colour=INK3, linewidth=.4) +
  geom_vline(xintercept=.70, linetype="dotted", colour=INK3, linewidth=.5) +
  geom_linerange(aes(xmin=ci.lower, xmax=ci.upper, colour=facet), linewidth=1.1, lineend="round") +
  geom_point(aes(x=est.std, colour=facet), size=3.2) +
  geom_point(aes(x=est.std), size=1.2, colour="white") +
  geom_text(aes(x=.93, label=d2(est.std)), hjust=0, size=3.1, colour=INK, fontface="bold") +
  scale_colour_manual(values=facet_col, name="Intended facet") +
  scale_x_continuous(limits=c(.35,1.0), breaks=c(.4,.5,.6,.7,.8), labels=d2) +
  scale_y_continuous(breaks=seq_along(ord), labels=NULL, expand=expansion(add=.6)) +
  labs(x="Standardized loading (95% CI)", y=NULL, subtitle="Standardized loading") +
  theme_pub() + theme(panel.grid.major.y=element_blank(), legend.position="bottom",
                      plot.subtitle=element_text(face="bold", colour=INK, size=11)) +
  guides(colour=guide_legend(ncol=1, title.position="top", override.aes=list(size=2.6)))

p2 <- (p2a | p2b) + plot_layout(widths=c(2.5,1)) +
  plot_annotation(
    title="AI career anxiety items: how students responded and how strongly each item loads",
    subtitle=sprintf("Items sorted by agreement. Agreement ranges from %d%% to %d%% across items, and loadings are similar across the three\nintended facets, which is consistent with a single underlying dimension.", round(100*min(agree)), round(100*max(agree))),
    caption=paste0("Note. N = ", N, ". Percentages left and right of each bar are total disagreement (1–2) and agreement (4–5); the neutral share is split at zero.\n",
                   "Loadings are from the four-factor measurement model (Figure 1). Dashed line = .50, dotted line = .70. One-factor model for the 9 items: CFI = ", d3(fu["cfi.robust"]),
                   ", RMSEA = ", d3(fu["rmsea.robust"]), "."),
    theme=theme_pub())
save_fig(p2, "Fig2_items_and_loadings", 12.5, 7.6)

# =====================================================================================
# FIGURE 3. Discriminant validity: HTMT heatmap + latent correlations with 95% CI
# =====================================================================================
Hm <- as.matrix(H); cons <- rownames(Hm)
hd <- do.call(rbind, lapply(seq_along(cons), function(i) lapply(seq_along(cons), function(j)
  if(i > j) data.frame(row=cons[i], col=cons[j], v=Hm[i,j]) else NULL)) |> unlist(recursive=FALSE))
hd$row <- factor(hd$row, levels=rev(cons)); hd$col <- factor(hd$col, levels=cons)
hd$flag <- hd$v >= .85
p3a <- ggplot(hd, aes(col, row)) +
  geom_tile(aes(fill=v), colour="white", linewidth=1.6) +
  geom_tile(data=hd[hd$flag,], fill=NA, colour=C_ORANGE, linewidth=1.3, width=.94, height=.94) +
  geom_text(aes(label=d2(v), colour=v > .55, fontface=ifelse(flag,"bold","plain")), size=4.2) +
  scale_fill_gradient(low="#eef4fc", high=BLUE_DK, limits=c(0,1), breaks=c(0,.5,1), labels=d2, name="HTMT") +
  scale_colour_manual(values=c(`FALSE`=INK, `TRUE`="white"), guide="none") +
  scale_x_discrete(position="top") +
  labs(x=NULL, y=NULL, subtitle="(a) HTMT ratios (orange outline ≥ .85)") +
  theme_pub() + theme(panel.grid.major=element_blank(), legend.position="right",
                      axis.text=element_text(face="bold", colour=INK, size=10.5),
                      plot.subtitle=element_text(face="bold", colour=INK, size=11))

lc <- ps_all; lc$pair <- paste(lc$lhs, "–", lc$rhs); lc$hl <- lc$est.std >= .85
lc <- lc[order(lc$est.std), ]; lc$pair <- factor(lc$pair, levels=lc$pair)
p3b <- ggplot(lc, aes(y=pair)) +
  annotate("rect", xmin=.85, xmax=1.0, ymin=-Inf, ymax=Inf, fill="#fdf0ea") +
  geom_vline(xintercept=1, colour=INK2, linewidth=.5) +
  geom_vline(xintercept=.85, colour=C_ORANGE, linetype="dashed", linewidth=.45) +
  geom_linerange(aes(xmin=ci.lower, xmax=ci.upper, colour=hl), linewidth=1.1, lineend="round") +
  geom_point(aes(x=est.std, colour=hl), size=3.2) +
  geom_text(aes(x=1.1, label=sprintf("%s  [%s, %s]", d2(est.std), d2(ci.lower), d2(ci.upper)), fontface=ifelse(hl,"bold","plain")),
            hjust=0, size=3.2, colour=INK) +
  annotate("text", x=1.1, y=nrow(lc)+.55, label="r  [95% CI]", hjust=0, size=3, colour=INK2) +
  scale_colour_manual(values=c(`FALSE`=C_BLUE, `TRUE`=C_ORANGE), guide="none") +
  scale_x_continuous(limits=c(-.15,1.5), breaks=c(0,.25,.5,.75,1), labels=d2) +
  scale_y_discrete(expand=expansion(add=c(.6,1))) + coord_cartesian(clip="off") +
  labs(x="Latent correlation (95% CI)", y=NULL, subtitle="(b) Latent factor correlations") +
  theme_pub() + theme(panel.grid.major.y=element_blank(),
                      axis.text.y=element_text(face="bold", colour=INK, size=10),
                      plot.subtitle=element_text(face="bold", colour=INK, size=11))

an <- as.data.frame(anova(f1, f2))
sb <- sprintf("Δχ²(%d) = %.2f, %s", as.integer(an[2,"Df diff"]), an[2,"Chisq diff"], fp(an[2,"Pr(>Chisq)"]))
rr <- lc[lc$pair=="ANX – THR", ]
p3 <- (p3a | p3b) + plot_layout(widths=c(1,1.25)) +
  plot_annotation(
    title="Discriminant validity: AI career anxiety is distinct from, but close to, threat appraisal",
    subtitle=paste0("ANX and THR share a high latent correlation (r = ", d2(rr$est.std), ", 95% CI [", d2(rr$ci.lower), ", ", d2(rr$ci.upper),
                    "]), yet the interval excludes 1.0 and a two-factor\nmodel fits better than a single merged factor (", sb, ")."),
    caption="Note. HTMT computed from absolute item correlations (Henseler et al., 2015). The shaded band in panel (b) marks the .85 caution zone; a confidence\ninterval that excludes 1.0 supports discriminant validity. Dashed orange line = .85, solid line = 1.0.",
    theme=theme_pub())
save_fig(p3, "Fig3_discriminant_validity", 12, 5.6)

# =====================================================================================
# FIGURE 4. Measurement invariance (separate panels; no dual axis)
# =====================================================================================
iv <- B6; iv$Step <- factor(tools::toTitleCase(iv$Step), levels=c("Configural","Metric","Scalar"))
short <- c("Discipline (IT vs Non-IT)"="Discipline", "Year (Junior vs Senior)"="Year", "Gender (Male vs Female)"="Gender")
iv$Grouping <- factor(iv$Grouping, levels=unique(iv$Grouping))
gcol <- setNames(c(C_BLUE, C_ORANGE, C_AQUA), levels(iv$Grouping))
pd <- position_dodge(width=.14)
mk_inv <- function(var, ttl, ylim, good, goodlab, dirn){
  endp <- iv[iv$Step=="Scalar", ]
  ggplot(iv, aes(Step, .data[[var]], colour=Grouping, group=Grouping)) +
    scale_x_discrete(expand=expansion(add=c(.35,.75))) + geom_blank() +
    annotate("rect", xmin=-Inf, xmax=Inf,
             ymin=if(dirn=="above") good else -Inf, ymax=if(dirn=="above") Inf else good, fill="#eef7f2") +
    geom_hline(yintercept=good, linetype="dashed", colour=INK3, linewidth=.4) +
    annotate("text", x=.7, y=good, label=goodlab, hjust=0, vjust=if(dirn=="above") 1.6 else -.7, size=3, colour=INK2) +
    geom_line(linewidth=1, position=pd) +
    geom_point(size=3.2, position=pd) + geom_point(size=1.2, colour="white", position=pd) +
    geom_text(data=endp, aes(label=short[as.character(Grouping)]), position=position_nudge(x=.16),
              hjust=0, size=3.1, show.legend=FALSE) +
    scale_colour_manual(values=gcol, name=NULL) +
    scale_y_continuous(limits=ylim, labels=d3) +
    labs(x=NULL, y=NULL, subtitle=ttl) +
    theme_pub() + theme(plot.subtitle=element_text(face="bold", colour=INK, size=11))
}
p4a <- mk_inv("cfi.robust", "(a) Robust CFI (higher is better)", c(.94, 1.004), .95, "Good fit ≥ .95", "above")
p4b <- mk_inv("rmsea.robust", "(b) Robust RMSEA (lower is better)", c(0, .08), .06, "Good fit ≤ .06", "below")
dc <- iv[iv$Step!="Configural", c("Grouping","Step","dCFI")]
dc$Step <- ifelse(dc$Step=="Metric", "Configural → Metric", "Metric → Scalar")
dc$row <- paste0(short[as.character(dc$Grouping)], "\n", dc$Step)
dc$row <- factor(dc$row, levels=rev(dc$row))
p4c <- ggplot(dc, aes(dCFI, row, colour=Grouping)) +
  scale_y_discrete() + geom_blank() +
  annotate("rect", xmin=-Inf, xmax=.010, ymin=-Inf, ymax=Inf, fill="#eef7f2") +
  geom_vline(xintercept=0, colour=INK3, linewidth=.4) +
  geom_vline(xintercept=.010, linetype="dashed", colour=C_ORANGE, linewidth=.5) +
  annotate("text", x=.0105, y=.6, label="Cutoff\n.010", hjust=0, vjust=0, size=2.9, colour=C_ORANGE, lineheight=.9) +
  geom_segment(aes(x=0, xend=dCFI, yend=row), linewidth=1) +
  geom_point(size=3.2) + geom_point(size=1.2, colour="white") +
  geom_text(aes(x=pmax(dCFI,0)+.0012, label=d3(dCFI)), hjust=0, size=3, colour=INK) +
  scale_colour_manual(values=gcol, guide="none") +
  scale_x_continuous(limits=c(-.008,.016), breaks=c(-.005,0,.005,.010), labels=d3) +
  labs(x="Drop in CFI after adding constraints", y=NULL, subtitle="(c) ΔCFI between steps") +
  theme_pub() + theme(plot.subtitle=element_text(face="bold", colour=INK, size=11),
                      panel.grid.major.y=element_blank(), axis.text.y=element_text(colour=INK, size=8.8, lineheight=.95))
p4 <- (p4a | p4b | p4c) + plot_layout(widths=c(1,1,1.05), guides="collect") & theme(legend.position="bottom")
p4 <- p4 + plot_annotation(
  title="Measurement invariance of the AI career anxiety scale across student groups",
  subtitle="Fit stays good as loadings (metric) and then intercepts (scalar) are held equal, and no constraint lowers CFI by more than .010.\nThe scale therefore measures the same construct on the same metric in every group, so group means can be compared.",
  caption=paste0("Note. Multi-group CFA with MLR on the 9-item one-factor model. Groups: Discipline = IT vs non-IT; Year = Years 1–2 vs Years 3+; Gender = male vs female\n",
                 "(respondents reporting another gender or preferring not to say excluded). ΔCFI cutoff follows Chen (2007); negative values mean fit improved."),
  theme=theme_pub() + theme(legend.position="bottom"))
save_fig(p4, "Fig4_measurement_invariance", 13, 5.8)

# =====================================================================================
# FIGURE 5. Group comparisons of the anxiety composite (raincloud-style)
# =====================================================================================
gm <- mean(dat$ANX_mean); gsd <- sd(dat$ANX_mean)
yr <- tr(dat$nam_hoc, map_nam); yr[yr %in% c("Year 4", ">4 years")] <- "Year 4+"
gd <- tr(dat$gioi_tinh, map_gioi); gd[!gd %in% c("Male","Female")] <- NA
ds <- tr(dat$khoi_nganh, map_nganh)
long <- rbind(
  data.frame(panel="Year of study", g=yr, s=dat$ANX_mean),
  data.frame(panel="Gender", g=gd, s=dat$ANX_mean),
  data.frame(panel="Discipline", g=ds, s=dat$ANX_mean))
long <- long[!is.na(long$g), ]
lev <- list(`Year of study`=c("Year 1","Year 2","Year 3","Year 4+"), Gender=c("Female","Male"),
            Discipline=c("IT","Business/Management","Design/Media","Other"))
test_txt <- function(d){
  d$g <- factor(d$g)
  if(nlevels(d$g)==2){
    t <- t.test(s ~ g, data=d); a <- split(d$s, d$g)
    sp <- sqrt(((length(a[[1]])-1)*var(a[[1]]) + (length(a[[2]])-1)*var(a[[2]]))/(length(d$s)-2))
    sprintf("t(%.0f) = %.2f, %s, d = %s", t$parameter, abs(t$statistic), fp(t$p.value), d2(abs(diff(sapply(a,mean)))/sp))
  } else {
    w <- oneway.test(s ~ g, data=d); av <- summary(aov(s ~ g, data=d))[[1]]
    eta <- av[1,"Sum Sq"]/sum(av[,"Sum Sq"])
    sprintf("F(%d, %.0f) = %.2f, %s, η² = %s", as.integer(w$parameter[1]), w$parameter[2], w$statistic, fp(w$p.value), d2(eta))
  }
}
mk_grp <- function(pn){
  d <- long[long$panel==pn, ]; d$g <- factor(d$g, levels=lev[[pn]])
  sm <- do.call(rbind, lapply(split(d$s, d$g), function(x)
    data.frame(m=mean(x), lo=mean(x)-qt(.975,length(x)-1)*sd(x)/sqrt(length(x)),
               hi=mean(x)+qt(.975,length(x)-1)*sd(x)/sqrt(length(x)), n=length(x))))
  sm$g <- factor(rownames(sm), levels=lev[[pn]])
  xl <- setNames(paste0(gsub("/", "/\n", levels(d$g)), "\nn = ", sm$n[match(levels(d$g), sm$g)]), levels(d$g))
  set.seed(7)
  ggplot(d, aes(g, s)) +
    geom_hline(yintercept=gm, linetype="dashed", colour=INK3, linewidth=.4) +
    geom_violin(fill=BLUE_LT, colour=NA, width=.85, alpha=.9, scale="width") +
    geom_jitter(width=.13, height=0, size=.75, alpha=.35, colour=BLUE_DK) +
    geom_boxplot(width=.16, outlier.shape=NA, fill="white", colour=INK2, linewidth=.4, alpha=.85) +
    geom_linerange(data=sm, aes(x=g, ymin=lo, ymax=hi), inherit.aes=FALSE, colour=C_ORANGE, linewidth=1.2,
                   position=position_nudge(x=.24)) +
    geom_point(data=sm, aes(x=g, y=m), inherit.aes=FALSE, colour=C_ORANGE, size=2.8, shape=18,
               position=position_nudge(x=.24)) +
    geom_label(data=sm, aes(x=g, y=m, label=sprintf("%.2f", m)), inherit.aes=FALSE, size=2.9, colour=INK,
               fontface="bold", position=position_nudge(x=.33), hjust=0, label.size=0, fill="white",
               label.padding=unit(.08,"lines")) +
    scale_x_discrete(labels=xl, expand=expansion(add=c(.55,.7))) +
    scale_y_continuous(limits=c(1,5.3), breaks=1:5) +
    labs(x=NULL, y=if(pn=="Year of study") "AI career anxiety (mean of 9 items, 1–5)" else NULL,
         subtitle=pn, title=NULL, caption=NULL, tag=NULL) +
    annotate("text", x=.45, y=5.22, label=test_txt(d), hjust=0, size=3, colour=INK2) +
    theme_pub() + theme(panel.grid.major.x=element_blank(),
                        plot.subtitle=element_text(face="bold", colour=INK, size=10.5, lineheight=1.25),
                        axis.text.x=element_text(colour=INK, size=9, lineheight=.95))
}
p5 <- (mk_grp("Year of study") | mk_grp("Gender") | mk_grp("Discipline")) + plot_layout(widths=c(4,2.6,4)) +
  plot_annotation(
    title="AI career anxiety across student groups",
    subtitle=sprintf("Overall M = %.2f, SD = %.2f (dashed line). Scalar invariance holds for each grouping, so these mean comparisons are interpretable.", gm, gsd),
    caption=paste0("Note. Violins show the distribution, boxes the median and IQR, dots individual students; orange diamonds and bars are means with 95% CI.\n",
                   "Tests are Welch's ANOVA (F) or Welch's t test (t) on the composite score; η² and Cohen's d are effect sizes. ",
                   "Year 4 and beyond-year-4 students\nare combined (n = ", sum(dat$nam_hoc=="Trên 4 năm"), " in the latter). ",
                   "The gender panel excludes ", sum(is.na(gd)), " respondents reporting another gender or preferring not to say."),
    theme=theme_pub())
save_fig(p5, "Fig5_group_comparisons", 12.5, 6.4)

message(">>> Figures saved to ", DIR_HINH)
print(list.files(DIR_HINH))
