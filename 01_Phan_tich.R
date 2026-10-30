# =============== 01. SCALE VALIDATION + EXPORT TABLES/FIGURES ===============
options(repos = c(CRAN = "https://cloud.r-project.org"))
need <- function(p){
  if(!requireNamespace(p, quietly=TRUE)){
    message("   installing ", p, " ...")
    tryCatch(install.packages(p, type="binary", quiet=TRUE),
             error=function(e) install.packages(p, quiet=TRUE))
  }
}
for(p in c("lavaan","semTools","psych","ggplot2","patchwork","ragg")) need(p)
suppressMessages({library(lavaan); library(semTools); library(psych); library(ggplot2)})

dir.create(DIR_BANG, recursive=TRUE, showWarnings=FALSE)
dir.create(DIR_HINH, recursive=TRUE, showWarnings=FALSE)

EST <- "MLR"
idx <- c("chisq.scaled","df.scaled","pvalue.scaled","cfi.robust","tli.robust","rmsea.robust","srmr")
sf  <- function(f) tryCatch(round(fitMeasures(f, idx),3),
                            error=function(e) setNames(rep(NA,length(idx)), idx))
wcsv <- function(x, name, rn=TRUE)
  write.csv(x, file.path(DIR_BANG, name), row.names=rn, fileEncoding="UTF-8")
navy <- "#1F3864"

# ---- English recode maps (display only; analysis keeps original values) ----
tr <- function(x, m){ y<-m[as.character(x)]; y[is.na(y)]<-as.character(x)[is.na(y)]; unname(y) }
map_nganh <- c("Công nghệ thông tin"="IT","Kinh doanh / Quản trị"="Business/Management",
               "Thiết kế / Truyền thông"="Design/Media","Khác"="Other")
map_nam   <- c("Năm 1"="Year 1","Năm 2"="Year 2","Năm 3"="Year 3","Năm 4"="Year 4","Trên 4 năm"=">4 years")
map_gioi  <- c("Nam"="Male","Nữ"="Female","Khác"="Other","Không muốn nêu"="Prefer not to say")
map_ts    <- c("Chưa từng"="Never","Hiếm khi"="Rarely","Thỉnh thoảng"="Sometimes","Thường xuyên"="Often","Hằng ngày"="Daily")
map_hoc   <- c("Chưa"="No","Có, không chính thức"="Yes, informal","Có, chính thức trong chương trình"="Yes, formal")

ANX<-paste0("ANX",1:9); THR<-paste0("THR",1:3); ASE<-paste0("ASE",1:3); CSE<-paste0("CSE",1:3)
dat <- read.csv(DATA_FILE, fileEncoding="UTF-8-BOM")
message(">>> Read ", DATA_FILE, " | N = ", nrow(dat))

sink(file.path(OUT_TXT,"Ket_qua_phan_tich.txt"), split=TRUE)
cat("RESULTS: VALIDATION OF THE 'AI CAREER ANXIETY' SCALE\n")
cat("N =", nrow(dat), "| estimator MLR | date", format(Sys.time(),"%Y-%m-%d %H:%M"),"\n")
cat("Data file:", DATA_FILE, "\n")

## ---- TABLE 1: sample description ----
demo_tab <- function(v, label, m=NULL){
  x <- if(is.null(m)) dat[[v]] else tr(dat[[v]], m)
  t<-table(x); data.frame(Variable=label, Group=names(t),
    Count=as.integer(t), Percent=round(100*as.integer(t)/sum(t),1))
}
B1 <- do.call(rbind, list(
  demo_tab("khoi_nganh","Discipline",map_nganh), demo_tab("nam_hoc","Year of study",map_nam),
  demo_tab("gioi_tinh","Gender",map_gioi), demo_tab("tan_suat_AI","AI use frequency",map_ts),
  demo_tab("da_hoc_AI","Studied AI",map_hoc)))
wcsv(B1,"Table1_sample.csv",rn=FALSE); cat("\n== TABLE 1: Sample description ==\n"); print(B1,row.names=FALSE)

## ---- 9-item models + TABLE 3 ----
m_uni<-'ANX=~ANX1+ANX2+ANX3+ANX4+ANX5+ANX6+ANX7+ANX8+ANX9'
m_3f <-'OBS=~ANX1+ANX2+ANX3
        OPP=~ANX4+ANX5+ANX6
        CMP=~ANX7+ANX8+ANX9'
m_2nd<-paste(m_3f,'ANX=~OBS+OPP+CMP',sep="\n")
f_uni<-cfa(m_uni,dat,estimator=EST); f_3f<-cfa(m_3f,dat,estimator=EST); f_2nd<-cfa(m_2nd,dat,estimator=EST)
m_bif<-'g=~ANX1+ANX2+ANX3+ANX4+ANX5+ANX6+ANX7+ANX8+ANX9
        sO=~ANX1+ANX2+ANX3
        sP=~ANX4+ANX5+ANX6
        sC=~ANX7+ANX8+ANX9'
f_bif<-suppressWarnings(tryCatch(cfa(m_bif,dat,estimator=EST,orthogonal=TRUE), error=function(e) NULL))
# admissibility checks: improper solutions are reported, not interpreted
adm_chk <- function(f){
  if(is.null(f) || !lavInspect(f,"converged")) return("No (did not converge)")
  psi <- lavInspect(f,"est")$psi; th <- lavInspect(f,"est")$theta
  if(any(diag(psi) < 0) || any(diag(th) < 0)) return("No (negative variance)")
  cl <- lavInspect(f,"cor.lv"); if(nrow(cl)>1 && max(abs(cl[lower.tri(cl)])) >= 1) return("No (latent r >= 1)")
  "Yes"
}
cl3 <- lavInspect(f_3f,"cor.lv")
B3<-rbind(Unidimensional=sf(f_uni), ThreeFactor=sf(f_3f), SecondOrder=sf(f_2nd),
          Bifactor=if(is.null(f_bif) || !lavInspect(f_bif,"converged")) setNames(rep(NA,length(idx)),idx) else sf(f_bif))
B3 <- data.frame(Model=rownames(B3), B3,
                 Admissible=c(adm_chk(f_uni), adm_chk(f_3f), adm_chk(f_2nd), adm_chk(f_bif)), row.names=NULL)
wcsv(B3,"Table3_model_comparison.csv",rn=FALSE); cat("\n== TABLE 3: Competing models (9 items) ==\n"); print(B3,row.names=FALSE)
cat("Facet correlations in the three-factor model (values near or above 1 mean the facets are not distinct):\n")
print(round(cl3,3))

## ---- full measurement model ----
meas<-'ANX=~ANX1+ANX2+ANX3+ANX4+ANX5+ANX6+ANX7+ANX8+ANX9
       THR=~THR1+THR2+THR3
       ASE=~ASE1+ASE2+ASE3
       CSE=~CSE1+CSE2+CSE3'
f_meas<-cfa(meas,dat,estimator=EST)
cat("\n== Measurement model fit ==\n"); print(sf(f_meas))

## ---- TABLE 2: item descriptives + loadings ----
items<-c(ANX,THR,ASE,CSE)
desc<-psych::describe(dat[,items])
itc<-c()
for(s in list(ANX,THR,ASE,CSE)){
  a<-suppressWarnings(suppressMessages(psych::alpha(dat[,s])))
  itc<-c(itc,setNames(a$item.stats$r.drop,s))
}
std_all<-standardizedSolution(f_meas); std_all<-std_all[std_all$op=="=~",]
load<-setNames(std_all$est.std,std_all$rhs)
B2<-data.frame(Item=items, Mean=round(desc$mean,2), SD=round(desc$sd,2),
  Skew=round(desc$skew,2), Kurtosis=round(desc$kurtosis,2),
  ItemTotal=round(as.numeric(itc[items]),2),
  StdLoading=round(as.numeric(load[items]),2))
wcsv(B2,"Table2_items.csv",rn=FALSE); cat("\n== TABLE 2: Item descriptives and loadings ==\n"); print(B2,row.names=FALSE)

## ---- TABLE 4: reliability, CR, AVE ----
crave<-function(f){l<-std_all$est.std[std_all$lhs==f]; c(CR=sum(l)^2/(sum(l)^2+sum(1-l^2)), AVE=mean(l^2))}
al<-function(v) suppressWarnings(suppressMessages(psych::alpha(dat[,v])$total$raw_alpha))
B4<-data.frame(Scale=c("ANX","THR","ASE","CSE"), Items=c(9,3,3,3),
  Alpha=round(c(al(ANX),al(THR),al(ASE),al(CSE)),3),
  CR=round(c(crave("ANX")[1],crave("THR")[1],crave("ASE")[1],crave("CSE")[1]),3),
  AVE=round(c(crave("ANX")[2],crave("THR")[2],crave("ASE")[2],crave("CSE")[2]),3))
wcsv(B4,"Table4_reliability.csv",rn=FALSE); cat("\n== TABLE 4: Reliability, CR, AVE ==\n"); print(B4,row.names=FALSE)

## ---- TABLE 5: HTMT ----
H<-round(htmt(meas,dat,htmt2=FALSE),3)
wcsv(as.data.frame(unclass(H)),"Table5_HTMT.csv"); cat("\n== TABLE 5: HTMT ==\n"); print(H)

## ---- discriminant ANX vs THR ----
one<-'F=~ANX1+ANX2+ANX3+ANX4+ANX5+ANX6+ANX7+ANX8+ANX9+THR1+THR2+THR3'
two<-'ANX=~ANX1+ANX2+ANX3+ANX4+ANX5+ANX6+ANX7+ANX8+ANX9
      THR=~THR1+THR2+THR3'
f1<-cfa(one,dat,estimator=EST); f2<-cfa(two,dat,estimator=EST)
cat("\n== Discriminant ANX vs THR (1- vs 2-factor) ==\n"); print(rbind(OneFactor=sf(f1),TwoFactor=sf(f2)))
cat("\nSatorra-Bentler scaled difference:\n"); print(anova(f1,f2))
cat("Latent correlation ANX ~ THR:\n")
print(subset(standardizedSolution(f2),op=="~~"&lhs!=rhs)[,c("lhs","rhs","est.std","se","ci.lower","ci.upper")])

## ---- TABLE 7: composite correlations ----
comp<-data.frame(ANX=rowMeans(dat[,ANX]),THR=rowMeans(dat[,THR]),ASE=rowMeans(dat[,ASE]),CSE=rowMeans(dat[,CSE]))
B7<-round(cor(comp,use="pairwise"),3)
wcsv(as.data.frame(B7),"Table7_composite_correlations.csv"); cat("\n== TABLE 7: Composite correlations ==\n"); print(B7)
dat$ANX_mean<-comp$ANX

## ---- TABLE 6: measurement invariance ----
dat$grp_disc <- ifelse(dat$khoi_nganh=="Công nghệ thông tin","IT","Non-IT")
dat$grp_year <- ifelse(dat$nam_hoc %in% c("Năm 1","Năm 2"),"Junior","Senior")
dat$grp_gender <- ifelse(dat$gioi_tinh=="Nam","Male", ifelse(dat$gioi_tinh=="Nữ","Female",NA))
inv<-function(g,d){
  r<-rbind(configural=sf(cfa(m_uni,d,estimator=EST,group=g)),
           metric    =sf(cfa(m_uni,d,estimator=EST,group=g,group.equal="loadings")),
           scalar    =sf(cfa(m_uni,d,estimator=EST,group=g,group.equal=c("loadings","intercepts"))))
  cbind(r, dCFI=round(c(NA,r[1,"cfi.robust"]-r[2,"cfi.robust"],r[2,"cfi.robust"]-r[3,"cfi.robust"]),3))
}
B6<-NULL
labs<-c(grp_disc="Discipline (IT vs Non-IT)",grp_year="Year (Junior vs Senior)",grp_gender="Gender (Male vs Female)")
for(g in c("grp_disc","grp_year","grp_gender")){
  d<-if(g=="grp_gender") subset(dat,!is.na(grp_gender)) else dat
  tb<-try(inv(g,d),silent=TRUE)
  cat("\n== Invariance:",labs[[g]],"==\n"); print(tb)
  if(!inherits(tb,"try-error")) B6<-rbind(B6,data.frame(Grouping=labs[[g]],Step=rownames(tb),tb,row.names=NULL))
}
if(!is.null(B6)) wcsv(B6,"Table6_invariance.csv",rn=FALSE)


## ---- SPLIT-SAMPLE EFA / HOLDOUT CFA, COMMON-METHOD CHECK, CORRELATES ----
set.seed(2026); sidx <- sample(nrow(dat)); hA <- dat[sidx[1:floor(nrow(dat)/2)],]; hB <- dat[sidx[(floor(nrow(dat)/2)+1):nrow(dat)],]
kmo <- psych::KMO(hA[,ANX])$MSA; bt <- psych::cortest.bartlett(cor(hA[,ANX]), n=nrow(hA))
ev <- eigen(cor(hA[,ANX]))$values
set.seed(1); pa <- suppressMessages(psych::fa.parallel(hA[,ANX], fa="fa", fm="ml", plot=FALSE, n.iter=500))
efa1 <- psych::fa(hA[,ANX], nfactors=1, fm="ml")
cat("\n== EFA (random half A, n =", nrow(hA), ") ==\n")
cat("KMO =", round(kmo,3), "| Bartlett chi2(", bt$df, ") =", round(bt$chisq,2), "p =", signif(bt$p.value,3), "\n")
cat("Eigenvalues:", round(ev,2), "\nParallel analysis factors:", pa$nfact, "\n")
cat("One-factor loadings:", round(efa1$loadings[,1],2), "| variance explained:", round(efa1$Vaccounted[2,1],3), "\n")
efa3 <- suppressMessages(psych::fa(hA[,ANX], nfactors=3, fm="ml", rotate="oblimin"))
cat("Three-factor oblimin loadings:\n"); print(round(unclass(efa3$loadings),2)); print(round(efa3$Phi,2))
fhB <- cfa(m_uni, hB, estimator=EST)
cat("\n== Holdout CFA (random half B, n =", nrow(hB), ") ==\n")
print(round(fitMeasures(fhB, c(idx,"rmsea.ci.lower.robust","rmsea.ci.upper.robust")),3))
B8 <- data.frame(Item=ANX, EFA_halfA=round(as.numeric(efa1$loadings[,1]),2),
                 CFA_halfB=round(standardizedSolution(fhB)$est.std[standardizedSolution(fhB)$op=="=~"],2))
wcsv(B8, "Table8_split_sample_loadings.csv", rn=FALSE)
cat("\n== Single-factor test for common-method variance (18 items) ==\n")
all18 <- c(THR, ANX, ASE, CSE)
f18 <- cfa(paste("G=~", paste(all18, collapse="+")), dat, estimator=EST); print(sf(f18))
cat("\n== RMSEA 90% CI ==\n")
print(round(rbind(Unidimensional=fitMeasures(f_uni, c("rmsea.ci.lower.robust","rmsea.ci.upper.robust")),
                  Measurement=fitMeasures(f_meas, c("rmsea.ci.lower.robust","rmsea.ci.upper.robust"))),3))
cat("\n== Correlates of the anxiety composite ==\n")
use_lev <- c("Chưa từng","Hiếm khi","Thỉnh thoảng","Thường xuyên","Hằng ngày")
ct <- cor.test(dat$ANX_mean, match(dat$tan_suat_AI, use_lev), method="spearman", exact=FALSE)
cat("AI use frequency: Spearman rho =", round(ct$estimate,3), "p =", round(ct$p.value,3), "\n")
ow <- oneway.test(ANX_mean ~ da_hoc_AI, dat); av <- summary(aov(ANX_mean ~ da_hoc_AI, dat))[[1]]
cat("Studied AI: Welch F(", ow$parameter[1], ",", round(ow$parameter[2],1), ") =", round(ow$statistic,2), "p =", signif(ow$p.value,3),
    "eta2 =", round(av[1,"Sum Sq"]/sum(av[,"Sum Sq"]),3), "\n"); print(round(tapply(dat$ANX_mean, dat$da_hoc_AI, mean),2))
print(pairwise.t.test(dat$ANX_mean, dat$da_hoc_AI, pool.sd=FALSE, p.adjust.method="holm"))
yr4 <- dat$nam_hoc; yr4[yr4 %in% c("Năm 4","Trên 4 năm")] <- "Năm 4+"
print(pairwise.t.test(dat$ANX_mean, yr4, pool.sd=FALSE, p.adjust.method="holm"))
print(pairwise.t.test(dat$ANX_mean, dat$khoi_nganh, pool.sd=FALSE, p.adjust.method="holm"))

cat("\n===== DONE =====\n"); sink()
message(">>> DONE. Tables: ", DIR_BANG, " | Figures: ", DIR_HINH)
cat("\nTables:\n"); print(list.files(DIR_BANG))
cat("Figures:\n"); print(list.files(DIR_HINH))
