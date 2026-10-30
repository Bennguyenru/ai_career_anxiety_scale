# ===================================================================
#  CHAY TAT CA: mo file nay trong RStudio roi bam nut Source (goc tren ben phai)
# ===================================================================
get_this_dir <- function(){
  if (file.exists("00_CONFIG.R") && file.exists("01_Phan_tich.R")) return(getwd())   # da o dung thu muc
  of <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
  if (!is.null(of) && file.exists(file.path(dirname(of), "00_CONFIG.R"))) return(dirname(normalizePath(of)))
  p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error=function(e) "")
  if (!is.null(p) && nzchar(p)) return(dirname(p))
  a <- commandArgs(FALSE); f <- sub("^--file=","", a[grep("^--file=",a)])
  if (length(f)) return(dirname(normalizePath(f)))
  of <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
  if (!is.null(of)) return(dirname(normalizePath(of)))
  getwd()
}
setwd(get_this_dir())
message(">>> Thu muc lam viec: ", getwd())
source("00_CONFIG.R")
if (file.exists(RAW_FILE)) source("00_Lam_sach.R") else message(">>> Khong thay RAW_FILE, dung truc tiep ", DATA_FILE)
source("01_Phan_tich.R")
source("02_Figures.R")
source("03_Figures_journal.R")
source("04_Phan_tich_bo_sung.R")
source("05_Phan_tich_chat_luong.R")
source("06_Phan_tich_phan_bien.R")
source("07_Kiem_tra_truoc_khi_cong_bo.R")   # ra soat truoc khi day len GitHub
