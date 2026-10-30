# ============ 07. RA SOAT TRUOC KHI CONG BO LEN GITHUB ============
# Chay sau 01-06. Kiem tra 4 nhom loi va DUNG chuong trinh neu co loi nghiem trong:
#   A. Du lieu: dung cot, khong co thong tin ca nhan, gia tri hop le, dung N.
#   B. Quy trinh lam sach: bang luong mau khop voi du lieu.
#   C. So lieu trong bai: tinh lai doc lap va so voi gia tri da bao cao.
#   D. Kho ma: khong lo email, so dien thoai, duong dan may ca nhan, file du lieu tho.
# Ket qua ghi vao ket_qua/Kiem_tra_truoc_khi_cong_bo.txt
suppressMessages({ library(lavaan); library(psych) })
dat <- read.csv(DATA_FILE, fileEncoding = "UTF-8-BOM")
ANX <- paste0("ANX", 1:9); THR <- paste0("THR", 1:3)
ITEMS <- c(THR, ANX, "ATTN", paste0("ASE", 1:3), paste0("CSE", 1:3))
log <- data.frame(Nhom = character(), Kiem_tra = character(), Ket_qua = character(), Chi_tiet = character())
chk <- function(nhom, ten, ok, chi_tiet = "", nghiem_trong = TRUE)
  log[nrow(log) + 1, ] <<- list(nhom, ten, if (ok) "DAT" else if (nghiem_trong) "LOI" else "CANH BAO", chi_tiet)

## ---- A. Du lieu ----
cot_dung <- c("ID", "khoi_nganh", "nam_hoc", "gioi_tinh", "tan_suat_AI", "da_hoc_AI", ITEMS)
chk("A", "Dung 25 cot theo bang hoi thang 10 (khong co cot them)", identical(names(dat), cot_dung),
    paste(setdiff(names(dat), cot_dung), collapse = ", "))
pii <- grep("mail|time|date|ngay|phone|sdt|ten|name|mssv|student|ip", names(dat), ignore.case = TRUE, value = TRUE)
chk("A", "Khong co cot dinh danh (email, thoi gian, ten, MSSV...)", length(pii) == 0, paste(pii, collapse = ", "))
txt <- unlist(lapply(dat[sapply(dat, is.character)], unique))
chk("A", "Khong co email/so dien thoai trong noi dung o", !any(grepl("@|\\b0[0-9]{9}\\b", txt)))
chk("A", "N = 464", nrow(dat) == 464, paste("N =", nrow(dat)))
chk("A", "ID lien tuc 1..N, khong trung", identical(dat$ID, seq_len(nrow(dat))))
chk("A", "Khong co o trong", sum(is.na(dat)) == 0, paste(sum(is.na(dat)), "o trong"))
chk("A", "Moi cau tra loi nam trong 1-5", all(as.matrix(dat[ITEMS]) %in% 1:5))
chk("A", "Tat ca qua cau chu y (ATTN = 2)", all(dat$ATTN == 2))
sl <- apply(dat[setdiff(ITEMS, "ATTN")], 1, function(r) length(unique(r)) == 1)
chk("A", "Khong con phieu straightlining", !any(sl), paste(sum(sl), "phieu"))
dup <- sum(duplicated(dat[-1]))
chk("A", "Khong co dong trung hoan toan", dup == 0, paste(dup, "dong"), nghiem_trong = FALSE)
muc <- list(khoi_nganh = c("Kinh doanh / Quản trị", "Công nghệ thông tin", "Thiết kế / Truyền thông", "Khác"),
            nam_hoc = c("Năm 1", "Năm 2", "Năm 3", "Năm 4", "Trên 4 năm"),
            gioi_tinh = c("Nam", "Nữ", "Khác", "Không muốn nêu"),
            tan_suat_AI = c("Chưa từng", "Hiếm khi", "Thỉnh thoảng", "Thường xuyên", "Hằng ngày"),
            da_hoc_AI = c("Chưa", "Có, không chính thức", "Có, chính thức trong chương trình"))
for (v in names(muc)) { la <- setdiff(unique(dat[[v]]), muc[[v]])
  chk("A", paste("Gia tri hop le:", v), length(la) == 0, paste(la, collapse = ", ")) }

## ---- B. Luong mau ----
fl <- read.csv(file.path(DIR_BANG, "Table0_participant_flow.csv"))
chk("B", "Bang luong mau: 583 phieu nhan", fl$Remaining[1] == 583, paste(fl$Remaining[1]))
chk("B", "Bang luong mau: so loai cong lai khop", fl$Remaining[1] - sum(fl$Removed, na.rm = TRUE) == nrow(dat),
    paste(fl$Remaining[1], "-", sum(fl$Removed, na.rm = TRUE), "=", fl$Remaining[1] - sum(fl$Removed, na.rm = TRUE)))
chk("B", "Nhom gioi tinh trong MIMIC = 442 (loai 22 'Khac'/'Khong muon neu')",
    sum(dat$gioi_tinh %in% c("Nam", "Nữ")) == 442)

## ---- C. So lieu bao cao trong bai (tinh lai doc lap) ----
near <- function(x, y, tol = .0015) abs(x - y) <= tol
a <- psych::alpha(dat[ANX], warnings = FALSE)$total$raw_alpha
chk("C", "Alpha ANX = .878", near(a, .878), sprintf("%.4f", a))
set.seed(2026); s <- sample(nrow(dat)); hB <- dat[s[(floor(nrow(dat) / 2) + 1):nrow(dat)], ]
fB <- cfa(paste("ANX =~", paste(ANX, collapse = " + ")), hB, estimator = "MLR")
fm <- fitMeasures(fB, c("cfi.robust", "rmsea.robust"))
chk("C", "CFA holdout: CFI = .975", near(fm[1], .975), sprintf("%.4f", fm[1]))
chk("C", "CFA holdout: RMSEA = .057", near(fm[2], .057), sprintf("%.4f", fm[2]))
f2 <- cfa(paste("ANX =~", paste(ANX, collapse = " + "), "\nTHR =~", paste(THR, collapse = " + ")), dat, estimator = "MLR")
r <- lavInspect(f2, "std")$psi["ANX", "THR"]
chk("C", "Tuong quan tiem an ANX-THR = .87", near(r, .87, .005), sprintf("%.3f", r))
cl <- tryCatch(jsonlite::fromJSON(file.path(OUT_TXT, "so_lieu_bai_bao.json")), error = function(e) NULL)
if (!is.null(cl)) {
  sl_irt <- range(unlist(cl$chat_luong$irt$a_logit))
  chk("C", "Do doc IRT 1.47-1.84", near(sl_irt[1], 1.47, .005) && near(sl_irt[2], 1.84, .005), paste(sl_irt, collapse = "-"))
  ir <- cl$chat_luong$irt
  chk("C", "Vi tri IRT -2.59 den 1.51 (Table S7)", near(ir$b_min, -2.59, .005) && near(ir$b_max, 1.51, .005), paste(ir$b_min, ir$b_max, sep = " den "))
  chk("C", "Do tin cay co dieu kien >= .80 tu theta -3.3 den 2.1", near(ir$rel80_lo, -3.3, .001) && near(ir$rel80_hi, 2.1, .001), paste(ir$rel80_lo, ir$rel80_hi, sep = " den "))
  chk("C", "Ty le mau trong khoang do (99.1%)", near(ir$share_in_range, .991, .0005), sprintf("%.3f", ir$share_in_range))
  L <- cl$phan_bien$esem$loadings
  x14 <- abs(L$THR_f[L$item %in% paste0("ANX", 1:4)]); x59 <- abs(L$THR_f[L$item %in% paste0("ANX", 5:9)])
  chk("C", "ESEM: ANX1-4 tai cheo .24-.28 len THR", min(x14) >= .24 && max(x14) <= .28, paste(range(x14), collapse = "-"))
  chk("C", "ESEM: ANX5-9 |tai cheo| <= .11", max(x59) <= .11, sprintf("max = %.2f", max(x59)))
} else chk("C", "Doc so_lieu_bai_bao.json", FALSE, "chua chay 05/06?")

## ---- D. Kho ma ----
files <- list.files(".", recursive = TRUE, all.files = TRUE, full.names = TRUE)
files <- files[!grepl("^\\./\\.git/", files)]
tho <- grep("\\.xlsx?$|Tong_hop|raw|Form Responses|\\.Rhistory|\\.DS_Store", files, ignore.case = TRUE, value = TRUE)
chk("D", "Khong co file du lieu tho/rac trong kho", length(tho) == 0, paste(tho, collapse = ", "))
txtf <- files[grepl("\\.(R|md|txt|csv|json|cff)$", files, ignore.case = TRUE) & !grepl("07_Kiem_tra", files)]  # bo qua chinh file nay (chua mau regex)
bad <- character()
for (f in txtf) { z <- readLines(f, warn = FALSE, encoding = "UTF-8")
  hit <- grepl("[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[a-z]{2,}|/Users/|C:\\\\Users|\\b0[0-9]{9}\\b", z)
  if (any(hit)) bad <- c(bad, paste0(f, ":", which(hit)[1])) }
chk("D", "Khong lo email/so dien thoai/duong dan may ca nhan", length(bad) == 0, paste(bad, collapse = "; "))
big <- files[file.size(files) > 50e6]
chk("D", "Khong co file > 50 MB", length(big) == 0, paste(big, collapse = ", "))

## ---- Bao cao ----
out <- file.path(OUT_TXT, "Kiem_tra_truoc_khi_cong_bo.txt")
op <- options(width = 160); sink(out, split = TRUE)
cat("RA SOAT TRUOC KHI CONG BO -", format(Sys.time(), "%Y-%m-%d %H:%M"), "\nDu lieu:", DATA_FILE, "\n\n")
print(log, row.names = FALSE, right = FALSE)
n_loi <- sum(log$Ket_qua == "LOI")
cat("\nTong:", nrow(log), "kiem tra |", sum(log$Ket_qua == "DAT"), "dat |", n_loi, "loi |",
    sum(log$Ket_qua == "CANH BAO"), "canh bao\n")
sink(); options(op)
if (n_loi > 0) stop("Co ", n_loi, " loi. Xem ", out, " truoc khi day len GitHub.")
message(">>> 07 xong: khong co loi nghiem trong.")
