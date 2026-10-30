# ============================ CAU HINH ============================
# QUY TRINH: dat file Excel tong hop (xuat tu Google Form) vao data/ roi chay CHAY_TAT_CA.R.
#   Buoc 1 (00_Lam_sach.R): doc RAW_FILE, khu trung email, bo phieu thieu, loc cau chu y,
#           loai phieu tra loi mot muc cho moi cau (straightlining) -> ghi DATA_FILE.
#   Buoc 2 (01_Phan_tich.R): phan tich + xuat bang.
#   Buoc 3 (02_Figures.R / 03_Figures_journal.R): ve hinh.
# Neu khong co RAW_FILE, chuong trinh dung thang DATA_FILE (ban da lam sach san).
RAW_FILE   <- "data/Tong_hop_form.xlsx"   # file Excel tong hop
RAW_SHEET  <- "Raw Data"                 # ten sheet chua du lieu tho
EXCLUDE_STRAIGHTLINE <- TRUE             # TRUE = loai phieu chon cung mot muc cho ca 18 cau
DATA_FILE  <- "data/ai_career_anxiety_oct2026_n464.csv"   # 464 phieu hop le, thu ngay 27-30/10/2026 (sau phe duyet dao duc 27/10/2026); file nay KHONG dua len GitHub
OUT_TXT    <- "ket_qua"
DIR_BANG   <- "ket_qua/bang"
DIR_HINH   <- "ket_qua/hinh"
DIR_HINH_J <- "ket_qua/hinh_journal"
JOURNAL_WIDTH_MM <- 165                  # do rong hinh toan trang (mm): 165 vua kho Letter va A4
