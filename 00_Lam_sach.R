# ================== 00. LAM SACH DU LIEU TU FILE EXCEL TONG HOP ==================
options(repos = c(CRAN = "https://cloud.r-project.org"))
if(!requireNamespace("readxl", quietly=TRUE))
  tryCatch(install.packages("readxl", type="binary"), error=function(e) install.packages("readxl"))

codes <- c("THR1","THR2","THR3", paste0("ANX",1:9), "ATTN", "ASE1","ASE2","ASE3","CSE1","CSE2","CSE3")
raw <- as.data.frame(readxl::read_excel(RAW_FILE, sheet=RAW_SHEET))
stopifnot(ncol(raw) >= 26)
raw <- raw[, 1:26]
names(raw) <- c("Timestamp","Email","khoi_nganh","nam_hoc","gioi_tinh","tan_suat_AI","da_hoc_AI", codes)
raw[raw == "#N/A"] <- NA                                            # o loi tu cong thuc Excel
raw <- raw[!(is.na(raw$Email) & rowSums(!is.na(raw[, codes])) == 0), ] # bo dong trong
for(cc in codes) raw[[cc]] <- suppressWarnings(as.numeric(raw[[cc]]))
raw$Email <- tolower(trimws(as.character(raw$Email)))
to_time <- function(x){
  if(inherits(x, "POSIXct")) return(x)
  x <- as.character(x); num <- suppressWarnings(as.numeric(x))
  out <- as.POSIXct(num*86400, origin="1899-12-30", tz="UTC")       # so serial Excel
  x[x %in% c("#N/A","N/A","NA","")] <- NA
  miss <- is.na(num) & !is.na(x)
  if(any(miss)) out[miss] <- as.POSIXct(x[miss], tz="UTC", tryFormats=c("%Y-%m-%d %H:%M:%OS","%m/%d/%Y %H:%M:%S","%d/%m/%Y %H:%M:%S"))
  out
}
raw$Timestamp <- to_time(raw$Timestamp)

flow <- data.frame(Step=character(), Removed=integer(), Remaining=integer())
add <- function(step, before, after) flow[nrow(flow)+1, ] <<- list(step, before-after, after)
n0 <- nrow(raw); flow[1, ] <- list("Responses received", NA, n0)

d <- raw[order(raw$Timestamp), ]
d <- d[!duplicated(d$Email, fromLast=TRUE), ]                     # giu lan nop cuoi cung
add("Duplicate email (kept latest submission)", n0, nrow(d)); n1 <- nrow(d)
d <- d[complete.cases(d[, codes]), ]
add("Incomplete item responses", n1, nrow(d)); n2 <- nrow(d)
d <- d[d$ATTN == 2, ]
add("Failed attention check (item 13, correct = 2)", n2, nrow(d)); n3 <- nrow(d)
subst <- setdiff(codes, "ATTN")
d$straightline <- apply(d[, subst], 1, function(r) length(unique(r)) == 1)
if(EXCLUDE_STRAIGHTLINE){
  d <- d[!d$straightline, ]
  add("Straightlining (identical answer on all 18 items)", n3, nrow(d))
}
flow$Removed[1] <- NA

out <- d[, c("khoi_nganh","nam_hoc","gioi_tinh","tan_suat_AI","da_hoc_AI", codes)]
out <- cbind(ID=seq_len(nrow(out)), out)
write.csv(out, DATA_FILE, row.names=FALSE, fileEncoding="UTF-8")
dir.create(DIR_BANG, recursive=TRUE, showWarnings=FALSE)
write.csv(flow, file.path(DIR_BANG, "Table0_participant_flow.csv"), row.names=FALSE, fileEncoding="UTF-8")
message(">>> Lam sach xong: ", n0, " phieu tho -> ", nrow(out), " phieu dung de phan tich")
print(flow, row.names=FALSE)
