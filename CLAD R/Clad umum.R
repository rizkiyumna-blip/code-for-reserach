# =====================================================================
# PENGATURAN  
# =====================================================================
file_data     <- "book2.xlsx"                 # nama file data Excel (satu folder dengan getwd())
nama_y        <- "Y"                          # nama kolom variabel dependen
nama_ara      <- "ARA"                        # nama kolom batas ARA per baris
var_x         <- c("X1", "X2", "Z", "X2_Z")   # semua variabel independen
# contoh 1 variabel : c("X1")
# contoh 2 variabel : c("X1", "X2")
var_strata    <- "X1"                         # dummy untuk bootstrap berstrata; NULL = tanpa strata
var_interaksi <- "X2_Z"                       # variabel interaksi (hanya untuk catatan VIF); NULL = tidak ada
# Catatan: tulis nama kolom persis seperti di Excel (huruf besar/kecil harus sama),

# =====================================================================

if (!require(readxl)) install.packages("readxl")
library(readxl)
data_skripsi <- read_excel(file_data)

# Pemeriksaan awal: nama kolom & sel kosong
nama_perlu <- unique(c(nama_y, nama_ara, var_x, var_strata))
hilang <- setdiff(nama_perlu, names(data_skripsi))
if (length(hilang) > 0)
  stop("Kolom ini tidak ditemukan di data: ", paste(hilang, collapse = ", "),
       ". Cek penulisan nama kolom (huruf besar/kecil harus sama).")
ada_na <- names(which(colSums(is.na(as.data.frame(data_skripsi)[nama_perlu])) > 0))
if (length(ada_na) > 0)
  stop("Ada sel kosong (NA) pada kolom: ", paste(ada_na, collapse = ", "),
       ". Lengkapi atau hapus baris tersebut di Excel.")
# =====================================================================
# CLAD dua sisi (kiri = 0, kanan = batas ARA per baris)
# Jalankan dari atas ke bawah. Syarat: objek `data_skripsi` ada di Environment.
# =====================================================================

# ---- 0. Persiapan --------------------------------------------------
if (!require(quantreg)) install.packages("quantreg")
library(quantreg)

# Bersihkan objek lama (data dan pengaturan dipertahankan)
rm(list = setdiff(ls(), c("data_skripsi", "file_data", "nama_y", "nama_ara",
                          "var_x", "var_strata", "var_interaksi")))

# ---- 1. Siapkan data -----------------------------------------------
# Kolom ARA berisi batas ARA nominal tiap baris (mis. 0.35 / 0.25)
d <- as.data.frame(data_skripsi)
d$Y <- d[[nama_y]]                                 # kolom Y dibakukan namanya di dalam skrip
toleransi <- 0.02                                  # toleransi pembulatan tick size

d$batas      <- d[[nama_ara]]
d$sentuh_ARA <- d$Y >= d$batas - toleransi         # TRUE = kena ARA
d$U          <- ifelse(d$sentuh_ARA, d$Y, d$batas) # batas kanan per observasi

# Ringkasan ARA
tabel_ara <- data.frame(
  Status = c("Tidak kena ARA", "Kena ARA"),
  Jumlah = as.vector(table(factor(d$sentuh_ARA, levels = c(FALSE, TRUE)))))
tabel_ara$Persen <- round(100 * tabel_ara$Jumlah / nrow(d), 2)
print(tabel_ara)
cat("Y melebihi batas nominal (harus 0):", sum(d$Y > d$batas), "\n")
cat("Y tepat bernilai 0 (sensor kiri)  :", sum(d$Y == 0), "\n\n")

# ---- 2. Matriks desain ---------------------------------------------
y <- d$Y
U <- d$U
X <- cbind(Intersep = 1, as.matrix(d[, var_x, drop = FALSE]))   # otomatis mengikuti var_x
n <- nrow(X); k <- ncol(X)

# Standardisasi prediktor (hanya untuk stabilitas optimasi)
mu  <- colMeans(X[, -1, drop = FALSE])
sdv <- apply(X[, -1, drop = FALSE], 2, sd)
Xs  <- X
Xs[, -1] <- scale(X[, -1, drop = FALSE], center = mu, scale = sdv)

# Fungsi untuk mengembalikan koefisien ke skala asli
ke_asli <- function(b) {
  out <- b
  out[-1] <- b[-1] / sdv
  out[1]  <- b[1] - sum(b[-1] * mu / sdv)
  out
}

# ---- 3. Fungsi CLAD ------------------------------------------------
# Fungsi objektif Powell dua sisi
obj_clad <- function(b, y, X, U) {
  xb   <- as.vector(X %*% b)
  pred <- pmin(U, pmax(0, xb))
  sum(abs(y - pred))
}

# Iterasi Powell: regresi median pada baris yang tidak tersensor
clad_iter <- function(y, X, U, b0, maxit = 200, tol = 1e-8) {
  b <- b0
  for (i in 1:maxit) {
    xb  <- as.vector(X %*% b)
    idx <- xb > 0 & xb < U
    if (sum(idx) <= ncol(X) + 2) break
    Xi <- X[idx, , drop = FALSE]
    if (qr(Xi)$rank < ncol(X)) break
    b_new <- tryCatch(
      suppressWarnings(coef(rq.fit.br(Xi, y[idx], tau = 0.5))),
      error = function(e) NULL)
    if (is.null(b_new) || any(!is.finite(b_new))) break
    if (max(abs(b_new - b)) < tol) { b <- b_new; break }
    b <- b_new
  }
  b
}

# Beberapa titik awal, ambil yang objektifnya paling kecil
fit_clad <- function(y, X, U, n_start = 5) {
  b_ols  <- coef(lm.fit(X, y))
  starts <- c(list(b_ols),
              lapply(seq_len(n_start), function(i) b_ols + rnorm(ncol(X), 0, 0.05)))
  hasil <- lapply(starts, function(s) {
    b <- clad_iter(y, X, U, s)
    list(par = b, value = obj_clad(b, y, X, U))
  })
  hasil[[which.min(sapply(hasil, `[[`, "value"))]]
}

# ---- 4. Estimasi utama ---------------------------------------------
set.seed(123)
fit      <- fit_clad(y, Xs, U, n_start = 10)
beta_std <- fit$par                      # skala standar
beta_hat <- ke_asli(fit$par)             # skala asli
cat("Nilai objektif minimum:", round(fit$value, 5), "\n")
cat("Koefisien (skala asli):\n"); print(round(beta_hat, 5))

# ---- 5. Bootstrap (berstrata jika var_strata diisi) ----------------
set.seed(123)
R <- 1000                                # uji coba: 100; hasil akhir: 1000
boot_std <- matrix(NA, R, k, dimnames = list(NULL, colnames(X)))
pesan <- character(0)
grup <- if (!is.null(var_strata)) split(seq_len(n), d[[var_strata]]) else list(seq_len(n))

for (r in 1:R) {
  idx <- unlist(lapply(grup, function(g) g[sample.int(length(g), length(g), replace = TRUE)]),
                use.names = FALSE)
  hasil_r <- tryCatch(fit_clad(y[idx], Xs[idx, ], U[idx], n_start = 2)$par,
                      error = function(e) { pesan <<- c(pesan, conditionMessage(e)); NULL })
  if (!is.null(hasil_r)) boot_std[r, ] <- hasil_r
}

ok <- complete.cases(boot_std)
cat("\nSampel bootstrap gagal:", sum(!ok), "dari", R,
    "(", round(100 * sum(!ok) / R, 2), "% )\n")
if (length(pesan) > 0) print(table(pesan))

# Kembalikan hasil bootstrap ke skala asli
boot_asli <- t(apply(boot_std, 1, ke_asli))
colnames(boot_asli) <- colnames(X)
se <- apply(boot_asli[ok, ], 2, sd)
ci <- t(apply(boot_asli[ok, ], 2, quantile, probs = c(0.025, 0.975)))


# ---- 6. Uji multikolinearitas (VIF & Tolerance) --------------------
Xp  <- X[, -1, drop = FALSE]
vif <- if (ncol(Xp) < 2) rep(1, ncol(Xp)) else sapply(seq_len(ncol(Xp)), function(j) {
  r2 <- summary(lm(Xp[, j] ~ Xp[, -j]))$r.squared
  1 / (1 - r2)
})
tabel_vif <- data.frame(VIF = vif, Tolerance = 1 / vif,
                        row.names = colnames(Xp))
cat("\n=== Uji Multikolinearitas ===\n")
print(round(tabel_vif, 4))

# ---- 7. Koefisien & Bootstrapped Standard Error --------------------
tabel_koef <- data.frame(Estimasi = beta_hat, Boot_SE = se,
                         CI_2.5 = ci[, 1], CI_97.5 = ci[, 2])
cat("\n=== Koefisien & Bootstrapped SE (skala asli) ===\n")
print(round(tabel_koef, 5))

# ---- 8. Uji signifikansi parsial (uji t) ---------------------------
t_hit <- beta_hat / se
p_t   <- 2 * pt(-abs(t_hit), df = n - k)
tabel_t <- data.frame(
  Estimasi = beta_hat, Boot_SE = se, t_hitung = t_hit, p_value = p_t,
  Sig = as.character(cut(p_t, c(-Inf, .01, .05, .10, Inf),
                         labels = c("***", "**", "*", "ns"))))
cat("\n=== Uji t (parsial) ===\n")
print(tabel_t, digits = 5)
cat("Keterangan: *** p<0.01, ** p<0.05, * p<0.10, ns = tidak signifikan\n")

# ---- 9. Uji signifikansi simultan (Wald / Chi-Square) --------------
b_slope <- beta_hat[-1]
V_boot  <- cov(boot_asli[ok, , drop = FALSE])[-1, -1, drop = FALSE]
wald    <- as.numeric(t(b_slope) %*% solve(V_boot) %*% b_slope)
df_wald <- length(b_slope)
p_wald  <- pchisq(wald, df = df_wald, lower.tail = FALSE)
cat("\n=== Uji Wald (Chi-Square) ===\n",
    "Chi-square =", round(wald, 4), "\n",
    "df         =", df_wald, "\n",
    "p-value    =", format.pval(p_wald, digits = 4), "\n")

# ---- 10. Pseudo R-Square -------------------------------------------
# Model nol (intersep saja): cari konstanta terbaik lewat grid
grid <- seq(0, max(U), length.out = 5000)
V0   <- min(sapply(grid, function(b) sum(abs(y - pmin(U, pmax(0, b))))))
V1   <- fit$value
pseudo_R2 <- 1 - V1 / V0

pred <- pmin(U, pmax(0, as.vector(X %*% beta_hat)))
R2_korelasi <- cor(y, pred)^2

cat("\n=== Pseudo R-Square ===\n",
    "Pseudo R2 (1 - V1/V0)           =", round(pseudo_R2, 4), "\n",
    "R2 korelasi (aktual vs prediksi) =", round(R2_korelasi, 4), "\n")

# ============================================================
# ---- 11. LAPORAN OUTPUT ------------------------------------
# ============================================================

nm_koef <- c("C", var_x)   # label gaya (C = konstanta)
lebar   <- 84                                # lebar garis laporan

# ---- 11.1 Fungsi bantu format -----------------------------------------
gr  <- function(ch = "-") cat(strrep(ch, lebar), "\n", sep = "")
jdl <- function(teks) { cat("\n"); gr("="); cat(teks, "\n", sep = ""); gr("=") }
f   <- function(x, dg = 6, w = 12) formatC(x, format = "f", digits = dg, width = w)
fs  <- function(s, w = 12) formatC(s, width = w)                # rata kanan
fl  <- function(s, w = 12) formatC(s, width = w, flag = "-")    # rata kiri
bintang <- function(p) ifelse(p < .01, "***", ifelse(p < .05, "**",
                                                     ifelse(p < .10, "*", "ns")))

# Cetak matriks (baris = nama, kolom = variabel) ala tabel EViews
cetak_matriks <- function(m, dg = 6, w = 13, lw = 16) {
  cat(fl("", lw), paste(fs(colnames(m), w), collapse = ""), "\n", sep = "")
  gr()
  for (i in seq_len(nrow(m))) {
    dgi <- if (rownames(m)[i] == "Observations") 0 else dg
    cat(fl(rownames(m)[i], lw), paste(f(m[i, ], dgi, w), collapse = ""),
        "\n", sep = "")
  }
  gr()
}

# ---- 11.2 Hitung tabel yang dibutuhkan --------------------------------
# (a) Statistik deskriptif (definisi momen sama seperti EViews)
vars_desk <- data.frame(d$Y, as.data.frame(X[, -1, drop = FALSE]), check.names = FALSE)
names(vars_desk)[1] <- nama_y

desk <- function(v) {
  m   <- mean(v); dv <- v - m; s2 <- mean(dv^2)
  S   <- mean(dv^3) / s2^1.5
  K   <- mean(dv^4) / s2^2
  JB  <- length(v) / 6 * (S^2 + (K - 3)^2 / 4)
  c("Mean" = m, "Median" = median(v), "Maximum" = max(v), "Minimum" = min(v),
    "Std. Dev." = sd(v), "Skewness" = S, "Kurtosis" = K,
    "Jarque-Bera" = JB, "Probability" = pchisq(JB, df = 2, lower.tail = FALSE),
    "Sum" = sum(v), "Sum Sq. Dev." = sum(dv^2), "Observations" = length(v))
}
tabel_desk <- sapply(vars_desk, desk)

# (b) Korelasi antar prediktor
tabel_kor <- cor(vars_desk[, -1, drop = FALSE])

# (c) Tabel hasil estimasi (angka, untuk dicetak & diekspor)
p_F   <- pf(wald / df_wald, df_wald, n - k, lower.tail = FALSE)
adj_R2 <- 1 - (1 - pseudo_R2) * (n - 1) / (n - k)

tabel_hasil <- data.frame(
  Variable    = nm_koef,
  Coefficient = as.numeric(beta_hat),
  Std_Error   = as.numeric(se),
  t_Statistic = as.numeric(t_hit),
  Prob        = as.numeric(p_t),
  Sig         = bintang(as.numeric(p_t)),
  CI_Lower    = as.numeric(ci[, 1]),
  CI_Upper    = as.numeric(ci[, 2]))

n_kiri  <- sum(y == 0)
n_kanan <- sum(d$sentuh_ARA)
n_tak   <- n - n_kiri - n_kanan
ket_boot <- if (!is.null(var_strata)) paste0("Stratified bootstrap (strata ", var_strata, ")") else "Bootstrap biasa (tanpa strata)"

# ---- 11.3 Fungsi pencetak laporan -------------------------------------
laporan_clad <- function() {
  
  # ---------- A. Informasi umum ----------
  jdl("ESTIMASI CENSORED LEAST ABSOLUTE DEVIATIONS (CLAD)")
  cat("Dependent Variable : ", nama_y, "\n", sep = "")
  cat("Method             : Censored LAD - Powell (sensor dua sisi)\n")
  cat("Date / Time        : ", format(Sys.time(), "%m/%d/%Y   %H:%M"), "\n", sep = "")
  cat("Sample             : 1 ", n, "\n", sep = "")
  cat("Included obs       : ", n, "\n", sep = "")
  cat("Left-censored      : ", n_kiri, " obs (Y = 0)\n", sep = "")
  cat("Right-censored     : ", n_kanan, " obs (Y >= batas ARA), ",
      round(100 * n_kanan / n, 2), "%\n", sep = "")
  cat("Uncensored         : ", n_tak, " obs\n", sep = "")
  cat("Std. errors        : ", ket_boot, ", ", R,
      " replikasi, valid = ", sum(ok), "\n", sep = "")
  gr()
  
  # ---------- B. Statistik deskriptif ----------
  jdl("TABEL 1. STATISTIK DESKRIPTIF")
  cetak_matriks(tabel_desk, dg = 6, w = 13, lw = 16)
  
  # ---------- C. Multikolinearitas ----------
  jdl("TABEL 2. UJI MULTIKOLINEARITAS")
  cat(fl("Variable", 14), fs("VIF", 12), fs("Tolerance", 12),
      fs("Keputusan", 24), "\n", sep = "")
  gr()
  for (i in seq_len(nrow(tabel_vif))) {
    ket <- if (tabel_vif$VIF[i] < 10) "Bebas multikolinearitas" else "Indikasi multikolinearitas"
    cat(fl(rownames(tabel_vif)[i], 14), f(tabel_vif$VIF[i], 4, 12),
        f(tabel_vif$Tolerance[i], 4, 12), fs(ket, 26), "\n", sep = "")
  }
  gr()
  cat("Kriteria: VIF < 10 dan Tolerance > 0.10 -> tidak ada multikolinearitas.\n")
  if (any(tabel_vif$VIF >= 10) && length(var_interaksi) > 0)
    cat("Catatan : VIF tinggi pada variabel interaksi (", paste(var_interaksi, collapse = ", "),
        ") dengan\n", "          komponen pembentuknya lazim terjadi (multikolinearitas non-esensial).\n", sep = "")
  if (ncol(Xp) < 2)
    cat("Catatan : hanya ada satu variabel independen, sehingga VIF = 1 dan Tolerance = 1.\n")
  
  if (ncol(tabel_kor) >= 2) {
    cat("\nMatriks Korelasi Antar Variabel\n")
    gr()
    cetak_matriks(tabel_kor, dg = 4, w = 12, lw = 10)
  }
  
  # ---------- D. Hasil estimasi & uji t ----------
  jdl("TABEL 3. HASIL ESTIMASI CLAD DAN UJI SIGNIFIKANSI PARSIAL (UJI t)")
  cat(fl("Variable", 14), fs("Coefficient", 14), fs("Std. Error", 14),
      fs("t-Statistic", 14), fs("Prob.", 10), fs("Sig.", 6), "\n", sep = "")
  gr()
  for (i in seq_len(k)) {
    cat(fl(nm_koef[i], 14), f(beta_hat[i], 6, 14), f(se[i], 6, 14),
        f(t_hit[i], 4, 14), f(p_t[i], 4, 10), fs(bintang(p_t[i]), 6),
        "\n", sep = "")
  }
  gr()
  
  # Ringkasan model (dua kolom, gaya EViews)
  baris2 <- function(l1, v1, l2, v2, dg = 6) {
    cat(fl(l1, 30), f(v1, dg, 12), "    ", fl(l2, 26), f(v2, dg, 12),
        "\n", sep = "")
  }
  baris2("Pseudo R-squared",           pseudo_R2,   "Mean dependent var",   mean(y))
  baris2("Adjusted pseudo R-squared",  adj_R2,      "S.D. dependent var",   sd(y))
  baris2("R-squared (korelasi)",       R2_korelasi, "Objective (V1)",       V1)
  baris2("Observations",               n,           "Restr. objective (V0)", V0)
  gr()
  cat("Keterangan: *** p<0.01, ** p<0.05, * p<0.10, ns = tidak signifikan.\n")
  cat("Prob. dihitung dari distribusi t dengan df = n - k = ", n - k, ".\n", sep = "")
  cat("Pseudo R-squared = 1 - V1/V0 ; V1 = objektif model, V0 = objektif model konstanta.\n")
  
  # Interval kepercayaan bootstrap
  cat("\nInterval Kepercayaan 95% (Bootstrap Percentile)\n")
  cat(fl("Variable", 14), fs("Coefficient", 14), fs("Lower 2.5%", 14),
      fs("Upper 97.5%", 14), fs("Mencakup 0?", 14), "\n", sep = "")
  gr()
  for (i in seq_len(k)) {
    nol <- if (ci[i, 1] <= 0 && ci[i, 2] >= 0) "Ya" else "Tidak"
    cat(fl(nm_koef[i], 14), f(beta_hat[i], 6, 14), f(ci[i, 1], 6, 14),
        f(ci[i, 2], 6, 14), fs(nol, 14), "\n", sep = "")
  }
  gr()
  
  # ---------- E. Uji Wald (simultan) ----------
  jdl("TABEL 4. UJI SIGNIFIKANSI SIMULTAN (WALD TEST)")
  cat("Equation: CLAD (matriks kovarians bootstrap)\n\n")
  cat(fl("Test Statistic", 20), fs("Value", 14), fs("df", 14),
      fs("Probability", 14), "\n", sep = "")
  gr()
  cat(fl("F-statistic", 20), f(wald / df_wald, 6, 14),
      fs(paste0("(", df_wald, ", ", n - k, ")"), 14), f(p_F, 4, 14), "\n", sep = "")
  cat(fl("Chi-square", 20), f(wald, 6, 14), fs(df_wald, 14),
      f(p_wald, 4, 14), "\n", sep = "")
  gr()
  cat("Null Hypothesis : ", paste(nm_koef[-1], collapse = " = "), " = 0\n", sep = "")
  cat("Keputusan (alpha 5%) : ",
      if (p_wald < 0.05)
        "H0 ditolak -> variabel independen berpengaruh signifikan secara simultan.\n"
      else
        "H0 diterima -> variabel independen tidak berpengaruh signifikan secara simultan.\n",
      sep = "")
  gr("=")
  cat("Sumber: hasil olahan data (CLAD, R). Std. error = bootstrap berstrata.\n")
}

# ---- 11.4 Jalankan & simpan -------------------------------------------
laporan_clad()                                                  # tampil di Console

# Simpan seluruh laporan ke file teks (buka dengan Notepad / tempel ke Word,
# gunakan font Courier New / Consolas agar kolom sejajar)
capture.output(laporan_clad(), file = "Output_CLAD_EViews.txt")

# Ekspor SEMUA tabel ke SATU file Excel (tampilan rapi seperti laporan .txt)
# ---- 11.5 Ekspor Excel ----------------------------------------------------
if (!require(openxlsx)) install.packages("openxlsx")
library(openxlsx)

wb    <- createWorkbook()
ws    <- "Laporan_CLAD"
addWorksheet(wb, ws, gridLines = FALSE)
kol_n <- max(8, ncol(tabel_desk) + 1)      # jumlah kolom yang dipakai laporan
baris <- 1                                 # penunjuk baris yang sedang ditulis

# ---- Gaya tampilan ----
fn           <- "Consolas"
st_judul     <- createStyle(fontName = fn, fontSize = 11, textDecoration = "bold",
                            fontColour = "#FFFFFF", fgFill = "#1F3864")
st_teks      <- createStyle(fontName = fn, fontSize = 10)
st_catat     <- createStyle(fontName = fn, fontSize = 9, textDecoration = "italic")
st_sub       <- createStyle(fontName = fn, fontSize = 10, textDecoration = "bold")
st_head      <- createStyle(fontName = fn, fontSize = 10, textDecoration = "bold",
                            border = c("top", "bottom"), borderStyle = "thin",
                            halign = "right")
st_head_kiri <- createStyle(fontName = fn, fontSize = 10, textDecoration = "bold",
                            border = c("top", "bottom"), borderStyle = "thin",
                            halign = "left")
st_bawah     <- createStyle(border = "bottom", borderStyle = "thin")

# ---- Fungsi bantu ----
judul <- function(teks) {
  if (baris > 1) baris <<- baris + 1
  writeData(wb, ws, teks, startCol = 1, startRow = baris)
  addStyle(wb, ws, st_judul, rows = baris, cols = 1:kol_n, gridExpand = TRUE)
  baris <<- baris + 2
}

sub_judul <- function(teks) {
  writeData(wb, ws, teks, startCol = 1, startRow = baris)
  addStyle(wb, ws, st_sub, rows = baris, cols = 1)
  baris <<- baris + 1
}

info <- function(label, nilai) {
  writeData(wb, ws, label, startCol = 1, startRow = baris)
  writeData(wb, ws, paste0(": ", nilai), startCol = 2, startRow = baris)
  addStyle(wb, ws, st_teks, rows = baris, cols = 1:2, gridExpand = TRUE)
  baris <<- baris + 1
}

catat <- function(teks) {
  for (t in teks) {
    writeData(wb, ws, t, startCol = 1, startRow = baris)
    addStyle(wb, ws, st_catat, rows = baris, cols = 1)
    baris <<- baris + 1
  }
}

# fmt: vektor per kolom ("TEKS", "TEKS_KANAN", atau format angka) atau matriks per sel
tulis_tabel <- function(df, fmt) {
  nr <- nrow(df); nc <- ncol(df)
  if (!is.matrix(fmt)) fmt <- matrix(rep(fmt, each = nr), nrow = nr)
  writeData(wb, ws, df, startCol = 1, startRow = baris,
            colNames = TRUE, headerStyle = st_head)
  addStyle(wb, ws, st_head_kiri, rows = baris, cols = 1)
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    fm  <- fmt[i, j]
    stl <- if (fm == "TEKS") {
      createStyle(fontName = fn, fontSize = 10, halign = "left")
    } else if (fm == "TEKS_KANAN") {
      createStyle(fontName = fn, fontSize = 10, halign = "right")
    } else {
      createStyle(fontName = fn, fontSize = 10, numFmt = fm, halign = "right")
    }
    addStyle(wb, ws, stl, rows = baris + i, cols = j)
  }
  addStyle(wb, ws, st_bawah, rows = baris + nr, cols = 1:nc,
           gridExpand = TRUE, stack = TRUE)
  baris <<- baris + nr + 2
}

# ---- A. Informasi umum ----
judul("ESTIMASI CENSORED LEAST ABSOLUTE DEVIATIONS (CLAD)")
info("Dependent Variable", nama_y)
info("Method",         "Censored LAD - Powell (sensor dua sisi)")
info("Date / Time",    format(Sys.time(), "%m/%d/%Y   %H:%M"))
info("Sample",         paste("1", n))
info("Included obs",   n)
info("Left-censored",  paste0(n_kiri, " obs (Y = 0)"))
info("Right-censored", paste0(n_kanan, " obs (Y >= batas ARA), ",
                              round(100 * n_kanan / n, 2), "%"))
info("Uncensored",     paste0(n_tak, " obs"))
info("Std. errors",    paste0(ket_boot, ", ", R, " replikasi, valid = ", sum(ok)))

# ---- B. Tabel 1: Statistik deskriptif ----
judul("TABEL 1. STATISTIK DESKRIPTIF")
df1 <- data.frame(Label = rownames(tabel_desk), tabel_desk, check.names = FALSE)
names(df1)[1] <- ""
fm1 <- matrix("0.000000", nrow(df1), ncol(df1))
fm1[, 1] <- "TEKS"
fm1[df1[[1]] == "Observations", -1] <- "0"
tulis_tabel(df1, fm1)

# ---- C. Tabel 2: Multikolinearitas ----
judul("TABEL 2. UJI MULTIKOLINEARITAS")
df2 <- data.frame(Variable  = rownames(tabel_vif),
                  VIF       = tabel_vif$VIF,
                  Tolerance = tabel_vif$Tolerance,
                  Keputusan = ifelse(tabel_vif$VIF < 10,
                                     "Bebas multikolinearitas",
                                     "Indikasi multikolinearitas"))
tulis_tabel(df2, c("TEKS", "0.0000", "0.0000", "TEKS"))
catat("Kriteria: VIF < 10 dan Tolerance > 0.10 -> tidak ada multikolinearitas.")
if (ncol(Xp) < 2) catat("Catatan : hanya ada satu variabel independen, sehingga VIF = 1 dan Tolerance = 1.")
baris <- baris + 1

if (ncol(tabel_kor) >= 2) {
  sub_judul("Matriks Korelasi Antar Variabel")
  df_k <- data.frame(Label = rownames(tabel_kor), tabel_kor, check.names = FALSE)
  names(df_k)[1] <- ""
  tulis_tabel(df_k, c("TEKS", rep("0.0000", ncol(tabel_kor))))
}

# ---- D. Tabel 3: Hasil estimasi & uji t ----
judul("TABEL 3. HASIL ESTIMASI CLAD DAN UJI SIGNIFIKANSI PARSIAL (UJI t)")
df3 <- data.frame(Variable      = nm_koef,
                  Coefficient   = as.numeric(beta_hat),
                  `Std. Error`  = as.numeric(se),
                  `t-Statistic` = as.numeric(t_hit),
                  `Prob.`       = as.numeric(p_t),
                  `Sig.`        = bintang(as.numeric(p_t)),
                  check.names = FALSE)
tulis_tabel(df3, c("TEKS", "0.000000", "0.000000", "0.0000", "0.0000", "TEKS_KANAN"))

# Ringkasan model (dua kolom, gaya EViews)
ringkas <- function(l1, v1, l2, v2, fm1 = "0.000000", fm2 = "0.000000") {
  writeData(wb, ws, l1, startCol = 1, startRow = baris)
  writeData(wb, ws, v1, startCol = 2, startRow = baris)
  writeData(wb, ws, l2, startCol = 4, startRow = baris)
  writeData(wb, ws, v2, startCol = 6, startRow = baris)
  addStyle(wb, ws, st_teks, rows = baris, cols = c(1, 4), gridExpand = TRUE)
  addStyle(wb, ws, createStyle(fontName = fn, fontSize = 10, numFmt = fm1,
                               halign = "right"), rows = baris, cols = 2)
  addStyle(wb, ws, createStyle(fontName = fn, fontSize = 10, numFmt = fm2,
                               halign = "right"), rows = baris, cols = 6)
  baris <<- baris + 1
}
ringkas("Pseudo R-squared",          pseudo_R2,   "Mean dependent var",    mean(y))
ringkas("Adjusted pseudo R-squared", adj_R2,      "S.D. dependent var",    sd(y))
ringkas("R-squared (korelasi)",      R2_korelasi, "Objective (V1)",        V1)
ringkas("Observations",              n,           "Restr. objective (V0)", V0, fm1 = "0")
addStyle(wb, ws, st_bawah, rows = baris - 1, cols = 1:6, gridExpand = TRUE, stack = TRUE)
baris <- baris + 1

catat(c("Keterangan: *** p<0.01, ** p<0.05, * p<0.10, ns = tidak signifikan.",
        paste0("Prob. dihitung dari distribusi t dengan df = n - k = ", n - k, "."),
        "Pseudo R-squared = 1 - V1/V0 ; V1 = objektif model, V0 = objektif model konstanta."))
baris <- baris + 1

sub_judul("Interval Kepercayaan 95% (Bootstrap Percentile)")
df_ci <- data.frame(Variable      = nm_koef,
                    Coefficient   = as.numeric(beta_hat),
                    `Lower 2.5%`  = as.numeric(ci[, 1]),
                    `Upper 97.5%` = as.numeric(ci[, 2]),
                    `Mencakup 0?` = ifelse(ci[, 1] <= 0 & ci[, 2] >= 0, "Ya", "Tidak"),
                    check.names = FALSE)
tulis_tabel(df_ci, c("TEKS", "0.000000", "0.000000", "0.000000", "TEKS_KANAN"))

# ---- E. Tabel 4: Uji Wald ----
judul("TABEL 4. UJI SIGNIFIKANSI SIMULTAN (WALD TEST)")
df4 <- data.frame(`Test Statistic` = c("F-statistic", "Chi-square"),
                  Value            = c(wald / df_wald, wald),
                  df               = c(paste0("(", df_wald, ", ", n - k, ")"),
                                       as.character(df_wald)),
                  Probability      = c(p_F, p_wald),
                  check.names = FALSE)
tulis_tabel(df4, c("TEKS", "0.000000", "TEKS_KANAN", "0.0000"))
info("Null Hypothesis", paste0(paste(nm_koef[-1], collapse = " = "), " = 0"))
info("Keputusan (alpha 5%)",
     if (p_wald < 0.05)
       "H0 ditolak -> variabel independen berpengaruh signifikan secara simultan."
     else
       "H0 diterima -> variabel independen tidak berpengaruh signifikan secara simultan.")
baris <- baris + 1
catat(paste0("Sumber: hasil olahan data (CLAD, R). Std. error = ", ket_boot, "."))

# ---- Lebar kolom & simpan ----
setColWidths(wb, ws, cols = 1:kol_n, widths = c(26, rep(16, kol_n - 1)))

nama_xlsx <- "Output_CLAD_Laporan.xlsx"
tryCatch(
  saveWorkbook(wb, nama_xlsx, overwrite = TRUE),
  error = function(e) stop("Gagal menyimpan. Tutup dulu file ", nama_xlsx,
                           " di Excel, lalu jalankan ulang bagian 11.5.")
)
cat("Excel tersimpan di:", file.path(getwd(), nama_xlsx), "\n")

