# ---- PENGATURAN ----------
toleransi <- 0      # toleransi pembulatan tick size untuk deteksi ARA
n_start   <- 50     # titik awal dekat OLS (estimasi utama)
n_random  <- 100    # titik awal acak yang disempurnakan (estimasi utama)
R         <- 1000   # jumlah ulangan bootstrap (uji coba cepat: 20 atau 100)

# ---- 0. Persiapan --------------------------------------------------
if (!require(readxl))   install.packages("readxl")
if (!require(quantreg)) install.packages("quantreg")
library(readxl)
library(quantreg)

data_skripsi <- read_excel("book1.xlsx")

# Bersihkan objek lama (data_skripsi dan pengaturan dipertahankan)
rm(list = setdiff(ls(), c("data_skripsi", "toleransi", "n_start",
                          "n_random", "R")))

# ---- 1. Siapkan data -----------------------------------------------
# Kolom ARA berisi batas ARA nominal tiap baris (mis. 0.35 / 0.25)
d <- as.data.frame(data_skripsi)
kolom_wajib <- c("Y", "ARA", "X1", "X2", "Z", "X2_Z")
if (!all(kolom_wajib %in% names(d)))
  stop("Kolom tidak ditemukan: ",
       paste(setdiff(kolom_wajib, names(d)), collapse = ", "))

# Buang baris yang punya data kosong (NA)
n_awal <- nrow(d)
d <- d[complete.cases(d[, kolom_wajib]), ]
cat("Baris dibuang karena NA:", n_awal - nrow(d), "dari", n_awal, "\n")

d$batas      <- d$ARA
d$sentuh_ARA <- d$Y >= d$batas - toleransi         # TRUE = kena ARA
d$U          <- ifelse(d$sentuh_ARA, d$Y, d$batas) # batas kanan per observasi

# Ringkasan ARA
tabel_ara <- data.frame(
  Status = c("Tidak kena ARA", "Kena ARA"),
  Jumlah = as.vector(table(factor(d$sentuh_ARA, levels = c(FALSE, TRUE)))))
tabel_ara$Persen <- round(100 * tabel_ara$Jumlah / nrow(d), 2)
cat("Toleransi ARA yang dipakai:", toleransi, "\n")
print(tabel_ara)
cat("Y melebihi batas nominal (harus 0):", sum(d$Y > d$batas), "\n")
cat("Y tepat bernilai 0 (sensor kiri)  :", sum(d$Y == 0), "\n")
cat("Y terkecil (harus >= 0)           :", min(d$Y), "\n\n")

# ---- 2. Matriks desain ---------------------------------------------
y <- d$Y
U <- d$U
X <- cbind(Intersep = 1, X1 = d$X1, X2 = d$X2, Z = d$Z, X2_Z = d$X2_Z)
n <- nrow(X); k <- ncol(X)

# Standardisasi prediktor (hanya untuk stabilitas optimasi)
mu  <- colMeans(X[, -1])
sdv <- apply(X[, -1], 2, sd)
Xs  <- X
Xs[, -1] <- scale(X[, -1], center = mu, scale = sdv)

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

# Beberapa titik awal dekat OLS, ambil yang objektifnya paling kecil
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

# Penyempurnaan (polish): Nelder-Mead dan iterasi Powell bergantian
polish <- function(b, y, X, U, rounds = 10) {
  f <- function(z) obj_clad(z, y, X, U)
  best_b <- b
  best_v <- f(b)
  for (r in 1:rounds) {
    v_awal <- best_v
    o <- optim(best_b, f, method = "Nelder-Mead",
               control = list(maxit = 5000, reltol = 1e-14))
    kandidat <- list(o$par, clad_iter(y, X, U, o$par))
    for (cb in kandidat) {
      v <- f(cb)
      if (v < best_v) { best_v <- v; best_b <- cb }
    }
    if (best_v > v_awal - 1e-10) break      # berhenti jika tidak membaik
  }
  list(par = best_b, value = best_v)
}

# Estimator utama: banyak titik awal (dekat dan jauh), semuanya di-polish
fit_clad2 <- function(y, X, U, n_start = 20, n_random = 20,
                      extra_starts = NULL) {
  p0     <- fit_clad(y, X, U, n_start = n_start)
  b_ols  <- coef(lm.fit(X, y))
  kk     <- ncol(X)
  sdg    <- c(0.2, 0.5, 1)                  # jarak titik awal acak
  starts <- c(list(p0$par, b_ols), extra_starts,
              lapply(seq_len(n_random),
                     function(i) p0$par + rnorm(kk, 0, sdg[i %% 3 + 1])))
  hasil <- lapply(starts, function(s) polish(s, y, X, U))
  hasil[[which.min(sapply(hasil, `[[`, "value"))]]
}

# ---- 4. Estimasi utama ---------------------------------------------
set.seed(123)
fit      <- fit_clad2(y, Xs, U, n_start = n_start, n_random = n_random)
beta_std <- fit$par                      # skala standar
beta_hat <- ke_asli(fit$par)             # skala asli

cat("\n=== ESTIMASI UTAMA ===\n")
cat("Nilai objektif minimum:", round(fit$value, 5), "\n")
cat("Koefisien (skala asli):\n"); print(round(beta_hat, 5))

# Proporsi observasi yang berkontribusi pada koefisien
xb_hat <- as.vector(Xs %*% beta_std)
cat("\nProporsi observasi aktif (prediksi di dalam rentang):",
    round(mean(xb_hat > 0 & xb_hat < U), 4), "\n")
cat("Proporsi prediksi di bawah 0     :", round(mean(xb_hat <= 0), 4), "\n")
cat("Proporsi prediksi di atas batas  :", round(mean(xb_hat >= U), 4), "\n")

# ---- 5. Bootstrap berstrata berdasarkan X1 -------------------------
set.seed(123)
boot_std <- matrix(NA, R, k, dimnames = list(NULL, colnames(X)))
pesan <- character(0)
g0 <- which(d$X1 == 0); g1 <- which(d$X1 == 1)
waktu_mulai <- Sys.time()

for (r in 1:R) {
  idx <- c(sample(g0, length(g0), replace = TRUE),
           sample(g1, length(g1), replace = TRUE))
  hasil_r <- tryCatch(
    fit_clad2(y[idx], Xs[idx, ], U[idx], n_start = 2, n_random = 3,
              extra_starts = list(beta_std))$par,
    error = function(e) { pesan <<- c(pesan, conditionMessage(e)); NULL })
  if (!is.null(hasil_r)) boot_std[r, ] <- hasil_r
  if (r %% 10 == 0) {
    lama <- as.numeric(difftime(Sys.time(), waktu_mulai, units = "mins"))
    cat("Bootstrap", r, "dari", R, "| berjalan", round(lama, 1), "menit",
        "| perkiraan total", round(lama * R / r, 1), "menit\n")
  }
}

ok <- complete.cases(boot_std)
cat("\nSampel bootstrap gagal:", sum(!ok), "dari", R,
    "(", round(100 * sum(!ok) / R, 2), "% )\n")
if (length(pesan) > 0) { cat("Pesan error:\n"); print(table(pesan)) }
if (sum(ok) < 10)
  stop("Terlalu banyak bootstrap gagal. Lihat tabel pesan error di atas.")

# Kembalikan hasil bootstrap ke skala asli
boot_asli <- t(apply(boot_std, 1, ke_asli))
colnames(boot_asli) <- colnames(X)
se <- apply(boot_asli[ok, ], 2, sd)
ci <- t(apply(boot_asli[ok, ], 2, quantile, probs = c(0.025, 0.975)))

# ---- 6. Uji multikolinearitas (VIF & Tolerance) --------------------
Xp  <- X[, -1]
vif <- sapply(seq_len(ncol(Xp)), function(j) {
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
V_boot  <- cov(boot_asli[ok, ])[-1, -1]
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

# ---- 11. OUTPUT ------------------------------
# Keluaran: tampil di console + file Output_CLAD.txt + Output_CLAD.xlsx
# (file tersimpan di folder kerja; lihat lokasinya dengan getwd()).

nama_file <- "Output_CLAD"      # nama dasar file keluaran

# ---- 11a. Fungsi bantu (tampilan) -----------------------------------
garis  <- function(ch = "-") cat(strrep(ch, 80), "\n", sep = "")
judul  <- function(teks) { cat("\n"); garis("-"); cat(teks, "\n", sep = ""); garis("-") }
bintang <- function(p) as.character(cut(p, c(-Inf, .01, .05, .10, Inf),
                                        labels = c("***", "**", "*", "ns")))
fmt_baris <- function(x, dg) mapply(function(v, dd)
  formatC(v, format = "f", digits = dd), x, dg)
cetak_tab <- function(df, dg = 4) {
  dg <- rep_len(dg, ncol(df))
  for (j in seq_along(df))
    if (is.numeric(df[[j]]))
      df[[j]] <- formatC(df[[j]], format = "f", digits = dg[j])
  print(df, right = TRUE)
}

# ---- 11b. Hitung semua tabel ---------------------------------------
bh <- setNames(as.numeric(beta_hat), colnames(X))   # koefisien bernama

# Tabel 1: Statistik deskriptif (definisi sama dengan EViews)
vars <- data.frame(Y = d$Y, X1 = d$X1, X2 = d$X2, Z = d$Z, X2_Z = d$X2_Z)
desk_fun <- function(v) {
  N <- length(v); dv <- v - mean(v)
  m2 <- mean(dv^2); m3 <- mean(dv^3); m4 <- mean(dv^4)
  S <- m3 / m2^1.5; K <- m4 / m2^2
  JB <- N / 6 * (S^2 + (K - 3)^2 / 4)
  c("Mean" = mean(v), "Median" = median(v), "Maximum" = max(v),
    "Minimum" = min(v), "Std. Dev." = sd(v), "Skewness" = S,
    "Kurtosis" = K, "Jarque-Bera" = JB,
    "Probability" = pchisq(JB, 2, lower.tail = FALSE),
    "Observations" = N)
}
tabel_desk <- sapply(vars, desk_fun)         # baris = statistik, kolom = variabel
dg_desk    <- c(4, 4, 4, 4, 4, 4, 4, 4, 4, 0)

# Tabel 2: Distribusi sensor
tabel_sensor <- data.frame(
  Jumlah = c(n, sum(y == 0), sum(d$sentuh_ARA), sum(y > 0 & !d$sentuh_ARA)),
  row.names = c("Total observasi", "Tersensor kiri (Y = 0)",
                "Tersensor kanan (kena ARA)", "Tidak tersensor"))
tabel_sensor$Persen <- 100 * tabel_sensor$Jumlah / n

# Tabel 3: Hasil estimasi (koefisien, SE bootstrap, uji t)
tabel_reg <- data.frame(
  unname(bh), as.numeric(se), as.numeric(t_hit), as.numeric(p_t),
  as.numeric(ci[, 1]), as.numeric(ci[, 2]), bintang(as.numeric(p_t)),
  row.names = colnames(X))
colnames(tabel_reg) <- c("Coefficient", "Std. Error", "t-Statistic", "Prob.",
                         "CI 2.5%", "CI 97.5%", "Sig.")
dg_reg <- c(6, 6, 4, 4, 6, 6, 0)

# Statistik kecocokan model
xb_akhir  <- as.vector(Xs %*% beta_std)
aktif_idx <- xb_akhir > 0 & xb_akhir < U
nilai_fit <- c(pseudo_R2, R2_korelasi, V1, V0, mean(y), sd(y),
               mean(aktif_idx), n, sum(ok))
nama_fit  <- c("Pseudo R-squared (1 - V1/V0)",
               "R-squared korelasi (aktual vs prediksi)",
               "Objective value model penuh (V1)",
               "Objective value model nol (V0)",
               "Mean dependent var", "S.D. dependent var",
               "Proporsi observasi aktif", "Observations",
               "Replikasi bootstrap berhasil")
dg_fit    <- c(4, 4, 4, 4, 4, 4, 4, 0, 0)
tabel_fit <- data.frame(Nilai = nilai_fit, row.names = nama_fit)

# Tabel 4: Uji Wald simultan
H0_teks <- paste(paste0(colnames(X)[-1], " = 0"), collapse = ", ")
F_wald  <- wald / df_wald
tabel_wald <- data.frame(
  "Value"       = c(F_wald, wald),
  "df"          = c(paste0("(", df_wald, ", ", n - k, ")"), as.character(df_wald)),
  "Probability" = c(pf(F_wald, df_wald, n - k, lower.tail = FALSE), p_wald),
  row.names = c("F-statistic", "Chi-square"), check.names = FALSE)

# Tabel 5: Efek kondisional (untuk interpretasi interaksi X2_Z)
efek_kond <- function(utama, moderator) {
  m   <- d[[moderator]]
  lev <- mean(m) + c(-1, 0, 1) * sd(m)
  bb  <- boot_asli[ok, , drop = FALSE]
  est <- as.numeric(bh[utama] + bh["X2_Z"] * lev)
  sim <- sapply(lev, function(v) bb[, utama] + bb[, "X2_Z"] * v)
  se_ <- apply(sim, 2, sd)
  tt  <- est / se_
  pp  <- 2 * pt(-abs(tt), df = n - k)
  out <- data.frame(
    lev, est, se_, tt, pp,
    apply(sim, 2, quantile, probs = 0.025),
    apply(sim, 2, quantile, probs = 0.975), bintang(pp),
    row.names = paste0(moderator, c(" = Mean - 1 SD", " = Mean", " = Mean + 1 SD")))
  colnames(out) <- c(paste("Nilai", moderator), paste("Efek", utama), "Std. Error",
                     "t-Statistic", "Prob.", "CI 2.5%", "CI 97.5%", "Sig.")
  out
}
tabel_efek_X2 <- efek_kond("X2", "Z")
tabel_efek_Z  <- efek_kond("Z", "X2")
dg_efek <- c(4, 6, 6, 4, 4, 6, 6, 0)

# Tabel 6: Multikolinearitas
tabel_vif2 <- data.frame(
  VIF = tabel_vif$VIF, Tolerance = tabel_vif$Tolerance,
  Keterangan = ifelse(tabel_vif$VIF < 5, "Rendah",
                      ifelse(tabel_vif$VIF < 10, "Sedang", "Tinggi")),
  row.names = rownames(tabel_vif))

# Tabel 7: Informasi bootstrap
tabel_boot <- data.frame(
  Nilai = c(R, sum(ok), sum(!ok), length(g0), length(g1)),
  row.names = c("Replikasi bootstrap", "Replikasi berhasil", "Replikasi gagal",
                "Observasi grup X1 = 0", "Observasi grup X1 = 1"))

# ---- 11c. Cetak semua ke console -----------------------------------
cetak_semua <- function() {
  garis("=")
  cat("HASIL ESTIMASI CLAD DUA SISI (CENSORED LEAST ABSOLUTE DEVIATIONS)\n")
  garis("=")
  cat("Dependent Variable   : Y\n")
  cat("Method               : Censored LAD dua sisi (Powell)\n")
  cat("                       sensor kiri = 0, sensor kanan = batas ARA\n")
  cat("Included observations: ", n, "\n", sep = "")
  cat("Toleransi ARA        : ", toleransi, "\n", sep = "")
  cat("Standard errors      : Bootstrap berstrata menurut X1 (",
      sum(ok), " dari ", R, " replikasi berhasil)\n", sep = "")
  
  judul("TABEL 1. STATISTIK DESKRIPTIF")
  fd <- apply(tabel_desk, 2, fmt_baris, dg = dg_desk)
  rownames(fd) <- rownames(tabel_desk)
  print(fd, quote = FALSE, right = TRUE)
  cat("Catatan: X1 bervariabel dummy (0/1) dan X2_Z adalah hasil kali X2 dan Z;\n",
      "         skewness, kurtosis, dan Jarque-Bera untuk X1 tidak bermakna.\n", sep = "")
  
  judul("TABEL 2. DISTRIBUSI SENSOR")
  cetak_tab(tabel_sensor, c(0, 2))
  
  judul("TABEL 3. HASIL ESTIMASI CLAD (UJI t PARSIAL)")
  cetak_tab(tabel_reg, dg_reg)
  cat("Keterangan: *** p<0.01, ** p<0.05, * p<0.10, ns = tidak signifikan\n")
  cat("Std. Error dan CI dari bootstrap; p-value memakai distribusi t (df = ",
      n - k, ").\n", sep = "")
  garis("-")
  print(data.frame(Nilai = fmt_baris(nilai_fit, dg_fit),
                   row.names = nama_fit), right = TRUE)
  
  judul("TABEL 4. UJI WALD (UJI SIMULTAN)")
  cat("Null Hypothesis: ", H0_teks, "\n", sep = "")
  cetak_tab(tabel_wald, 4)
  cat("Keputusan: ", ifelse(p_wald < 0.05, "H0 ditolak pada alpha 5%",
                            "H0 tidak ditolak pada alpha 5%"), "\n", sep = "")
  
  judul("TABEL 5. EFEK KONDISIONAL (INTERPRETASI INTERAKSI X2 x Z)")
  cat("(a) Efek X2 pada median variabel laten, pada tiga tingkat Z\n")
  cetak_tab(tabel_efek_X2, dg_efek)
  cat("\n(b) Efek Z pada median variabel laten, pada tiga tingkat X2\n")
  cetak_tab(tabel_efek_Z, dg_efek)
  
  judul("TABEL 6. UJI MULTIKOLINEARITAS (VIF)")
  cetak_tab(tabel_vif2, 4)
  cat("Kriteria: VIF < 5 rendah; 5-10 sedang; > 10 tinggi (Tolerance < 0.10).\n")
  cat("VIF tinggi pada variabel interaksi (X2_Z) lazim dan tidak dianggap masalah.\n")
  
  judul("TABEL 7. INFORMASI BOOTSTRAP")
  cetak_tab(tabel_boot, 0)
  
  cat("\n")
  garis("=")
  cat("Catatan interpretasi: koefisien CLAD adalah efek marginal terhadap median\n")
  cat("variabel laten (sebelum disensor), bukan terhadap Y yang teramati.\n")
  garis("=")
}

cetak_semua()

# ---- 11d. Simpan ke file teks ---------------------------------------
writeLines(capture.output(cetak_semua()), paste0(nama_file, ".txt"))

# ---- 11e. Simpan ke Excel (satu sheet per tabel) --------------------
ada_xlsx <- suppressWarnings(require(openxlsx, quietly = TRUE))
if (!ada_xlsx) {
  try(install.packages("openxlsx"), silent = TRUE)
  ada_xlsx <- suppressWarnings(require(openxlsx, quietly = TRUE))
}

daftar_tabel <- list(
  "Deskriptif"       = list("Tabel 1. Statistik Deskriptif", as.data.frame(tabel_desk)),
  "Sensor"           = list("Tabel 2. Distribusi Sensor", tabel_sensor),
  "Regresi_CLAD"     = list("Tabel 3. Hasil Estimasi CLAD (Uji t Parsial)", tabel_reg),
  "Kecocokan_Model"  = list("Statistik Kecocokan Model", tabel_fit),
  "Uji_Wald"         = list(paste0("Tabel 4. Uji Wald; H0: ", H0_teks), tabel_wald),
  "Efek_Kondisional_X2" = list("Tabel 5a. Efek X2 pada tingkat Z", tabel_efek_X2),
  "Efek_Kondisional_Z"  = list("Tabel 5b. Efek Z pada tingkat X2", tabel_efek_Z),
  "VIF"              = list("Tabel 6. Uji Multikolinearitas", tabel_vif2),
  "Bootstrap"        = list("Tabel 7. Informasi Bootstrap", tabel_boot))

if (ada_xlsx) {
  wb <- createWorkbook()
  gaya_judul <- createStyle(textDecoration = "bold")
  gaya_kepala <- createStyle(textDecoration = "bold", border = "bottom")
  for (nm in names(daftar_tabel)) {
    judul_s <- daftar_tabel[[nm]][[1]]
    df_s    <- daftar_tabel[[nm]][[2]]
    df_s[] <- lapply(df_s, function(x) if (is.numeric(x)) round(x, 6) else x)
    addWorksheet(wb, nm)
    writeData(wb, nm, judul_s, startRow = 1, startCol = 1)
    writeData(wb, nm, df_s, startRow = 3, startCol = 1, rowNames = TRUE)
    addStyle(wb, nm, gaya_judul, rows = 1, cols = 1)
    addStyle(wb, nm, gaya_kepala, rows = 3, cols = 1:(ncol(df_s) + 1),
             gridExpand = TRUE)
    setColWidths(wb, nm, cols = 1:(ncol(df_s) + 1), widths = "auto")
  }
  saveWorkbook(wb, paste0(nama_file, ".xlsx"), overwrite = TRUE)
  cat("\nFile tersimpan di:", getwd(), "\n",
      " -", paste0(nama_file, ".txt"), "\n",
      " -", paste0(nama_file, ".xlsx"), "\n")
} else {
  # Cadangan: simpan tiap tabel sebagai CSV jika openxlsx tidak tersedia
  for (nm in names(daftar_tabel))
    write.csv(daftar_tabel[[nm]][[2]], paste0(nama_file, "_", nm, ".csv"))
  cat("\nPaket openxlsx tidak tersedia; tabel disimpan sebagai CSV di:",
      getwd(), "\n")
}

