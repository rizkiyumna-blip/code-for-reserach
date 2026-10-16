# =====================================================================
# SIMULASI: apakah kode CLAD dua sisi menemukan kembali koefisien sebenarnya?
# Syarat: objek `data_skripsi` ada di Environment (kolom X1, X2, Z, ARA).
# Skrip ini berdiri sendiri; tidak perlu menjalankan clad_lengkap.R dulu.
# =====================================================================

if (!require(quantreg)) install.packages("quantreg")
library(quantreg)

asal <- data_skripsi    # dipakai hanya untuk sebaran X1, X2, Z, dan batas ARA

# ---- 1. PENGATURAN (boleh diubah) -----------------------------------
beta_lereng <- c(X1 = -0.05, X2 = 0.008, Z = 0.003, X2_Z = -0.0002)
sigma       <- 0.15     # simpangan baku galat acak
target_ARA  <- 0.56     # proporsi observasi menyentuh batas ARA yang diinginkan
n_besar     <- 2000     # ukuran sampel besar (membuktikan kode benar)
n_kecil     <- 107      # ukuran sampel seperti datamu
reps_besar  <- 100      # jumlah ulangan untuk sampel besar
reps_kecil  <- 500      # jumlah ulangan untuk sampel kecil
set.seed(2024)

# ---- 2. Fungsi CLAD (sama dengan clad_lengkap.R) --------------------
obj_clad <- function(b, y, X, U) {
  xb   <- as.vector(X %*% b)
  pred <- pmin(U, pmax(0, xb))
  sum(abs(y - pred))
}

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

# Estimasi CLAD dengan standardisasi, hasil dikembalikan ke skala asli
est_clad <- function(y, X, U, n_start = 3) {
  mu  <- colMeans(X[, -1]); sdv <- apply(X[, -1], 2, sd)
  Xs  <- X
  Xs[, -1] <- scale(X[, -1], center = mu, scale = sdv)
  b   <- fit_clad(y, Xs, U, n_start)$par
  out <- b
  out[-1] <- b[-1] / sdv
  out[1]  <- b[1] - sum(b[-1] * mu / sdv)
  out
}

# Pembanding: regresi median biasa (LAD) dan OLS, tanpa koreksi sensor
est_lad <- function(y, X) suppressWarnings(coef(rq.fit.br(X, y, tau = 0.5)))
est_ols <- function(y, X) coef(lm.fit(X, y))

# ---- 3. Pembuat data simulasi ---------------------------------------
# X1, X2, Z, dan batas ARA diambil ulang dari datamu; Y dibuat dari koefisien
ambil_X <- function(n) {
  i  <- sample(nrow(asal), n, replace = TRUE)
  X2 <- asal$X2[i]; Z <- asal$Z[i]
  list(X = cbind(Intersep = 1, X1 = asal$X1[i], X2 = X2, Z = Z, X2_Z = X2 * Z),
       U = asal$ARA[i])
}

buat_data <- function(n, beta) {
  a <- ambil_X(n)
  ystar <- as.vector(a$X %*% beta) + rnorm(n, 0, sigma)   # Y laten (tak teramati)
  y <- pmin(a$U, pmax(0, ystar))                          # Y teramati (tersensor)
  list(y = y, X = a$X, U = a$U, ystar = ystar)
}

# ---- 4. Kalibrasi intersep agar tingkat sensor kanan ~ target_ARA ---
kal <- ambil_X(20000)
m   <- as.vector(kal$X[, -1] %*% beta_lereng)
e   <- rnorm(20000, 0, sigma)
b0  <- as.numeric(quantile(kal$U - m - e, target_ARA))
beta_benar <- c(Intersep = b0, beta_lereng)

cat("Koefisien sebenarnya yang dipakai:\n"); print(round(beta_benar, 5))
cek <- buat_data(20000, beta_benar)
cat("\nKontrol data simulasi (n = 20000):\n",
    " Menyentuh batas ARA :", round(100 * mean(cek$y >= cek$U), 1), "%\n",
    " Menyentuh batas 0   :", round(100 * mean(cek$y <= 0), 1), "%\n",
    " Rentang Y           :", round(range(cek$y), 3), "\n\n")

# ---- 5. Jalankan simulasi -------------------------------------------
jalankan <- function(n, reps) {
  hasil_clad <- hasil_lad <- hasil_ols <-
    matrix(NA, reps, length(beta_benar), dimnames = list(NULL, names(beta_benar)))
  for (r in 1:reps) {
    dat <- buat_data(n, beta_benar)
    hasil_clad[r, ] <- tryCatch(est_clad(dat$y, dat$X, dat$U), error = function(e) NA)
    hasil_lad[r, ]  <- tryCatch(est_lad(dat$y, dat$X),         error = function(e) NA)
    hasil_ols[r, ]  <- tryCatch(est_ols(dat$y, dat$X),         error = function(e) NA)
  }
  list(CLAD = hasil_clad, LAD_naif = hasil_lad, OLS_naif = hasil_ols)
}

ringkas <- function(est, beta) {
  rata <- colMeans(est, na.rm = TRUE)
  data.frame(Benar = beta,
             Rata2 = rata,
             Bias  = rata - beta,
             SD    = apply(est, 2, sd, na.rm = TRUE),
             RMSE  = sqrt(colMeans(sweep(est, 2, beta)^2, na.rm = TRUE)))
}

laporkan <- function(hasil, n, reps) {
  cat("==============================================================\n")
  cat("Sampel n =", n, "| ulangan =", reps, "\n")
  cat("Ulangan CLAD yang gagal:", sum(!complete.cases(hasil$CLAD)), "\n")
  for (nm in names(hasil)) {
    cat("\n--- ", nm, " ---\n", sep = "")
    print(round(ringkas(hasil[[nm]], beta_benar), 5))
  }
  cat("\n")
}

cat("Menjalankan sampel besar...\n")
hasil_besar <- jalankan(n_besar, reps_besar)
laporkan(hasil_besar, n_besar, reps_besar)

cat("Menjalankan sampel kecil...\n")
hasil_kecil <- jalankan(n_kecil, reps_kecil)
laporkan(hasil_kecil, n_kecil, reps_kecil)

# ---- 6. Ringkasan satu tabel (opsional ekspor) ----------------------
# tabel_besar <- cbind(Variabel = names(beta_benar), ringkas(hasil_besar$CLAD, beta_benar))
# tabel_kecil <- cbind(Variabel = names(beta_benar), ringkas(hasil_kecil$CLAD, beta_benar))
# writexl::write_xlsx(list(CLAD_besar = tabel_besar, CLAD_kecil = tabel_kecil), "simulasi_clad.xlsx")