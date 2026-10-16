install.packages("readxl")
library(readxl)
data_skripsi <- read_excel("book1.xlsx")
# =====================================================================
# CLAD dua sisi (kiri = 0, kanan = batas ARA per baris)
# Jalankan dari atas ke bawah. Syarat: objek `data_skripsi` ada di Environment.
# =====================================================================

# ---- 0. Persiapan --------------------------------------------------
if (!require(quantreg)) install.packages("quantreg")
library(quantreg)

# Bersihkan objek lama (data_skripsi dipertahankan)
rm(list = setdiff(ls(), "data_skripsi"))

# ---- 1. Siapkan data -----------------------------------------------
# Kolom ARA berisi batas ARA nominal tiap baris (mis. 0.35 / 0.25)
d <- data_skripsi
toleransi <- 0.02                                  # toleransi pembulatan tick size

d$batas      <- d$ARA
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

# ---- 5. Bootstrap berstrata berdasarkan X1 -------------------------
set.seed(123)
R <- 1000                                # uji coba: 100; hasil akhir: 1000
boot_std <- matrix(NA, R, k, dimnames = list(NULL, colnames(X)))
pesan <- character(0)
g0 <- which(d$X1 == 0); g1 <- which(d$X1 == 1)

for (r in 1:R) {
  idx <- c(sample(g0, length(g0), replace = TRUE),
           sample(g1, length(g1), replace = TRUE))
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

# ---- 11. (Opsional) Ekspor ke CSV ----------------------------------
# write.csv(tabel_koef, "koefisien_clad.csv")
# write.csv(tabel_t,    "uji_t_clad.csv")
# write.csv(tabel_vif,  "vif_clad.csv")