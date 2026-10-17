import time
import numpy as np
import pandas as pd
from scipy import stats
from clad import clad, objective

# =====================================================
# BAGIAN YANG PERLU KAMU UBAH
# =====================================================
NAMA_FILE  = "Book1.xlsx"
NAMA_SHEET = 0

Y_COL      = "Y"
X_COLS     = ["X1", "X2", "Z", "X2_Z"]

BATAS_KIRI  = 0
BATAS_KANAN = "ARA"

JUMLAH_REPLIKASI = 1000   # setara R = 1000 di R
SEED = 2024               # agar hasil bisa direproduksi
UJI_COBA = False          # True = jalankan hanya 20 replikasi untuk tes cepat
# =====================================================

path = "data/" + NAMA_FILE
if NAMA_FILE.lower().endswith((".xlsx", ".xls")):
    df = pd.read_excel(path, sheet_name=NAMA_SHEET)
else:
    df = pd.read_csv(path, sep=None, engine="python")

kolom_pakai = [Y_COL] + X_COLS + [c for c in [BATAS_KIRI, BATAS_KANAN] if isinstance(c, str)]
df = df[kolom_pakai].dropna().reset_index(drop=True)

def ambil_batas(spec, n):
    if spec is None:
        return None
    if isinstance(spec, str):
        return df[spec].astype(float).to_numpy()
    return np.full(n, float(spec))

n = len(df)
y = df[Y_COL].astype(float).to_numpy()
X = np.column_stack([np.ones(n)] + [df[c].astype(float).to_numpy() for c in X_COLS])
nama = ["Konstanta"] + X_COLS
L = ambil_batas(BATAS_KIRI, n)
U = ambil_batas(BATAS_KANAN, n)
k = X.shape[1]

R = 20 if UJI_COBA else JUMLAH_REPLIKASI

# -----------------------------------------------------
# 1. Estimasi utama (pencarian penuh)
# -----------------------------------------------------
print("=" * 60)
print("1. ESTIMASI UTAMA, N =", n)
hasil = clad(y, X, L, U, n_acak=60)
b_hat = hasil["b"]
print(pd.DataFrame({"Koefisien": b_hat}, index=nama).round(6))

# -----------------------------------------------------
# 2. Bootstrap (pairs bootstrap)
# -----------------------------------------------------
print("\n" + "=" * 60)
print(f"2. BOOTSTRAP, R = {R} replikasi")
rng = np.random.default_rng(SEED)
B = np.zeros((R, k))
t0 = time.time()

for r in range(R):
    idx = rng.integers(0, n, n)
    yb, Xb = y[idx], X[idx]
    Lb = L[idx] if L is not None else None
    Ub = U[idx] if U is not None else None
    # titik awal: hasil utama; ditambah beberapa titik acak
    hb = clad(yb, Xb, Lb, Ub, n_acak=4, n_perturb=2,
              seed=int(rng.integers(0, 10**9)),
              titik_awal_tambahan=b_hat, polish=False)
    B[r] = hb["b"]
    if (r + 1) % 50 == 0 or r == R - 1:
        selesai = time.time() - t0
        sisa = selesai / (r + 1) * (R - r - 1)
        print(f"  {r+1}/{R} selesai, perkiraan sisa {sisa/60:.1f} menit")

print("Waktu total bootstrap (menit):", round((time.time() - t0) / 60, 2))

# -----------------------------------------------------
# 3. Tabel regresi: SE, z, p-value, CI
# -----------------------------------------------------
se = B.std(axis=0, ddof=1)
z = b_hat / se
p = 2 * (1 - stats.norm.cdf(np.abs(z)))
ci_lo = np.percentile(B, 2.5, axis=0)
ci_hi = np.percentile(B, 97.5, axis=0)

def bintang(pv):
    return "***" if pv < 0.01 else "**" if pv < 0.05 else "*" if pv < 0.10 else ""

tabel = pd.DataFrame({
    "Koefisien": b_hat, "SE bootstrap": se, "z": z, "p-value": p,
    "CI 2.5%": ci_lo, "CI 97.5%": ci_hi,
    "Sig.": [bintang(v) for v in p]
}, index=nama)

pd.set_option("display.width", 200)
print("\n" + "=" * 60)
print("3. TABEL REGRESI CLAD")
print(tabel.round(6))

# -----------------------------------------------------
# 4. Uji Wald simultan (semua slope = 0, konstanta tidak diuji)
# -----------------------------------------------------
V = np.cov(B, rowvar=False, ddof=1)
bs = b_hat[1:]
Vs = V[1:, 1:]
wald = float(bs @ np.linalg.solve(Vs, bs))
df_w = k - 1
p_wald = 1 - stats.chi2.cdf(wald, df_w)

print("\n" + "=" * 60)
print("4. UJI WALD SIMULTAN (H0: semua koefisien slope = 0)")
print("Statistik Wald (chi2):", round(wald, 4))
print("Derajat bebas        :", df_w)
print("p-value              :", round(p_wald, 6))

# -----------------------------------------------------
# 5. Diagnostik bootstrap
# -----------------------------------------------------
print("\n" + "=" * 60)
print("5. DIAGNOSTIK BOOTSTRAP")
print("Rata-rata bootstrap vs estimasi utama:")
print(pd.DataFrame({"Estimasi": b_hat, "Rata-rata bootstrap": B.mean(axis=0),
                    "Median bootstrap": np.median(B, axis=0)}, index=nama).round(6))

# Simpan hasil mentah agar tidak perlu mengulang bootstrap
import os
os.makedirs("output", exist_ok=True)
np.save("output/bootstrap_B.npy", B)
np.save("output/b_hat.npy", b_hat)
tabel.round(6).to_excel("output/tabel_regresi_sementara.xlsx")
print("\nHasil disimpan di folder output.")