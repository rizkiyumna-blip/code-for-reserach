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

SEED = 2024          # HARUS sama dengan SEED di 04_bootstrap.py
AMBANG = 3.0         # replikasi ditandai jika menyimpang > 3 skala robust dari median
UJI_COBA = False      # True = hanya 10 replikasi ditandai pertama (tes cepat)
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

B = np.load("output/bootstrap_B.npy")
b_hat = np.load("output/b_hat.npy")
R = B.shape[0]

# -----------------------------------------------------
# 1. Buat ulang sampel bootstrap yang sama
# -----------------------------------------------------
rng = np.random.default_rng(SEED)
IDX, SEEDS = [], []
for r in range(R):
    IDX.append(rng.integers(0, n, n))
    SEEDS.append(int(rng.integers(0, 10**9)))

# Verifikasi: tiap replikasi harus punya Q <= Q(b_hat) pada sampelnya sendiri
salah = 0
for r in range(R):
    i = IDX[r]
    Lb = L[i] if L is not None else None
    Ub = U[i] if U is not None else None
    if objective(B[r], y[i], X[i], Lb, Ub) > objective(b_hat, y[i], X[i], Lb, Ub) + 1e-7:
        salah += 1
print("=" * 60)
print("1. VERIFIKASI SAMPEL BOOTSTRAP")
print("Replikasi yang tidak cocok:", salah, "dari", R, "(harus 0)")
if salah > 0:
    raise SystemExit("Sampel tidak cocok. Periksa SEED dan jangan lanjut.")

# -----------------------------------------------------
# 2. Tandai replikasi ekstrem
# -----------------------------------------------------
med = np.median(B, axis=0)
skala = 1.4826 * np.median(np.abs(B - med), axis=0)
ditandai = np.where((np.abs(B - med) > AMBANG * skala).any(axis=1))[0]
print("\n" + "=" * 60)
print("2. REPLIKASI EKSTREM")
print(f"Ditandai (> {AMBANG} skala robust): {len(ditandai)} dari {R}")
if UJI_COBA:
    ditandai = ditandai[:10]
    print("Mode uji coba: hanya memeriksa", len(ditandai), "replikasi")

# -----------------------------------------------------
# 3. Pencarian teliti pada replikasi yang ditandai
# -----------------------------------------------------
print("\n" + "=" * 60)
print("3. PENCARIAN ULANG (teliti)")
B2 = B.copy()
membaik = 0
turun = []
t0 = time.time()
for m, r in enumerate(ditandai):
    i = IDX[r]
    yb, Xb = y[i], X[i]
    Lb = L[i] if L is not None else None
    Ub = U[i] if U is not None else None
    q_lama = objective(B[r], yb, Xb, Lb, Ub)
    h = clad(yb, Xb, Lb, Ub, n_acak=60, n_perturb=20,
             seed=SEEDS[r] + 1, titik_awal_tambahan=b_hat, polish=True)
    if h["q"] < q_lama - 1e-6:
        B2[r] = h["b"]
        membaik += 1
        turun.append(q_lama - h["q"])
    if (m + 1) % 10 == 0 or m == len(ditandai) - 1:
        sisa = (time.time() - t0) / (m + 1) * (len(ditandai) - m - 1)
        print(f"  {m+1}/{len(ditandai)} diperiksa, sisa sekitar {sisa/60:.1f} menit")

print("\nReplikasi yang ternyata nyangkut di minimum lokal:", membaik, "dari", len(ditandai))
if turun:
    print("Rata-rata penurunan Q:", round(float(np.mean(turun)), 6),
          "| terbesar:", round(float(np.max(turun)), 6))

# -----------------------------------------------------
# 4. Perbandingan sebelum dan sesudah
# -----------------------------------------------------
def tabel(Bx):
    se = Bx.std(axis=0, ddof=1)
    z = b_hat / se
    p = 2 * (1 - stats.norm.cdf(np.abs(z)))
    return se, p

se1, p1 = tabel(B)
se2, p2 = tabel(B2)
pd.set_option("display.width", 200)
print("\n" + "=" * 60)
print("4. SEBELUM VS SESUDAH PERBAIKAN")
print(pd.DataFrame({"Koefisien": b_hat, "SE lama": se1, "SE baru": se2,
                    "p lama": p1, "p baru": p2,
                    "CI 2.5% baru": np.percentile(B2, 2.5, axis=0),
                    "CI 97.5% baru": np.percentile(B2, 97.5, axis=0)},
                   index=nama).round(5))

V = np.cov(B2, rowvar=False, ddof=1)
w = float(b_hat[1:] @ np.linalg.solve(V[1:, 1:], b_hat[1:]))
print("\nWald simultan baru:", round(w, 4), "| df = 4 | p =", round(1 - stats.chi2.cdf(w, 4), 6))

if not UJI_COBA:
    np.save("output/bootstrap_B_diperbaiki.npy", B2)
    print("\nDisimpan: output/bootstrap_B_diperbaiki.npy")