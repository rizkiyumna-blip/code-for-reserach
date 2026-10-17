import numpy as np
import pandas as pd
from scipy import stats

# =====================================================
# BAGIAN YANG PERLU KAMU UBAH
# =====================================================
NAMA_FILE  = "Book1.xlsx"
NAMA_SHEET = 0
Y_COL      = "Y"
BATAS_KIRI  = 0
BATAS_KANAN = "ARA"
# =====================================================
# Urutan koefisien: Konstanta, X1, X2, Z, X2_Z

path = "data/" + NAMA_FILE
if NAMA_FILE.lower().endswith((".xlsx", ".xls")):
    df = pd.read_excel(path, sheet_name=NAMA_SHEET)
else:
    df = pd.read_csv(path, sep=None, engine="python")

kolom = [Y_COL, "X1", "X2", "Z", "X2_Z"] + [c for c in [BATAS_KIRI, BATAS_KANAN] if isinstance(c, str)]
df = df[kolom].dropna().reset_index(drop=True)
n = len(df)

def batas(spec):
    if isinstance(spec, str):
        return df[spec].astype(float).to_numpy()
    return np.full(n, float(spec))

L, U = batas(BATAS_KIRI), batas(BATAS_KANAN)
x1 = df["X1"].to_numpy(float)
x2 = df["X2"].to_numpy(float)
zz = df["Z"].to_numpy(float)
x2z = df["X2_Z"].to_numpy(float)
Xf = np.column_stack([np.ones(n), x1, x2, zz, x2z])

B = np.load("output/bootstrap_B.npy")
b_hat = np.load("output/b_hat.npy")
R = B.shape[0]

def ringkas(est, boot):
    se = boot.std(ddof=1)
    z = est / se
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return [est, se, z, p, np.percentile(boot, 2.5), np.percentile(boot, 97.5)]

kol = ["Estimasi", "SE bootstrap", "z", "p-value", "CI 2.5%", "CI 97.5%"]
pd.set_option("display.width", 200)

# -----------------------------------------------------
# A. Diagnostik distribusi bootstrap
# -----------------------------------------------------
print("=" * 60)
print("A. DIAGNOSTIK DISTRIBUSI BOOTSTRAP (R =", R, ")")
nama = ["Konstanta", "X1", "X2", "Z", "X2_Z"]
persen = [1, 5, 25, 50, 75, 95, 99]
dist = pd.DataFrame(np.percentile(B, persen, axis=0).T, index=nama,
                    columns=[f"P{p}" for p in persen])
print(dist.round(4))
print("\nSelisih P95-P5 dibanding SE (kalau SE jauh lebih besar, ada pencilan):")
print(pd.DataFrame({"SE": B.std(axis=0, ddof=1),
                    "P95-P5": np.percentile(B, 95, axis=0) - np.percentile(B, 5, axis=0)},
                   index=nama).round(4))
print("\nJumlah replikasi dengan |X1| < 0.001:", int((np.abs(B[:, 1]) < 0.001).sum()), "dari", R)
print("Jumlah replikasi dengan |X1| > 0.3   :", int((np.abs(B[:, 1]) > 0.3).sum()), "dari", R)

# -----------------------------------------------------
# B. Efek marginal pada variabel laten (y*)
# -----------------------------------------------------
print("\n" + "=" * 60)
print("B. EFEK MARGINAL PADA VARIABEL LATEN (y*)")
baris, idx = [], []
for nm, z0 in [("Z = rata-rata - 1 SD", zz.mean() - zz.std()),
               ("Z = rata-rata", zz.mean()),
               ("Z = rata-rata + 1 SD", zz.mean() + zz.std())]:
    est = b_hat[2] + b_hat[4] * z0
    boot = B[:, 2] + B[:, 4] * z0
    baris.append(ringkas(est, boot)); idx.append(f"Efek X2 pada {nm} ({z0:.2f})")
for nm, x0 in [("X2 = rata-rata - 1 SD", x2.mean() - x2.std()),
               ("X2 = rata-rata", x2.mean()),
               ("X2 = rata-rata + 1 SD", x2.mean() + x2.std())]:
    est = b_hat[3] + b_hat[4] * x0
    boot = B[:, 3] + B[:, 4] * x0
    baris.append(ringkas(est, boot)); idx.append(f"Efek Z pada {nm} ({x0:.2f})")
baris.append(ringkas(b_hat[1], B[:, 1])); idx.append("Efek X1 (0 ke 1)")
me_laten = pd.DataFrame(baris, index=idx, columns=kol)
print(me_laten.round(5))

# -----------------------------------------------------
# C. Efek marginal rata-rata (AME) pada median bersyarat y teramati
# -----------------------------------------------------
print("\n" + "=" * 60)
print("C. EFEK MARGINAL RATA-RATA PADA MEDIAN BERSYARAT Y TERAMATI")

def ame(theta):
    xb = Xf @ theta
    ind = ((xb > L) & (xb < U)).astype(float)
    return np.array([ind.mean() * theta[1],
                     np.mean(ind * (theta[2] + theta[4] * zz)),
                     np.mean(ind * (theta[3] + theta[4] * x2)),
                     ind.mean()])

a_hat = ame(b_hat)
a_boot = np.array([ame(b) for b in B])
me_obs = pd.DataFrame([ringkas(a_hat[j], a_boot[:, j]) for j in range(3)],
                      index=["X1", "X2", "Z"], columns=kol)
print(me_obs.round(5))
print("Proporsi observasi dengan prediksi di antara batas kiri dan kanan:", round(a_hat[3], 3))

# -----------------------------------------------------
# D. Uji Wald parsial
# -----------------------------------------------------
print("\n" + "=" * 60)
print("D. UJI WALD PARSIAL (bootstrap)")
V = np.cov(B, rowvar=False, ddof=1)
def wald(ix):
    bs = b_hat[ix]
    w = float(bs @ np.linalg.solve(V[np.ix_(ix, ix)], bs))
    return w, len(ix), 1 - stats.chi2.cdf(w, len(ix))
uji = {"X2 dan X2_Z = 0 (X2 berpengaruh?)": [2, 4],
       "Z dan X2_Z = 0 (Z berpengaruh?)": [3, 4],
       "X2, Z, dan X2_Z = 0": [2, 3, 4],
       "Semua slope = 0": [1, 2, 3, 4]}
tw = pd.DataFrame([wald(v) for v in uji.values()], index=list(uji.keys()),
                  columns=["Wald chi2", "df", "p-value"])
print(tw.round(4))

import os
os.makedirs("output", exist_ok=True)
with pd.ExcelWriter("output/efek_marginal_sementara.xlsx") as w:
    me_laten.round(5).to_excel(w, sheet_name="Efek laten")
    me_obs.round(5).to_excel(w, sheet_name="AME teramati")
    tw.round(4).to_excel(w, sheet_name="Wald parsial")
    dist.round(4).to_excel(w, sheet_name="Distribusi bootstrap")
print("\nHasil disimpan di output/efek_marginal_sementara.xlsx")