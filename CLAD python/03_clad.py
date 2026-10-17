import time
import numpy as np
import pandas as pd
from clad import clad, objective

# =====================================================
# BAGIAN YANG PERLU KAMU UBAH
# =====================================================
NAMA_FILE  = "Book1.xlsx"
NAMA_SHEET = 0

Y_COL      = "Y"
X_COLS     = ["X1", "X2", "Z", "X2_Z"]

BATAS_KIRI  = 0          # angka, nama kolom, atau None
BATAS_KANAN = "ARA"      # angka, nama kolom, atau None

TOLERANSI  = 0.0         # isi sesuai pilihanmu di langkah 5 (hanya untuk pelaporan jumlah sensor)
JUMLAH_TITIK_AWAL_ACAK = 60
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

print("=" * 60)
print("ESTIMASI CLAD, N =", n)
t0 = time.time()
hasil = clad(y, X, L, U, n_acak=JUMLAH_TITIK_AWAL_ACAK)
waktu = time.time() - t0
b, q = hasil["b"], hasil["q"]

print("\nKOEFISIEN CLAD")
print(pd.DataFrame({"Koefisien": b}, index=nama).round(6))

# Diagnostik titik awal
qs = hasil["semua_q"]
sama = int(np.sum(np.abs(qs - q) <= 1e-6 * max(abs(q), 1e-12)))
print("\n" + "=" * 60)
print("DIAGNOSTIK KONVERGENSI")
print("Jumlah titik awal         :", len(qs))
print("Mencapai nilai terbaik    :", sama)
print("Nilai fungsi tujuan terbaik:", round(q, 8))
print("Nilai terburuk dari semua titik awal:", round(qs.max(), 8))
print("Waktu satu estimasi (detik):", round(waktu, 2))
print("Nilai terbaik tahap 1 (titik awal)  :", round(hasil["q_tahap1"], 8))
print("Nilai terbaik setelah tahap 2 dan 3 :", round(hasil["q_tahap2"], 8))

# Pseudo-R2 (model hanya konstanta sebagai pembanding)
X0 = np.ones((n, 1))
q0 = clad(y, X0, L, U, n_acak=5)["q"]
print("\nPseudo-R2 (1 - Q_model / Q_konstanta):", round(1 - q / q0, 4))

# Jumlah observasi tersensor
kiri  = int((y <= L + TOLERANSI).sum()) if L is not None else 0
kanan = int((y >= U - TOLERANSI).sum()) if U is not None else 0
print("\nTersensor kiri / kanan / tidak:", kiri, "/", kanan, "/", n - kiri - kanan)

# Cek logika: apakah beta CLAD lebih baik dari OLS menurut fungsi tujuan CLAD?
ols = np.linalg.lstsq(X, y, rcond=None)[0]
print("\nSanity check, Q pada OLS :", round(objective(ols, y, X, L, U), 8))
print("Sanity check, Q pada CLAD:", round(q, 8), "(harus lebih kecil atau sama)")