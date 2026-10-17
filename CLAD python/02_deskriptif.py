import numpy as np
import pandas as pd

# =====================================================
# BAGIAN YANG PERLU KAMU UBAH
# =====================================================
NAMA_FILE  = "Book1.xlsx"      # ganti sesuai nama file datamu
NAMA_SHEET = 0

Y_COL      = "Y"                 # variabel dependen
X_COLS     = ["X1", "X2", "Z", "X2_Z"]   # variabel independen

BATAS_KIRI  = 0                  # angka tetap, atau nama kolom (contoh: "LAR"), atau None
BATAS_KANAN = "ARA"              # angka tetap, nama kolom, atau None

TOLERANSI  = 0.0                 # sementara 0; kita tentukan setelah lihat hasil
# =====================================================

path = "data/" + NAMA_FILE
if NAMA_FILE.lower().endswith((".xlsx", ".xls")):
    df = pd.read_excel(path, sheet_name=NAMA_SHEET)
else:
    df = pd.read_csv(path, sep=None, engine="python")

def ambil_batas(spec, n):
    """Ubah pengaturan batas menjadi deret angka per observasi."""
    if spec is None:
        return None
    if isinstance(spec, str):
        return df[spec].astype(float).to_numpy()
    return np.full(n, float(spec))

y = df[Y_COL].astype(float).to_numpy()
L = ambil_batas(BATAS_KIRI, len(df))
U = ambil_batas(BATAS_KANAN, len(df))

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", None)

# -----------------------------------------------------
# 1. CEK KUALITAS DATA
# -----------------------------------------------------
print("=" * 60)
print("1. CEK KUALITAS DATA")
semua = [Y_COL] + X_COLS + [c for c in [BATAS_KIRI, BATAS_KANAN] if isinstance(c, str)]
print("Data kosong per kolom:")
print(df[semua].isnull().sum())

if U is not None:
    print("\nObservasi dengan Y > batas kanan :", int((y > U + 1e-12).sum()))
if L is not None:
    print("Observasi dengan Y < batas kiri  :", int((y < L - 1e-12).sum()))

if "X2" in X_COLS and "Z" in X_COLS and "X2_Z" in X_COLS:
    selisih = np.abs(df["X2"] * df["Z"] - df["X2_Z"]).max()
    print("Selisih terbesar antara X2*Z dan X2_Z :", round(selisih, 6))

# -----------------------------------------------------
# 2. DIAGNOSTIK TOLERANSI SENSOR
# -----------------------------------------------------
print("\n" + "=" * 60)
print("2. JUMLAH OBSERVASI TERSENSOR MENURUT TOLERANSI")
baris = []
for tol in [0, 0.001, 0.002, 0.005, 0.01]:
    kiri  = int((y <= L + tol).sum()) if L is not None else 0
    kanan = int((y >= U - tol).sum()) if U is not None else 0
    baris.append({"Toleransi": tol, "Tersensor kiri": kiri,
                  "Tersensor kanan": kanan,
                  "Tidak tersensor": len(y) - kiri - kanan})
print(pd.DataFrame(baris).to_string(index=False))

# -----------------------------------------------------
# 3. STATISTIK DESKRIPTIF
# -----------------------------------------------------
print("\n" + "=" * 60)
print("3. STATISTIK DESKRIPTIF")
kolom_desk = [Y_COL] + X_COLS + ([BATAS_KANAN] if isinstance(BATAS_KANAN, str) else [])
desk = df[kolom_desk].agg(["count", "mean", "median", "std", "min", "max"]).T
desk.columns = ["N", "Mean", "Median", "Std. Dev.", "Min", "Max"]
print(desk.round(4))

kiri  = int((y <= L + TOLERANSI).sum()) if L is not None else 0
kanan = int((y >= U - TOLERANSI).sum()) if U is not None else 0
print(f"\nDengan TOLERANSI = {TOLERANSI}:")
print("  Tersensor kiri  :", kiri)
print("  Tersensor kanan :", kanan)
print("  Tidak tersensor :", len(y) - kiri - kanan)

import os
os.makedirs("output", exist_ok=True)
desk.round(4).to_excel("output/tabel_deskriptif_sementara.xlsx")
print("\nTabel sementara disimpan di output/tabel_deskriptif_sementara.xlsx")