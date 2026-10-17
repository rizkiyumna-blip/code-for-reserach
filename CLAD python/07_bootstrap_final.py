import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import time
import numpy as np
import pandas as pd
from scipy import stats
from concurrent.futures import ProcessPoolExecutor
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

JUMLAH_REPLIKASI = 1000
SEED = 2024            # HARUS sama dengan SEED di 04_bootstrap.py
N_PEKERJA = 3          # jumlah core; disarankan jumlah core komputermu dikurangi 1
UJI_COBA = False        # True = hanya 8 replikasi untuk tes kecepatan
# =====================================================

_G = {}

def _init(y, X, L, U, b_hat):
    _G.update(y=y, X=X, L=L, U=U, b_hat=b_hat)

def _kerja(arg):
    r, idx, seed = arg
    y, X, L, U, b_hat = _G["y"], _G["X"], _G["L"], _G["U"], _G["b_hat"]
    yb, Xb = y[idx], X[idx]
    Lb = L[idx] if L is not None else None
    Ub = U[idx] if U is not None else None
    h = clad(yb, Xb, Lb, Ub, n_acak=60, n_perturb=20, seed=seed,
             titik_awal_tambahan=b_hat, polish=True)
    return r, h["b"], h["q"]

def main():
    path = "data/" + NAMA_FILE
    if NAMA_FILE.lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(path, sheet_name=NAMA_SHEET)
    else:
        df = pd.read_csv(path, sep=None, engine="python")
    kolom = [Y_COL] + X_COLS + [c for c in [BATAS_KIRI, BATAS_KANAN] if isinstance(c, str)]
    df = df[kolom].dropna().reset_index(drop=True)

    def batas(spec, n):
        if spec is None:
            return None
        if isinstance(spec, str):
            return df[spec].astype(float).to_numpy()
        return np.full(n, float(spec))

    n = len(df)
    y = df[Y_COL].astype(float).to_numpy()
    X = np.column_stack([np.ones(n)] + [df[c].astype(float).to_numpy() for c in X_COLS])
    nama = ["Konstanta"] + X_COLS
    L, U = batas(BATAS_KIRI, n), batas(BATAS_KANAN, n)
    k = X.shape[1]

    b_hat = np.load("output/b_hat.npy")
    B_lama = np.load("output/bootstrap_B.npy")
    R = 8 if UJI_COBA else JUMLAH_REPLIKASI

    # Sampel bootstrap yang sama persis dengan 04_bootstrap.py
    rng = np.random.default_rng(SEED)
    tugas = []
    for r in range(JUMLAH_REPLIKASI):
        idx = rng.integers(0, n, n)
        seed = int(rng.integers(0, 10**9))
        tugas.append((r, idx, seed + 7))
    tugas = tugas[:R]

    file_prog = "output/progres_bootstrap_final.npz"
    B = np.full((R, k), np.nan)
    Q = np.full(R, np.nan)
    if (not UJI_COBA) and os.path.exists(file_prog):
        d = np.load(file_prog)
        if d["B"].shape[0] == R:
            B, Q = d["B"], d["Q"]
            print("Melanjutkan progres:", int(np.sum(~np.isnan(Q))), "replikasi sudah selesai")

    sisa = [t for t in tugas if np.isnan(Q[t[0]])]
    print("=" * 60)
    print(f"BOOTSTRAP TELITI: {len(sisa)} replikasi akan dijalankan, {N_PEKERJA} pekerja")
    t0 = time.time()
    selesai = 0
    ukuran_batch = max(N_PEKERJA * 4, 16)
    with ProcessPoolExecutor(max_workers=N_PEKERJA, initializer=_init,
                             initargs=(y, X, L, U, b_hat)) as ex:
        for s in range(0, len(sisa), ukuran_batch):
            for r, b, q in ex.map(_kerja, sisa[s:s + ukuran_batch]):
                B[r], Q[r] = b, q
            selesai = min(s + ukuran_batch, len(sisa))
            if not UJI_COBA:
                np.savez(file_prog, B=B, Q=Q)
            el = time.time() - t0
            print(f"  {selesai}/{len(sisa)} selesai, sisa sekitar {el/selesai*(len(sisa)-selesai)/60:.1f} menit")

    print("Waktu total (menit):", round((time.time() - t0) / 60, 2))

    # Verifikasi: Q tiap replikasi <= Q(b_hat) pada sampelnya
    salah = 0
    for r, idx, _ in tugas:
        Lb = L[idx] if L is not None else None
        Ub = U[idx] if U is not None else None
        if Q[r] > objective(b_hat, y[idx], X[idx], Lb, Ub) + 1e-7:
            salah += 1
    print("\nVerifikasi: replikasi lebih buruk dari b_hat:", salah, "(harus 0)")

    if UJI_COBA:
        print("Mode uji coba selesai. Ubah UJI_COBA = False untuk versi penuh.")
        return

    np.save("output/bootstrap_B_final.npy", B)

    se = B.std(axis=0, ddof=1)
    z = b_hat / se
    p = 2 * (1 - stats.norm.cdf(np.abs(z)))
    tabel = pd.DataFrame({"Koefisien": b_hat, "SE bootstrap": se, "z": z, "p-value": p,
                          "CI 2.5%": np.percentile(B, 2.5, axis=0),
                          "CI 97.5%": np.percentile(B, 97.5, axis=0)}, index=nama)
    pd.set_option("display.width", 200)
    print("\n" + "=" * 60)
    print("TABEL REGRESI CLAD (BOOTSTRAP FINAL)")
    print(tabel.round(5))

    V = np.cov(B, rowvar=False, ddof=1)
    w = float(b_hat[1:] @ np.linalg.solve(V[1:, 1:], b_hat[1:]))
    print("\nWald simultan:", round(w, 4), "| df = 4 | p =", round(1 - stats.chi2.cdf(w, 4), 6))

    print("\n" + "=" * 60)
    print("DIAGNOSTIK")
    print(pd.DataFrame({"Estimasi": b_hat, "Rata-rata boot": B.mean(axis=0),
                        "Median boot": np.median(B, axis=0)}, index=nama).round(5))
    med = np.median(B, axis=0)
    skala = 1.4826 * np.median(np.abs(B - med), axis=0)
    print("Replikasi ekstrem (>3 skala robust) pada hasil final:",
          int((np.abs(B - med) > 3 * skala).any(axis=1).sum()), "dari", R)
    print("Replikasi yang lebih baik daripada versi ringan awal:",
          int((np.abs(B - B_lama).max(axis=1) > 1e-6).sum()))

main_guard = True

if __name__ == "__main__":
    main()