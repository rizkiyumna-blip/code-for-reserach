import warnings
import numpy as np
import statsmodels.api as sm
from scipy.optimize import minimize


def sensor_clip(xb, L, U):
    """Terapkan sensor: prediksi dipotong pada batas kiri dan kanan."""
    z = xb.copy()
    if L is not None:
        z = np.maximum(z, L)
    if U is not None:
        z = np.minimum(z, U)
    return z


def objective(b, y, X, L, U):
    """Fungsi tujuan CLAD: jumlah selisih absolut."""
    return np.abs(y - sensor_clip(X @ b, L, U)).sum()


def _satu_titik_awal(y, X, L, U, b0, max_iter=300, tol=1e-8, sabar=25):
    k = X.shape[1]
    b = b0.copy()
    best_b, best_q = b.copy(), objective(b, y, X, L, U)
    tanpa_perbaikan = 0

    for _ in range(max_iter):
        xb = X @ b
        interior = np.ones(len(y), dtype=bool)
        if L is not None:
            interior &= xb > L
        if U is not None:
            interior &= xb < U
        if interior.sum() <= k:
            break
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = sm.QuantReg(y[interior], X[interior]).fit(q=0.5, max_iter=2000)
        except Exception:
            break

        b_baru = np.asarray(res.params, dtype=float)
        q_baru = objective(b_baru, y, X, L, U)

        if q_baru < best_q - 1e-12:
            best_q, best_b = q_baru, b_baru.copy()
            tanpa_perbaikan = 0
        else:
            tanpa_perbaikan += 1

        if np.max(np.abs(b_baru - b)) < tol or tanpa_perbaikan >= sabar:
            break
        b = b_baru

    return best_b, best_q


def clad(y, X, L=None, U=None, n_acak=60, n_perturb=20, seed=123,
         titik_awal_tambahan=None, polish=True):
    """
    Estimasi CLAD dengan sensor kiri dan/atau kanan.
    Kolom pertama X harus konstanta (kolom berisi angka 1).
    Tahap 1: banyak titik awal. Tahap 2: perturbasi dari solusi terbaik.
    Tahap 3: pemolesan Nelder-Mead.
    """
    n, k = X.shape

    # --- Standardisasi internal (hasil dikembalikan ke skala asli) ---
    T = np.eye(k)
    for j in range(1, k):
        mu = X[:, j].mean()
        sd = X[:, j].std()
        if sd == 0:
            sd = 1.0
        T[j, j] = 1.0 / sd
        T[0, j] = -mu / sd
    Xs = X @ T

    # --- Titik awal tahap 1 ---
    ols = np.linalg.lstsq(Xs, y, rcond=None)[0]
    starts = [ols]

    mask = np.ones(n, dtype=bool)
    if L is not None:
        mask &= y > L
    if U is not None:
        mask &= y < U
    if mask.sum() > k:
        starts.append(np.linalg.lstsq(Xs[mask], y[mask], rcond=None)[0])

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            starts.append(np.asarray(sm.QuantReg(y, Xs).fit(q=0.5, max_iter=2000).params, dtype=float))
    except Exception:
        pass

    if titik_awal_tambahan is not None:
        # titik awal tambahan diberikan pada skala asli, ubah ke skala standar
        starts.append(np.linalg.solve(T, np.asarray(titik_awal_tambahan, dtype=float)))

    rng = np.random.default_rng(seed)
    resid = y - Xs @ ols
    sigma2 = resid @ resid / max(n - k, 1)
    se = np.sqrt(np.diag(sigma2 * np.linalg.pinv(Xs.T @ Xs)))
    for _ in range(n_acak):
        starts.append(ols + rng.normal(size=k) * se * 3)

    hasil1 = [_satu_titik_awal(y, Xs, L, U, s) for s in starts]
    qs1 = np.array([h[1] for h in hasil1])

    # --- Tahap 2: perturbasi dari 5 solusi terbaik ---
    semua = list(hasil1)
    urut = np.argsort(qs1)[:5]
    for idx in urut:
        b_dasar = hasil1[idx][0]
        for _ in range(n_perturb):
            b_coba = b_dasar + rng.normal(size=k) * se * 0.5
            semua.append(_satu_titik_awal(y, Xs, L, U, b_coba))

    qs_all = np.array([h[1] for h in semua])
    i = int(np.argmin(qs_all))
    b_best, q_best = semua[i][0], semua[i][1]

    # --- Tahap 3: pemolesan ---
    if polish:
        try:
            r = minimize(lambda b: objective(b, y, Xs, L, U), b_best,
                         method="Nelder-Mead",
                         options={"maxiter": 3000, "xatol": 1e-10, "fatol": 1e-12})
            if r.fun < q_best:
                b_best, q_best = np.asarray(r.x, dtype=float), float(r.fun)
        except Exception:
            pass

    # --- Kembalikan ke skala asli ---
    b_asli = T @ b_best
    q_asli = objective(b_asli, y, X, L, U)
    return {"b": b_asli, "q": q_asli, "semua_q": qs1,
            "q_tahap1": float(qs1.min()), "q_tahap2": float(q_asli)}