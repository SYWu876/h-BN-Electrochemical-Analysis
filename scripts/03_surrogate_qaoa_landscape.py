"""Rebuild local surrogate and QAOA landscape tables for the hBN discrete EIS branch."""
from __future__ import annotations

from itertools import product
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "EIS" / "hbn_EIS_1.csv"
OUT = ROOT / "data" / "processed" / "EIS" / "qaoa_landscapes"
OUT.mkdir(parents=True, exist_ok=True)

LOG_IDX = {0, 1, 2, 3, 5}
LB_PHYS = np.array([1e-8, 0.1, 1.0, 1e-4, 0.3, 1e-4, 0.3], float)
UB_PHYS = np.array([1e-6, 10.0, 100.0, 1e-2, 0.99, 1e-2, 0.99], float)


def z_cpe(jw: np.ndarray, q: float, alpha: float) -> np.ndarray:
    return 1.0 / (q * (jw ** alpha))


def z_parallel(z1: np.ndarray | float, z2: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 / z1 + 1.0 / z2)


def z_model_phys(freq: np.ndarray, p: np.ndarray) -> np.ndarray:
    l_, rs, r1, q1, a1, q2, a2 = p
    jw = 1j * 2 * np.pi * freq
    z1 = z_parallel(r1, z_cpe(jw, q1, a1))
    z2 = z_cpe(jw, q2, a2)
    return rs + jw * l_ + z1 + z2


def encode_phys(p: np.ndarray) -> np.ndarray:
    z = np.empty(7, float)
    for i in range(7):
        if i in LOG_IDX:
            z[i] = (np.log10(p[i]) - np.log10(LB_PHYS[i])) / (np.log10(UB_PHYS[i]) - np.log10(LB_PHYS[i]))
        else:
            z[i] = (p[i] - LB_PHYS[i]) / (UB_PHYS[i] - LB_PHYS[i])
    return z


def decode_z(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, float)
    p = np.empty(7, float)
    for i in range(7):
        if i in LOG_IDX:
            p[i] = 10 ** (np.log10(LB_PHYS[i]) + z[i] * (np.log10(UB_PHYS[i]) - np.log10(LB_PHYS[i])))
        else:
            p[i] = LB_PHYS[i] + z[i] * (UB_PHYS[i] - LB_PHYS[i])
    return p


def main() -> None:
    df = pd.read_csv(RAW)
    f = df["Frequency (Hz)"].to_numpy(float)
    zre = df["Z' (Ohms)"].to_numpy(float)
    minus_zim = df['-Z" (Ohms)'].to_numpy(float)
    zexp = zre - 1j * minus_zim

    def z_model_fit(freq: np.ndarray, params: np.ndarray) -> np.ndarray:
        l_, rs, log_r1, log_q1, a1, log_q2, a2 = params
        p = np.array([l_, rs, 10**log_r1, 10**log_q1, a1, 10**log_q2, a2], float)
        return z_model_phys(freq, p)

    def residuals_fit(params: np.ndarray) -> np.ndarray:
        zfit = z_model_fit(f, params)
        return np.concatenate([zfit.real - zexp.real, zfit.imag - zexp.imag])

    x0 = [1.83e-7, 1.33, np.log10(34.0), np.log10(3.31e-3), 0.874, np.log10(2.31e-3), 0.915]
    lb = [0.0, 0.0, -3.0, -8.0, 0.20, -8.0, 0.20]
    ub = [1e-3, 100.0, 6.0, 2.0, 1.00, 2.0, 1.00]
    res = least_squares(residuals_fit, x0, bounds=(lb, ub), max_nfev=50000, xtol=1e-12, ftol=1e-12, gtol=1e-12)
    p_class = np.array([res.x[0], res.x[1], 10**res.x[2], 10**res.x[3], res.x[4], 10**res.x[5], res.x[6]], float)
    z0 = encode_phys(p_class)

    def sse_phys(p: np.ndarray) -> float:
        zfit = z_model_phys(f, p)
        d = zfit - zexp
        return float(np.sum(d.real**2 + d.imag**2))

    def sse_z(z: np.ndarray) -> float:
        return sse_phys(decode_z(np.clip(np.asarray(z, float), 0.0, 1.0)))

    # Local quadratic surrogate
    n = 7
    h = 0.015
    s0 = sse_z(z0)
    g = np.zeros(n)
    H = np.zeros((n, n))
    for i in range(n):
        e = np.zeros(n); e[i] = 1
        sp = sse_z(z0 + h * e)
        sm = sse_z(z0 - h * e)
        g[i] = (sp - sm) / (2 * h)
        H[i, i] = (sp - 2 * s0 + sm) / (h**2)
    for i in range(n):
        for j in range(i + 1, n):
            ei = np.zeros(n); ej = np.zeros(n)
            ei[i] = 1; ej[j] = 1
            spp = sse_z(z0 + h * ei + h * ej)
            spm = sse_z(z0 + h * ei - h * ej)
            smp = sse_z(z0 - h * ei + h * ej)
            smm = sse_z(z0 - h * ei - h * ej)
            hij = (spp - spm - smp + smm) / (4 * h * h)
            H[i, j] = hij
            H[j, i] = hij

    def sse_surrogate(z: np.ndarray) -> float:
        dz = np.asarray(z, float) - z0
        return float(s0 + g @ dz + 0.5 * dz @ H @ dz)

    # Panel-like data products
    i1, i2 = 3, 4
    eta = 0.08
    xg = np.linspace(np.clip(z0[i1] - eta, 0, 1), np.clip(z0[i1] + eta, 0, 1), 180)
    yg = np.linspace(np.clip(z0[i2] - eta, 0, 1), np.clip(z0[i2] + eta, 0, 1), 180)
    XX, YY = np.meshgrid(xg, yg)
    ZZ = np.zeros_like(XX)
    for r in range(XX.shape[0]):
        for c in range(XX.shape[1]):
            z = z0.copy()
            z[i1] = XX[r, c]
            z[i2] = YY[r, c]
            ZZ[r, c] = sse_surrogate(z)
    pd.DataFrame({"z_Q1": XX.ravel(), "z_alpha1": YY.ravel(), "surrogate_SSE": ZZ.ravel()}).to_csv(OUT / "Figure_12a_hBN_surrogate_slice.csv", index=False)

    bitstrings = np.array(list(product([0, 1], repeat=6)), dtype=int)
    weights = np.array([4, 2, 1])

    def bits_to_z_pair(bits: np.ndarray) -> np.ndarray:
        l1 = bits[:3] @ weights
        l2 = bits[3:] @ weights
        return np.array([
            np.clip(z0[i1] + eta * (2 * l1 / 7 - 1), 0.0, 1.0),
            np.clip(z0[i2] + eta * (2 * l2 / 7 - 1), 0.0, 1.0),
        ])

    energies = []
    for b in bitstrings:
        z = z0.copy()
        pair = bits_to_z_pair(b)
        z[i1], z[i2] = pair
        energies.append(sse_surrogate(z))
    energies = np.array(energies)

    pairs = [(i, j) for i in range(6) for j in range(i + 1, 6)]
    Phi = np.ones((len(bitstrings), 1 + 6 + len(pairs)))
    Phi[:, 1:7] = bitstrings
    for k, (i, j) in enumerate(pairs):
        Phi[:, 7 + k] = bitstrings[:, i] * bitstrings[:, j]
    coef, *_ = np.linalg.lstsq(Phi, energies, rcond=None)
    cost_diag = Phi @ coef

    N = 6
    dim = 2 ** N
    plus = np.ones(dim, dtype=complex) / np.sqrt(dim)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    I = np.eye(2, dtype=complex)

    def mixer_unitary(beta: float) -> np.ndarray:
        rx = np.cos(beta) * I - 1j * np.sin(beta) * X
        U = rx
        for _ in range(N - 1):
            U = np.kron(U, rx)
        return U

    gamma_coarse = np.linspace(0, np.pi, 41)
    beta_coarse = np.linspace(0, np.pi / 2, 41)
    EEc = np.zeros((len(beta_coarse), len(gamma_coarse)))
    for i, b in enumerate(beta_coarse):
        Ub = mixer_unitary(b)
        for j, gma in enumerate(gamma_coarse):
            Ug = np.exp(-1j * gma * cost_diag)
            psi = Ub @ (Ug * plus)
            EEc[i, j] = float(np.real(np.vdot(psi, cost_diag * psi)))
    jc = np.unravel_index(np.argmin(EEc), EEc.shape)
    beta_star_c = beta_coarse[jc[0]]
    gamma_star_c = gamma_coarse[jc[1]]
    pd.DataFrame({
        "gamma": np.repeat(gamma_coarse[None, :], len(beta_coarse), axis=0).ravel(),
        "beta": np.repeat(beta_coarse[:, None], len(gamma_coarse), axis=1).ravel(),
        "expected_energy": EEc.ravel(),
    }).to_csv(OUT / "Figure_12b_hBN_qaoa_coarse_landscape.csv", index=False)

    gamma_ref = np.linspace(max(0, gamma_star_c - np.pi / 8), min(np.pi, gamma_star_c + np.pi / 8), 61)
    beta_ref = np.linspace(max(0, beta_star_c - np.pi / 10), min(np.pi / 2, beta_star_c + np.pi / 10), 61)
    EEr = np.zeros((len(beta_ref), len(gamma_ref)))
    for i, b in enumerate(beta_ref):
        Ub = mixer_unitary(b)
        for j, gma in enumerate(gamma_ref):
            Ug = np.exp(-1j * gma * cost_diag)
            psi = Ub @ (Ug * plus)
            EEr[i, j] = float(np.real(np.vdot(psi, cost_diag * psi)))
    pd.DataFrame({
        "gamma": np.repeat(gamma_ref[None, :], len(beta_ref), axis=0).ravel(),
        "beta": np.repeat(beta_ref[:, None], len(gamma_ref), axis=1).ravel(),
        "expected_energy": EEr.ravel(),
    }).to_csv(OUT / "Figure_12c_hBN_qaoa_refined_landscape.csv", index=False)

    rng = np.random.default_rng(7)
    samples = np.clip(z0 + rng.uniform(-eta, eta, size=(200, 7)), 0.0, 1.0)
    true_sse = np.array([sse_z(z) for z in samples])
    pred_sse = np.array([sse_surrogate(z) for z in samples])
    pd.DataFrame({"surrogate_SSE": pred_sse, "true_SSE": true_sse}).to_csv(OUT / "Figure_12d_hBN_surrogate_fidelity.csv", index=False)

    summary = pd.DataFrame({
        "Item": ["finite_difference_step_h", "trust_region_eta", "bits_per_parameter", "levels_per_parameter", "shots", "gamma_coarse_opt", "beta_coarse_opt"],
        "Value": [0.015, 0.08, 3, 8, 4096, gamma_star_c, beta_star_c],
    })
    summary.to_csv(OUT / "hBN_qaoa_summary.csv", index=False)

    print("Wrote surrogate / QAOA tables to", OUT)


if __name__ == "__main__":
    main()
