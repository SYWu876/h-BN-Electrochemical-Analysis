"""Rebuild Figure 10 data products for classical / continuous / discrete hBN EIS routes."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "EIS" / "hbn_EIS_1.csv"
OUT = ROOT / "data" / "processed" / "EIS" / "quantum_branches"
OUT.mkdir(parents=True, exist_ok=True)


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


def z_model_fit(freq: np.ndarray, params: np.ndarray) -> np.ndarray:
    l_, rs, log_r1, log_q1, a1, log_q2, a2 = params
    p = np.array([l_, rs, 10**log_r1, 10**log_q1, a1, 10**log_q2, a2], float)
    return z_model_phys(freq, p)


def residuals_fit(params: np.ndarray, f: np.ndarray, zexp: np.ndarray) -> np.ndarray:
    zfit = z_model_fit(f, params)
    return np.concatenate([zfit.real - zexp.real, zfit.imag - zexp.imag])


def sse(zfit: np.ndarray, zexp: np.ndarray) -> float:
    return float(np.sum((zfit.real - zexp.real) ** 2 + (zfit.imag - zexp.imag) ** 2))


def main() -> None:
    df = pd.read_csv(RAW)
    f = df["Frequency (Hz)"].to_numpy(float)
    zre = df["Z' (Ohms)"].to_numpy(float)
    minus_zim = df['-Z" (Ohms)'].to_numpy(float)
    zexp = zre - 1j * minus_zim
    idx = np.argsort(f)

    x0 = [1.83e-7, 1.33, np.log10(34.0), np.log10(3.31e-3), 0.874, np.log10(2.31e-3), 0.915]
    lb = [0.0, 0.0, -3.0, -8.0, 0.20, -8.0, 0.20]
    ub = [1e-3, 100.0, 6.0, 2.0, 1.00, 2.0, 1.00]
    res = least_squares(residuals_fit, x0, bounds=(lb, ub), args=(f, zexp), max_nfev=50000, xtol=1e-12, ftol=1e-12, gtol=1e-12)
    p_class = np.array([res.x[0], res.x[1], 10**res.x[2], 10**res.x[3], res.x[4], 10**res.x[5], res.x[6]], float)

    p_cont = np.array([1.86949039e-07, 1.35875136e+00, 3.38462246e+01, 3.29950223e-03, 8.76651483e-01, 2.31450321e-03, 9.15092829e-01], float)
    p_disc = np.array([1.85226350e-07, 1.27792466e+00, 3.53612945e+01, 3.53824936e-03, 8.60561941e-01, 2.34558071e-03, 9.13259822e-01], float)

    z_class = z_model_phys(f, p_class)
    z_cont = z_model_phys(f, p_cont)
    z_disc = z_model_phys(f, p_disc)

    param_df = pd.DataFrame({
        "Parameter": ["L", "Rs", "R1", "Q1", "alpha1", "Q2", "alpha2"],
        "Classical": p_class,
        "Continuous_branch": p_cont,
        "Discrete_branch": p_disc,
        "Unit": ["H", "Ohm", "Ohm", "S·s^alpha", "-", "S·s^alpha", "-"],
    })
    param_df.to_csv(OUT / "hBN_EIS_quantum_branch_parameters.csv", index=False)

    metrics_df = pd.DataFrame({
        "Route": ["Classical", "Continuous_branch", "Discrete_branch"],
        "SSE_complex": [sse(z_class, zexp), sse(z_cont, zexp), sse(z_disc, zexp)],
    })
    metrics_df.to_csv(OUT / "hBN_EIS_quantum_branch_metrics.csv", index=False)

    zexp_s = zexp[idx]
    overlay_df = pd.DataFrame({
        "Frequency_Hz": f[idx],
        "Zprime_data_Ohm": zexp_s.real,
        "minus_Zdoubleprime_data_Ohm": -zexp_s.imag,
        "Zmag_data_Ohm": np.abs(zexp_s),
        "Phase_data_deg": np.angle(zexp_s, deg=True),
        "Zprime_classical_Ohm": z_class[idx].real,
        "minus_Zdoubleprime_classical_Ohm": -z_class[idx].imag,
        "Zmag_classical_Ohm": np.abs(z_class[idx]),
        "Phase_classical_deg": np.angle(z_class[idx], deg=True),
        "Zprime_continuous_Ohm": z_cont[idx].real,
        "minus_Zdoubleprime_continuous_Ohm": -z_cont[idx].imag,
        "Zmag_continuous_Ohm": np.abs(z_cont[idx]),
        "Phase_continuous_deg": np.angle(z_cont[idx], deg=True),
        "Zprime_discrete_Ohm": z_disc[idx].real,
        "minus_Zdoubleprime_discrete_Ohm": -z_disc[idx].imag,
        "Zmag_discrete_Ohm": np.abs(z_disc[idx]),
        "Phase_discrete_deg": np.angle(z_disc[idx], deg=True),
    })
    overlay_df.to_csv(OUT / "hBN_EIS_quantum_branch_overlays.csv", index=False)

    print("Wrote quantum-branch comparison outputs to", OUT)


if __name__ == "__main__":
    main()
