"""Classical hBN EIS fit using L-Rs-(R1||Q1)-Q2.

Exports:
- hBN_EIS_final_fitting_curve.csv
- hBN_EIS_final_fit_parameters.csv
- hBN_EIS_residuals_final_model.csv
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "EIS" / "hbn_EIS_1.csv"
OUT = ROOT / "data" / "processed" / "EIS" / "classical_fit"
OUT.mkdir(parents=True, exist_ok=True)


def z_cpe(jw: np.ndarray, q: float, alpha: float) -> np.ndarray:
    return 1.0 / (q * (jw ** alpha))


def z_parallel(z1: np.ndarray | float, z2: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 / z1 + 1.0 / z2)


def z_model(freq: np.ndarray, params: np.ndarray) -> np.ndarray:
    l_, rs, log_r1, log_q1, a1, log_q2, a2 = params
    jw = 1j * 2 * np.pi * freq
    r1 = 10.0 ** log_r1
    q1 = 10.0 ** log_q1
    q2 = 10.0 ** log_q2
    z1 = z_parallel(r1, z_cpe(jw, q1, a1))
    z2 = z_cpe(jw, q2, a2)
    return rs + jw * l_ + z1 + z2


def residuals(params: np.ndarray, freq: np.ndarray, zexp: np.ndarray) -> np.ndarray:
    zfit = z_model(freq, params)
    return np.concatenate([zfit.real - zexp.real, zfit.imag - zexp.imag])


def main() -> None:
    df = pd.read_csv(RAW)
    f = df["Frequency (Hz)"].to_numpy(float)
    zre = df["Z' (Ohms)"].to_numpy(float)
    minus_zim = df['-Z" (Ohms)'].to_numpy(float)
    zexp = zre - 1j * minus_zim

    x0 = [1.83e-7, 1.33, np.log10(34.0), np.log10(3.31e-3), 0.874, np.log10(2.31e-3), 0.915]
    lb = [0.0, 0.0, -3.0, -8.0, 0.20, -8.0, 0.20]
    ub = [1e-3, 100.0, 6.0, 2.0, 1.00, 2.0, 1.00]

    res = least_squares(
        residuals, x0, bounds=(lb, ub), args=(f, zexp),
        max_nfev=50000, xtol=1e-12, ftol=1e-12, gtol=1e-12
    )

    zfit = z_model(f, res.x)
    idx = np.argsort(f)

    fit_df = pd.DataFrame({
        "Frequency_Hz": f[idx],
        "Zprime_data_Ohm": zexp.real[idx],
        "minus_Zdoubleprime_data_Ohm": (-zexp.imag)[idx],
        "Zmag_data_Ohm": np.abs(zexp)[idx],
        "Phase_data_deg": np.angle(zexp, deg=True)[idx],
        "Zprime_fit_Ohm": zfit.real[idx],
        "minus_Zdoubleprime_fit_Ohm": (-zfit.imag)[idx],
        "Zmag_fit_Ohm": np.abs(zfit)[idx],
        "Phase_fit_deg": np.angle(zfit, deg=True)[idx],
    })
    fit_df.to_csv(OUT / "hBN_EIS_final_fitting_curve.csv", index=False)

    l_, rs, log_r1, log_q1, a1, log_q2, a2 = res.x
    param_df = pd.DataFrame({
        "Parameter": ["L", "Rs", "R1", "Q1", "alpha1", "Q2", "alpha2"],
        "Value": [l_, rs, 10**log_r1, 10**log_q1, a1, 10**log_q2, a2],
        "Unit": ["H", "Ohm", "Ohm", "S·s^alpha", "-", "S·s^alpha", "-"],
    })
    param_df.to_csv(OUT / "hBN_EIS_final_fit_parameters.csv", index=False)

    zfit_s = zfit[idx]
    zexp_s = zexp[idx]
    norm = np.abs(zexp_s)
    residual_df = pd.DataFrame({
        "Frequency_Hz": f[idx],
        "Residual_ReZ_data_minus_fit_Ohm": zexp_s.real - zfit_s.real,
        "Residual_ImZ_data_minus_fit_Ohm": zexp_s.imag - zfit_s.imag,
        "Residual_ReZ_over_absZ": (zexp_s.real - zfit_s.real) / norm,
        "Residual_ImZ_over_absZ": (zexp_s.imag - zfit_s.imag) / norm,
    })
    residual_df.to_csv(OUT / "hBN_EIS_residuals_final_model.csv", index=False)

    print("Wrote classical fit outputs to", OUT)


if __name__ == "__main__":
    main()
