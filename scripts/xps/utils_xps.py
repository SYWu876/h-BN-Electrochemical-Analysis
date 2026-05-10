
from __future__ import annotations
import os
import shutil
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

TABLE_X_ROWS = [
    ("B1s", "B-def", 189.87, 0.2315, "lower-BE defect-related B environment, possibly B-rich / under-coordinated B"),
    ("B1s", "B–N",   190.37, 0.5727, "main lattice B–N component of h-BN"),
    ("B1s", "B-ox",  191.09, 0.1958, "high-BE defect / oxidized boron shoulder, such as B–Ox-related surface species"),
    ("N1s", "N-def", 397.70, 0.1785, "defect-rich N / perturbed nitrogen environment"),
    ("N1s", "N–B",   398.00, 0.6855, "main lattice N–B component of h-BN"),
    ("N1s", "N-ox",  398.85, 0.1360, "surface / oxidized N shoulder"),
    ("C1s", "C–C",   284.59, 0.7161, "C–C / C–H"),
    ("C1s", "C–O",   285.50, 0.0967, "C–N / C–O"),
    ("C1s", "C=O",   286.35, 0.0923, "C=O-related"),
    ("C1s", "O–C=O", 288.43, 0.0949, "O–C=O"),
    ("O1s", "O-lat", 530.37, 0.0041, "lattice / B–O-related O"),
    ("O1s", "O-ads", 532.60, 0.7565, "surface hydroxyl / adsorbed O"),
    ("O1s", "H2O/O", 533.34, 0.2394, "adsorbed H2O / weak oxygen species"),
]

REGION_SPECS = {
    "B1s": dict(sheet="B1s Scan A", window=(188.0, 193.5), plot_xlim=(193.5, 188.0),
                centers=[189.87, 190.37, 191.09], labels=["B-def", "B–N", "B-ox"], main=190.37,
                amp0=[2500, 7000, 1800], sig0=[0.35, 0.40, 0.55], eta0=[0.35, 0.20, 0.45], ampmax=[20000, 30000, 15000], sigmax=[1.6, 1.6, 2.0]),
    "N1s": dict(sheet="N1s Scan A", window=(395.5, 401.5), plot_xlim=(401.5, 395.5),
                centers=[397.70, 398.00, 398.85], labels=["N-def", "N–B", "N-ox"], main=398.00,
                amp0=[3000, 26000, 3500], sig0=[0.45, 0.38, 0.55], eta0=[0.25, 0.20, 0.40], ampmax=[15000, 50000, 15000], sigmax=[1.8, 1.5, 2.0]),
    "C1s": dict(sheet="C1s Scan A", window=(280.0, 297.2), fit_window=(282.0, 292.0), plot_xlim=(297.2, 280.0),
                centers=[284.59, 285.50, 286.35, 288.43], labels=["C–C", "C–O", "C=O", "O–C=O"], main=284.59,
                amp0=[850, 140, 100, 35], sig0=[0.35, 0.35, 0.50, 1.10], eta0=[0.20, 0.40, 0.40, 0.50], ampmax=[2500, 1200, 1200, 600], sigmax=[1.4, 1.4, 1.8, 3.0]),
    "O1s": dict(sheet="O1s Scan A", window=(526.0, 540.0), fit_window=(528.0, 535.5), plot_xlim=(540.0, 526.0),
                centers=[530.37, 532.60, 533.34], labels=["O-lat", "O-ads", "H2O/O"], main=532.60,
                amp0=[120, 2400, 500], sig0=[0.25, 0.60, 0.65], eta0=[0.20, 0.20, 0.35], ampmax=[2000, 12000, 4000], sigmax=[1.4, 2.0, 2.2]),
}

MARKER_MAP = {"B1s": "o", "N1s": "s", "C1s": "^", "O1s": "D"}

def ensure_xlsx(raw_file: str | Path) -> Path:
    raw_file = Path(raw_file)
    if raw_file.suffix.lower() == ".xlsx":
        return raw_file
    xlsx = raw_file.with_suffix(".xlsx")
    if xlsx.exists():
        return xlsx
    try:
        subprocess.run(
            ["libreoffice", "--headless", "--convert-to", "xlsx", "--outdir", str(raw_file.parent), str(raw_file)],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
    except Exception as exc:
        raise RuntimeError(
            "Failed to convert legacy .xls to .xlsx. Please install LibreOffice or provide an .xlsx file."
        ) from exc
    if not xlsx.exists():
        raise FileNotFoundError(f"Expected converted file not found: {xlsx}")
    return xlsx

def load_region(raw_file: str | Path, sheet_name: str) -> pd.DataFrame:
    raw_file = Path(raw_file)
    try:
        raw = pd.read_excel(raw_file, sheet_name=sheet_name, header=None)
    except Exception:
        xlsx = ensure_xlsx(raw_file)
        raw = pd.read_excel(xlsx, sheet_name=sheet_name, engine="openpyxl", header=None)
    df = raw.iloc[16:, [0, 2, 4]].copy()
    df.columns = ["BE", "counts", "background"]
    df = df.dropna()
    df = df[pd.to_numeric(df["BE"], errors="coerce").notna()]
    df = df.astype(float).reset_index(drop=True)
    df = df.iloc[1:].sort_values("BE").reset_index(drop=True)
    df["net"] = df["counts"] - df["background"]
    return df

def pseudo_voigt(x, amp, sigma, eta, center):
    gauss = np.exp(-((x - center) ** 2) / (2 * sigma ** 2))
    lorentz = 1.0 / (1.0 + ((x - center) / sigma) ** 2)
    return amp * (eta * lorentz + (1.0 - eta) * gauss)

def numerical_fwhm(x, y):
    ymax = np.max(y)
    if ymax <= 0:
        return np.nan
    half = ymax / 2.0
    idx_peak = int(np.argmax(y))
    left_idx = idx_peak
    while left_idx > 0 and y[left_idx] > half:
        left_idx -= 1
    right_idx = idx_peak
    while right_idx < len(y) - 1 and y[right_idx] > half:
        right_idx += 1
    return abs(x[right_idx] - x[left_idx])

def fit_fixed_centers(x, y, centers, amp0, sig0, eta0, ampmax, sigmax):
    n = len(centers)
    def model(x, *params):
        total = np.zeros_like(x, dtype=float)
        for i, c in enumerate(centers):
            total += pseudo_voigt(x, params[3*i], params[3*i+1], params[3*i+2], c)
        return total
    p0, lower, upper = [], [], []
    for i in range(n):
        p0 += [amp0[i], sig0[i], eta0[i]]
        lower += [0.0, 0.05, 0.0]
        upper += [ampmax[i], sigmax[i], 1.0]
    params, _ = curve_fit(model, x, y, p0=p0, bounds=(lower, upper), maxfev=800000)
    comps = []
    for i, c in enumerate(centers):
        amp = params[3*i]
        sigma = params[3*i+1]
        eta = params[3*i+2]
        yy = pseudo_voigt(x, amp, sigma, eta, c)
        comps.append({"center": c, "amp": amp, "sigma": sigma, "eta": eta, "y": yy})
    return comps

def table_x_dataframe():
    return pd.DataFrame(TABLE_X_ROWS, columns=["Region", "Label", "Center_eV", "Area_fraction", "Tentative_assignment"])

def build_descriptor_matrix(raw_file: str | Path):
    table_df = table_x_dataframe()
    records = []
    for region, spec in REGION_SPECS.items():
        df = load_region(raw_file, spec["sheet"])
        x = df["BE"].to_numpy()
        y = df["net"].to_numpy()
        fit_window = spec.get("fit_window", spec["window"])
        mask = (x >= fit_window[0]) & (x <= fit_window[1])
        xf, yf = x[mask], y[mask]
        comps = fit_fixed_centers(xf, yf, spec["centers"], spec["amp0"], spec["sig0"], spec["eta0"], spec["ampmax"], spec["sigmax"])
        for label, center, comp in zip(spec["labels"], spec["centers"], comps):
            area_fraction = float(table_df[(table_df.Region == region) & (table_df.Label == label)]["Area_fraction"].iloc[0])
            records.append({
                "Region": region,
                "Label": label,
                "Center_eV": center,
                "FWHM_eV": numerical_fwhm(xf, comp["y"]),
                "Area_fraction": area_fraction,
                "Eta": comp["eta"],
                "DeltaE_eV": center - spec["main"],
            })
    return pd.DataFrame(records)

def ensure_dir(path: str | Path):
    Path(path).mkdir(parents=True, exist_ok=True)
