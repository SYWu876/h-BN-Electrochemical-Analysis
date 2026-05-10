
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import ListedColormap

SCAN_RATES = np.array([5, 10, 20, 30, 40, 50], dtype=float)  # mV s^-1

def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]

def raw_cv_path() -> Path:
    return repo_root() / "data" / "raw" / "CV" / "CV_total.csv"

def processed_dir() -> Path:
    out = repo_root() / "data" / "processed" / "CV"
    out.mkdir(parents=True, exist_ok=True)
    return out

def load_cv_dataframe() -> pd.DataFrame:
    return pd.read_csv(raw_cv_path())

def extract_all_curves(df: pd.DataFrame | None = None):
    if df is None:
        df = load_cv_dataframe()
    curves = []
    for i, rate in enumerate(SCAN_RATES):
        sub = df.iloc[:, [2*i, 2*i + 1]].dropna()
        E = sub.iloc[:, 0].to_numpy(dtype=float)
        I_A = sub.iloc[:, 1].to_numpy(dtype=float)
        curves.append({"scan_rate_mV_s": float(rate), "E_V": E, "I_A": I_A, "I_mA": I_A * 1000.0})
    return curves

def extract_anodic_branches(df: pd.DataFrame | None = None):
    curves = extract_all_curves(df)
    branches = []
    for c in curves:
        E = c["E_V"]
        I_mA = c["I_mA"]
        turn_idx = int(np.argmin(E))
        branches.append({
            "scan_rate_mV_s": c["scan_rate_mV_s"],
            "E_V": E[turn_idx + 1:],
            "I_mA": I_mA[turn_idx + 1:]
        })
    return branches

def get_common_anodic_grid(branches=None, n_points: int = 700):
    if branches is None:
        branches = extract_anodic_branches()
    E_min = max(b["E_V"][0] for b in branches)
    E_max = min(b["E_V"][-1] for b in branches)
    E_grid = np.linspace(E_min, E_max, n_points)
    I_grid = np.vstack([
        np.interp(E_grid, b["E_V"], b["I_mA"])
        for b in branches
    ])
    return E_grid, I_grid

def peak_summary(df: pd.DataFrame | None = None) -> pd.DataFrame:
    curves = extract_all_curves(df)
    rows = []
    for c in curves:
        idx = int(np.argmax(c["I_mA"]))
        rows.append({
            "scan_rate_mV_s": c["scan_rate_mV_s"],
            "sqrt_scan_rate_sqrt_mV_s": np.sqrt(c["scan_rate_mV_s"]),
            "peak_potential_V": c["E_V"][idx],
            "peak_current_mA": c["I_mA"][idx]
        })
    return pd.DataFrame(rows)

def linear_fit(x: np.ndarray, y: np.ndarray):
    coef = np.polyfit(x, y, 1)
    y_fit = np.polyval(coef, x)
    ss_res = float(np.sum((y - y_fit) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = np.nan if np.isclose(ss_tot, 0.0) else 1.0 - ss_res / ss_tot
    return coef, y_fit, r2

def b_value_profile(E_grid: np.ndarray | None = None, I_grid: np.ndarray | None = None):
    if E_grid is None or I_grid is None:
        branches = extract_anodic_branches()
        E_grid, I_grid = get_common_anodic_grid(branches)
    valid = np.all(I_grid > 0, axis=0)
    E_valid = E_grid[valid]
    I_valid = I_grid[:, valid]
    log_v = np.log(SCAN_RATES)
    b_vals = np.empty(E_valid.size)
    r2_vals = np.empty(E_valid.size)
    for j in range(E_valid.size):
        y = np.log(I_valid[:, j])
        coef, y_fit, r2 = linear_fit(log_v, y)
        b_vals[j] = coef[0]
        r2_vals[j] = r2
    window = 31 if E_valid.size >= 31 else max(5, (E_valid.size // 2) * 2 + 1)
    b_smooth = savgol_filter(b_vals, window_length=window, polyorder=3)
    return pd.DataFrame({
        "Potential_E_V": E_valid,
        "b_raw": b_vals,
        "b_smooth": b_smooth,
        "fit_R2": r2_vals
    })

def dunn_partition(E_grid: np.ndarray | None = None, I_grid: np.ndarray | None = None):
    if E_grid is None or I_grid is None:
        branches = extract_anodic_branches()
        E_grid, I_grid = get_common_anodic_grid(branches)
    sqrt_v = np.sqrt(SCAN_RATES)
    k1 = np.empty(E_grid.size)
    k2 = np.empty(E_grid.size)
    fit_r2 = np.empty(E_grid.size)
    for j in range(E_grid.size):
        y = I_grid[:, j] / sqrt_v
        coef, y_fit, r2 = linear_fit(sqrt_v, y)
        k1[j] = coef[0]
        k2[j] = coef[1]
        fit_r2[j] = r2
    f_cap = []
    for idx, v in enumerate(SCAN_RATES):
        i_cap = k1 * v
        i_tot = I_grid[idx]
        f = np.trapezoid(np.abs(i_cap), E_grid) / np.trapezoid(np.abs(i_tot), E_grid)
        f_cap.append(f)
    frac = pd.DataFrame({
        "scan_rate_mV_s": SCAN_RATES,
        "fractional_capacitive_fcap": np.array(f_cap)
    })
    prof = pd.DataFrame({
        "Potential_E_V": E_grid,
        "k1_mA_per_mV_per_s": k1,
        "k2_mA_per_sqrt_mV_per_s": k2,
        "local_fit_R2": fit_r2
    })
    return prof, frac

def robust_z(x: np.ndarray) -> np.ndarray:
    med = np.median(x)
    mad = np.median(np.abs(x - med)) + 1e-12
    return (x - med) / (1.4826 * mad)

def ising_segmentation(E_grid: np.ndarray | None = None, I_grid: np.ndarray | None = None):
    if E_grid is None or I_grid is None:
        branches = extract_anodic_branches()
        E_grid, I_grid = get_common_anodic_grid(branches, n_points=500)
    bdf = b_value_profile(E_grid, I_grid)
    b_on_grid = np.interp(E_grid, bdf["Potential_E_V"].to_numpy(), bdf["b_smooth"].to_numpy())
    dunn_df, frac = dunn_partition(E_grid, I_grid)
    k1 = dunn_df["k1_mA_per_mV_per_s"].to_numpy()
    k2 = dunn_df["k2_mA_per_sqrt_mV_per_s"].to_numpy()
    phi_local = np.mean([
        np.abs(k1 * v) / (np.abs(k1 * v) + np.abs(k2 * np.sqrt(v)) + 1e-12)
        for v in SCAN_RATES
    ], axis=0)
    h_field = 0.65 * robust_z(phi_local - 0.5) + 0.35 * robust_z(b_on_grid - 0.75)
    h_field = savgol_filter(h_field, 31, 3)
    J = 0.35
    n = E_grid.size
    unary = np.vstack([+h_field, -h_field]).T  # 0=DD, 1=CD
    cost = np.zeros((n, 2))
    back = np.zeros((n, 2), dtype=int)
    cost[0] = unary[0]
    for i in range(1, n):
        for s in range(2):
            vals = cost[i - 1] + J * (np.arange(2) != s)
            back[i, s] = int(np.argmin(vals))
            cost[i, s] = unary[i, s] + vals[back[i, s]]
    states = np.zeros(n, dtype=int)
    states[-1] = int(np.argmin(cost[-1]))
    for i in range(n - 2, -1, -1):
        states[i] = back[i + 1, states[i + 1]]
    return pd.DataFrame({
        "Potential_V": E_grid,
        "local_field_h": h_field,
        "b_value": b_on_grid,
        "phi_local_avg": phi_local,
        "Ising_state": states
    })

def qkpca_coordinates(branches=None):
    if branches is None:
        branches = extract_anodic_branches()
    E_min = max(b["E_V"][0] for b in branches)
    E_max = min(b["E_V"][-1] for b in branches)
    E_grid = np.linspace(E_min, E_max, 80)
    X = np.vstack([np.interp(E_grid, b["E_V"], b["I_mA"]) for b in branches])
    X_scaled = (X - X.min()) / (X.max() - X.min() + 1e-12) * np.pi
    m = X_scaled.shape[0]
    K = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            vals = np.cos((X_scaled[i] - X_scaled[j]) / 2.0) ** 2
            K[i, j] = np.exp(np.mean(np.log(vals + 1e-12)))
    one = np.ones((m, m)) / m
    Kc = K - one @ K - K @ one + one @ K @ one
    eigvals, eigvecs = np.linalg.eigh(Kc)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    coords = eigvecs[:, :2] * np.sqrt(np.maximum(eigvals[:2], 0))
    return pd.DataFrame({
        "scan_rate_mV_s": SCAN_RATES,
        "QKPCA1": coords[:, 0],
        "QKPCA2": coords[:, 1]
    })

def set_nc_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 14,
        "axes.linewidth": 1.4,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.minor.width": 0.8,
        "ytick.minor.width": 0.8,
    })

def style_axes(ax, minor=True):
    ax.tick_params(direction="in", top=True, right=True, length=5.5, width=1.2, pad=6)
    if minor:
        ax.tick_params(which="minor", direction="in", top=True, right=True, length=3.0, width=0.8)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

def ising_cmap():
    return ListedColormap(["white", "#e69f00"])
