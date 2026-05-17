from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from cv_utils import SCAN_RATES, extract_all_curves, load_cv_dataframe, processed_dir, set_nc_style, style_axes


def _ensure_increasing(E_V: np.ndarray, I_mA: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(E_V)
    E_sorted = E_V[order]
    I_sorted = I_mA[order]
    unique_E, inverse = np.unique(E_sorted, return_inverse=True)
    if unique_E.size == E_sorted.size:
        return E_sorted, I_sorted

    I_unique = np.empty(unique_E.size, dtype=float)
    for idx in range(unique_E.size):
        I_unique[idx] = I_sorted[inverse == idx].mean()
    return unique_E, I_unique


def split_cv_branches(curve: dict[str, np.ndarray | float]) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    E_V = np.asarray(curve["E_V"], dtype=float)
    I_mA = np.asarray(curve["I_mA"], dtype=float)
    turn_idx = int(np.argmin(E_V))

    cathodic_E, cathodic_I = _ensure_increasing(E_V[: turn_idx + 1], I_mA[: turn_idx + 1])
    anodic_E, anodic_I = _ensure_increasing(E_V[turn_idx:], I_mA[turn_idx:])
    return {
        "cathodic": (cathodic_E, cathodic_I),
        "anodic": (anodic_E, anodic_I),
    }


def common_potential_grid(curves, branch_name: str, n_points: int = 700) -> np.ndarray:
    branches = [split_cv_branches(curve)[branch_name] for curve in curves]
    E_min = max(branch[0].min() for branch in branches)
    E_max = min(branch[0].max() for branch in branches)
    if not E_min < E_max:
        raise ValueError(f"No overlapping potential window for {branch_name} CV branches")
    return np.linspace(E_min, E_max, n_points)


def interpolate_branch_matrix(curves, branch_name: str, E_grid: np.ndarray) -> np.ndarray:
    rows = []
    for curve in curves:
        E_branch, I_branch = split_cv_branches(curve)[branch_name]
        rows.append(np.interp(E_grid, E_branch, I_branch))
    return np.vstack(rows)


def scan_rate_edges(scan_rates: np.ndarray = SCAN_RATES) -> np.ndarray:
    rates = np.asarray(scan_rates, dtype=float)
    if rates.ndim != 1 or rates.size < 2:
        raise ValueError("At least two scan rates are required for heatmap edges")
    midpoints = (rates[:-1] + rates[1:]) / 2.0
    first = rates[0] - (midpoints[0] - rates[0])
    last = rates[-1] + (rates[-1] - midpoints[-1])
    return np.concatenate([[first], midpoints, [last]])


def grid_edges(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size < 2:
        raise ValueError("At least two grid values are required for heatmap edges")
    midpoints = (values[:-1] + values[1:]) / 2.0
    first = values[0] - (midpoints[0] - values[0])
    last = values[-1] + (values[-1] - midpoints[-1])
    return np.concatenate([[first], midpoints, [last]])


def branch_rms_intensity_matrix(
    cathodic_matrix: np.ndarray,
    anodic_matrix: np.ndarray,
) -> np.ndarray:
    if cathodic_matrix.shape != anodic_matrix.shape:
        raise ValueError("Cathodic and anodic matrices must have the same shape")
    return np.sqrt((cathodic_matrix**2 + anodic_matrix**2) / 2.0)


def main() -> None:
    set_nc_style()
    curves = extract_all_curves(load_cv_dataframe())
    out = processed_dir()

    E_grid = common_potential_grid(curves, "anodic")
    cathodic_matrix = interpolate_branch_matrix(curves, "cathodic", E_grid)
    anodic_matrix = interpolate_branch_matrix(curves, "anodic", E_grid)
    heatmap_matrix = branch_rms_intensity_matrix(cathodic_matrix, anodic_matrix)
    x_edges = grid_edges(E_grid)
    y_edges = scan_rate_edges()

    fig, ax = plt.subplots(figsize=(6.4, 4.8), dpi=300)
    mesh = ax.pcolormesh(x_edges, y_edges, heatmap_matrix, cmap="viridis")
    ax.set_xlabel("Potential E(V)", fontsize=16)
    ax.set_ylabel(r"Scan rate $\nu$ (mV/s)", fontsize=16)
    style_axes(ax)

    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(r"$|I(E)|$ (mA)", fontsize=15)
    fig.tight_layout()
    fig.savefig(out / "CV_current_heatmap.png", bbox_inches="tight")
    fig.savefig(out / "CV_current_heatmap.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
