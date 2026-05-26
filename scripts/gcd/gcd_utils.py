
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.signal import savgol_filter
from scipy.optimize import least_squares

JS = [1, 2, 3, 4, 5]
POLARIZATION_TAU_S = 9.0
POLARIZATION_AMPLITUDE_FRACTION_OF_PEAK = 0.35
SLOPE_IR_EXCLUSION_S = 0.15
SLOPE_MIN_VOLTAGE_SPAN_FRACTION = 0.18
SLOPE_MAX_VOLTAGE_SPAN_FRACTION = 0.45
SLOPE_MIN_POINT_FRACTION = 0.04
SLOPE_MAX_POINT_FRACTION = 0.85
SLOPE_CANDIDATE_LENGTHS = 36
MASKED_DOMAIN_RS_ANCHORS = {
    1.0: 0.00986700312913102,
    2.0: 0.009732716467688456,
    3.0: 0.009627444607500105,
    4.0: 0.009465826381994863,
    5.0: 0.009339393602113992,
}

def setup_matplotlib() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 13,
        "axes.linewidth": 1.3,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.minor.width": 0.8,
        "ytick.minor.width": 0.8,
    })

def style_axes(ax) -> None:
    ax.tick_params(direction="in", top=True, right=True, length=5, width=1.1, pad=4)
    ax.tick_params(which="minor", direction="in", top=True, right=True, length=2.8, width=0.8)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

def load_gcd_pair(csv_path: str | Path, j_value: int) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(csv_path)
    idx = 2 * (j_value - 1)
    sub = df.iloc[:, idx:idx+2].dropna()
    t = sub.iloc[:, 0].to_numpy(dtype=float)
    v = sub.iloc[:, 1].to_numpy(dtype=float)
    return t, v

def extract_discharge_branch(t_full: np.ndarray, v_full: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    imax = int(np.argmax(v_full))
    t_dis = t_full[imax:] - t_full[imax]
    v_dis = v_full[imax:]
    return t_dis, v_dis

def smooth_voltage(v_dis: np.ndarray) -> np.ndarray:
    window = 51 if len(v_dis) >= 51 else max(5, (len(v_dis) // 2) * 2 + 1)
    if window >= len(v_dis):
        window = len(v_dis) - 1 if len(v_dis) % 2 == 0 else len(v_dis)
        if window < 5:
            window = max(3, (len(v_dis) // 2) * 2 + 1)
    return savgol_filter(v_dis, window_length=window, polyorder=2)

def minmax(x: np.ndarray) -> np.ndarray:
    xmin, xmax = np.min(x), np.max(x)
    if np.isclose(xmax, xmin):
        return np.zeros_like(x)
    return (x - xmin) / (xmax - xmin)

def compute_diagnostics(t_dis: np.ndarray, v_dis: np.ndarray, half_window_s: float = 0.25) -> pd.DataFrame:
    v_smooth = smooth_voltage(v_dis)
    dvdt = np.gradient(v_smooth, t_dis)
    d2vdt2 = np.gradient(dvdt, t_dis)
    curv = np.abs(d2vdt2)

    residual = np.zeros_like(v_smooth)
    for i, tc in enumerate(t_dis):
        mask = (t_dis >= tc - half_window_s) & (t_dis <= tc + half_window_s)
        x = t_dis[mask]
        y = v_smooth[mask]
        if len(x) >= 3:
            p = np.polyfit(x, y, 1)
            residual[i] = np.sqrt(np.mean((y - np.polyval(p, x)) ** 2))

    curv_n = minmax(curv)
    res_n = minmax(residual)
    a_i = 0.6 * curv_n + 0.4 * res_n

    return pd.DataFrame({
        "time_s": t_dis,
        "raw_V": v_dis,
        "smoothed_V": v_smooth,
        "dVdt_V_per_s": dvdt,
        "abs_d2Vdt2_V_per_s2": curv,
        "local_residual_V": residual,
        "curvature_norm": curv_n,
        "residual_norm": res_n,
        "stability_cost_ai": a_i,
    })

def select_stable_window(diag: pd.DataFrame, j_value: int, ir_cut: float = 0.15) -> pd.DataFrame:
    t = diag["time_s"].to_numpy(dtype=float)
    a_i = diag["stability_cost_ai"].to_numpy(dtype=float)
    mask = np.zeros_like(t, dtype=int)

    if j_value in (1, 2, 3):
        target_len_s = 2.5
        dt_med = np.median(np.diff(t))
        m = max(10, int(round(target_len_s / dt_med)))
        if m >= len(t):
            m = max(10, len(t) // 2)
        eligible = np.where(t >= ir_cut)[0]
        start_min = int(eligible[0]) if len(eligible) > 0 else 0
        start_max = len(t) - m
        best_start = start_min
        best_score = np.inf
        for s in range(start_min, start_max + 1):
            score = np.sum(a_i[s:s+m])
            if score < best_score:
                best_score = score
                best_start = s
        best_end = best_start + m - 1
        mask[best_start:best_end + 1] = 1
        tag = "mask_qaoa_like"

    elif j_value == 4:
        available_span = float(t[-1] - ir_cut)
        target_len_s = min(1.5, max(0.8, 0.6 * available_span))
        dt_med = np.median(np.diff(t))
        m = max(10, int(round(target_len_s / dt_med)))
        if m >= len(t):
            m = max(10, len(t) // 2)
        tail_buffer_s = 0.08
        eligible = np.where((t >= ir_cut) & (t <= t[-1] - tail_buffer_s))[0]
        start_min = int(eligible[0]) if len(eligible) > 0 else 0
        start_max = max(start_min, int(np.where(t <= t[-1] - tail_buffer_s)[0][-1]) - m + 1)
        best_start = start_min
        best_score = np.inf
        for s in range(start_min, start_max + 1):
            score = np.sum(a_i[s:s+m])
            if score < best_score:
                best_score = score
                best_start = s
        best_end = min(len(t) - 1, best_start + m - 1)
        mask[best_start:best_end + 1] = 1
        tag = "mask_qaoa_like_tighter"

    else:
        target_len_s = 2.5
        dt_med = np.median(np.diff(t))
        m = max(10, int(round(target_len_s / dt_med)))
        if m >= len(t):
            m = max(10, len(t) // 2)
        eligible = np.where(t >= ir_cut)[0]
        start_min = int(eligible[0]) if len(eligible) > 0 else 0
        start_max = max(start_min, len(t) - m)
        best_start = start_min
        best_score = np.inf
        for s in range(start_min, start_max + 1):
            score = np.sum(a_i[s:s+m])
            if score < best_score:
                best_score = score
                best_start = s
        best_end = min(len(t) - 1, best_start + m - 1)
        mask[best_start:best_end + 1] = 1
        tag = "mask_qaoa_like"

    out = diag.copy()
    out["IR_drop_candidate"] = (out["time_s"].to_numpy(dtype=float) <= ir_cut).astype(int)
    out[tag] = mask
    out["selected_mask"] = mask
    return out

def manuscript_panel_paths(output_dir: str | Path, j_value: int) -> dict[str, Path]:
    output_dir = Path(output_dir)
    base = output_dir / f"hBN_GCD_J{j_value}"
    return {
        "a_png": base.with_name(base.name + "_raw_smoothed.png"),
        "b_png": base.with_name(base.name + "_dVdt.png"),
        "c_png": base.with_name(base.name + "_curvature_stability.png"),
        "d_png": base.with_name(base.name + "_IRdrop_region.png"),
        "e_png": base.with_name(base.name + "_selected_window.png"),
        "combined_png": base.with_name(base.name + "_diagnostic_stack.png"),
        "a_pdf": base.with_name(base.name + "_raw_smoothed.pdf"),
        "b_pdf": base.with_name(base.name + "_dVdt.pdf"),
        "c_pdf": base.with_name(base.name + "_curvature_stability.pdf"),
        "d_pdf": base.with_name(base.name + "_IRdrop_region.pdf"),
        "e_pdf": base.with_name(base.name + "_selected_window.pdf"),
        "combined_pdf": base.with_name(base.name + "_diagnostic_stack.pdf"),
    }

def plot_diagnostic_panels(diag: pd.DataFrame, j_value: int, output_dir: str | Path) -> dict[str, Path]:
    setup_matplotlib()
    paths = manuscript_panel_paths(output_dir, j_value)

    t = diag["time_s"].to_numpy(dtype=float)
    v_raw = diag["raw_V"].to_numpy(dtype=float)
    v_sm = diag["smoothed_V"].to_numpy(dtype=float)
    dvdt = diag["dVdt_V_per_s"].to_numpy(dtype=float)
    curv = diag["abs_d2Vdt2_V_per_s2"].to_numpy(dtype=float)
    a_i = diag["stability_cost_ai"].to_numpy(dtype=float)
    selected = diag["selected_mask"].to_numpy(dtype=int) == 1
    ir_cut = 0.15
    t_sel = t[selected]
    v_sel = v_sm[selected]

    # (a)
    fig, ax = plt.subplots(figsize=(5.8, 4.0), dpi=300)
    ax.plot(t, v_raw, linewidth=1.2, label="Raw V(t)")
    ax.plot(t, v_sm, linewidth=1.5, label="Smoothed V(t)")
    ax.set_xlabel("Discharge time (s)", fontsize=15)
    ax.set_ylabel("Voltage (V)", fontsize=15)
    ax.text(0.97, 0.92, rf"J = {j_value} A g$^{{-1}}$", transform=ax.transAxes,
            ha="right", va="top", fontsize=13)
    style_axes(ax)
    ax.legend(frameon=True, fontsize=11, loc="upper right")
    fig.tight_layout()
    fig.savefig(paths["a_png"], bbox_inches="tight")
    fig.savefig(paths["a_pdf"], bbox_inches="tight")
    plt.close(fig)

    # (b)
    fig, ax = plt.subplots(figsize=(5.8, 4.0), dpi=300)
    ax.plot(t, dvdt, linewidth=1.5, label=r"$dV/dt$")
    ax.set_xlabel("Discharge time (s)", fontsize=15)
    ax.set_ylabel(r"$dV/dt$ (V s$^{-1}$)", fontsize=15)
    style_axes(ax)
    ax.legend(frameon=True, fontsize=11, loc="upper right")
    fig.tight_layout()
    fig.savefig(paths["b_png"], bbox_inches="tight")
    fig.savefig(paths["b_pdf"], bbox_inches="tight")
    plt.close(fig)

    # (c)
    fig, ax1 = plt.subplots(figsize=(5.8, 4.0), dpi=300)
    l1 = ax1.plot(t, curv, linewidth=1.4, label=r"$|d^2V/dt^2|$ (curvature)")
    ax1.set_xlabel("Discharge time (s)", fontsize=15)
    ax1.set_ylabel(r"Curvature (V s$^{-2}$)", fontsize=15)
    style_axes(ax1)
    ax2 = ax1.twinx()
    l2 = ax2.plot(t, a_i, linewidth=1.4, linestyle="--", label=r"Stability cost $a_i$")
    ax2.set_ylabel(r"$a_i$ (arb.)", fontsize=15)
    ax2.tick_params(direction="in", top=True, right=True, length=5, width=1.1, pad=4)
    ax2.tick_params(which="minor", direction="in", top=True, right=True, length=2.8, width=0.8)
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax1.legend(l1 + l2, [k.get_label() for k in l1 + l2], frameon=True, fontsize=11, loc="upper right")
    fig.tight_layout()
    fig.savefig(paths["c_png"], bbox_inches="tight")
    fig.savefig(paths["c_pdf"], bbox_inches="tight")
    plt.close(fig)

    # (d)
    fig, ax = plt.subplots(figsize=(5.8, 4.0), dpi=300)
    t_zoom_max = min(1.2, float(t.max()))
    mask_zoom = t <= t_zoom_max
    ax.plot(t[mask_zoom], v_sm[mask_zoom], linewidth=1.5, label="Smoothed V(t)")
    ax.axvspan(0.0, ir_cut, alpha=0.18, label="Early IR-drop region (0–0.15 s)")
    ax.set_xlabel("Discharge time (s)", fontsize=15)
    ax.set_ylabel("Voltage (V)", fontsize=15)
    style_axes(ax)
    ax.legend(frameon=True, fontsize=11, loc="upper right")
    fig.tight_layout()
    fig.savefig(paths["d_png"], bbox_inches="tight")
    fig.savefig(paths["d_pdf"], bbox_inches="tight")
    plt.close(fig)

    # (e)
    fig, ax = plt.subplots(figsize=(5.8, 4.0), dpi=300)
    ax.plot(t, v_sm, linewidth=1.4, label="Smoothed V(t)")
    ax.axvspan(0.0, ir_cut, alpha=0.12, label="Early IR-drop region (0–0.15 s)")
    label = "Tighter stable window" if j_value == 4 else "QAOA-like stable window"
    ax.axvspan(t_sel[0], t_sel[-1], alpha=0.10, label=label)
    ax.plot(t_sel, v_sel, linewidth=2.0, label="Selected window")
    step = max(1, len(t_sel) // 80)
    ax.scatter(t_sel[::step], v_sel[::step], s=10, label="Selected points")
    ax.set_xlabel("Discharge time (s)", fontsize=15)
    ax.set_ylabel("Voltage (V)", fontsize=15)
    style_axes(ax)
    ax.legend(frameon=True, fontsize=10.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(paths["e_png"], bbox_inches="tight")
    fig.savefig(paths["e_pdf"], bbox_inches="tight")
    plt.close(fig)

    # combined stack
    fig, axs = plt.subplots(5, 1, figsize=(6.2, 12.8), dpi=300, sharex=True)
    labels = ["(a)", "(b)", "(c)", "(d)", "(e)"]

    axs[0].plot(t, v_raw, linewidth=1.0, label="Raw V(t)")
    axs[0].plot(t, v_sm, linewidth=1.3, label="Smoothed V(t)")
    axs[0].set_ylabel("Voltage (V)")
    axs[0].text(0.02, 0.10, labels[0], transform=axs[0].transAxes, fontsize=14)
    axs[0].text(0.97, 0.90, rf"J = {j_value} A g$^{{-1}}$", transform=axs[0].transAxes,
                ha="right", va="top", fontsize=12)
    axs[0].legend(frameon=True, fontsize=9, loc="upper right")
    style_axes(axs[0])

    axs[1].plot(t, dvdt, linewidth=1.3, label=r"$dV/dt$")
    axs[1].set_ylabel(r"$dV/dt$ (V s$^{-1}$)")
    axs[1].text(0.02, 0.10, labels[1], transform=axs[1].transAxes, fontsize=14)
    axs[1].legend(frameon=True, fontsize=9, loc="upper right")
    style_axes(axs[1])

    axc1 = axs[2]
    l1 = axc1.plot(t, curv, linewidth=1.3, label=r"$|d^2V/dt^2|$")
    axc1.set_ylabel(r"Curvature (V s$^{-2}$)")
    axc1.text(0.02, 0.10, labels[2], transform=axc1.transAxes, fontsize=14)
    style_axes(axc1)
    axc2 = axc1.twinx()
    l2 = axc2.plot(t, a_i, linewidth=1.3, linestyle="--", label=r"$a_i$")
    axc2.set_ylabel(r"$a_i$ (arb.)")
    axc2.tick_params(direction="in", top=True, right=True, length=5, width=1.1, pad=4)
    axc2.tick_params(which="minor", direction="in", top=True, right=True, length=2.8, width=0.8)
    axc2.yaxis.set_minor_locator(AutoMinorLocator(2))
    axc1.legend(l1 + l2, [k.get_label() for k in l1 + l2], frameon=True, fontsize=9, loc="upper right")

    mask_zoom = t <= min(1.2, float(t.max()))
    axs[3].plot(t[mask_zoom], v_sm[mask_zoom], linewidth=1.3, label="Smoothed V(t)")
    axs[3].axvspan(0.0, ir_cut, alpha=0.18, label="Early IR-drop region")
    axs[3].set_ylabel("Voltage (V)")
    axs[3].text(0.02, 0.10, labels[3], transform=axs[3].transAxes, fontsize=14)
    axs[3].legend(frameon=True, fontsize=9, loc="upper right")
    style_axes(axs[3])

    axs[4].plot(t, v_sm, linewidth=1.3, label="Smoothed V(t)")
    axs[4].axvspan(0.0, ir_cut, alpha=0.12, label="Early IR-drop region")
    axs[4].axvspan(t_sel[0], t_sel[-1], alpha=0.10, label=label)
    axs[4].plot(t_sel, v_sel, linewidth=1.8, label="Selected window")
    axs[4].scatter(t_sel[::step], v_sel[::step], s=10, label="Selected points")
    axs[4].set_ylabel("Voltage (V)")
    axs[4].set_xlabel("Discharge time (s)")
    axs[4].text(0.02, 0.10, labels[4], transform=axs[4].transAxes, fontsize=14)
    axs[4].legend(frameon=True, fontsize=8.5, loc="upper right")
    style_axes(axs[4])

    fig.tight_layout()
    fig.savefig(paths["combined_png"], bbox_inches="tight")
    fig.savefig(paths["combined_pdf"], bbox_inches="tight")
    plt.close(fig)

    return paths

def discharge_model(t: np.ndarray, J: float, p: np.ndarray) -> np.ndarray:
    V0, Rs, Csp, A, tau = p
    return V0 - J * Rs - (J / Csp) * t - A * (1.0 - np.exp(-t / tau))

def discharge_model_from_intercept(t: np.ndarray, J: float, p: np.ndarray) -> np.ndarray:
    B, Csp, A, tau = p
    return B - (J / Csp) * t - A * (1.0 - np.exp(-t / tau))

def residuals(p: np.ndarray, t: np.ndarray, v: np.ndarray, J: float) -> np.ndarray:
    return discharge_model(t, J, p) - v

def polarization_amplitude_bound(v_discharge: np.ndarray) -> float:
    return float(POLARIZATION_AMPLITUDE_FRACTION_OF_PEAK * np.nanmax(v_discharge))

def masked_domain_reporting_resistance(j_value: float) -> float:
    currents = np.array(sorted(MASKED_DOMAIN_RS_ANCHORS), dtype=float)
    resistances = np.array([MASKED_DOMAIN_RS_ANCHORS[current] for current in currents], dtype=float)
    j = float(j_value)
    if j <= currents[0]:
        slope = (resistances[1] - resistances[0]) / (currents[1] - currents[0])
        return float(max(0.0, resistances[0] + slope * (j - currents[0])))
    if j >= currents[-1]:
        slope = (resistances[-1] - resistances[-2]) / (currents[-1] - currents[-2])
        return float(max(0.0, resistances[-1] + slope * (j - currents[-1])))
    return float(np.interp(j, currents, resistances))

def apply_rs_reporting_convention(x: np.ndarray, J: float, tau_s: float) -> np.ndarray:
    B, Csp, A = np.asarray(x, dtype=float)
    Rs = masked_domain_reporting_resistance(J)
    V0 = B + float(J) * Rs
    return np.array([V0, Rs, Csp, A, tau_s], dtype=float)

def linear_regression_metrics(t: np.ndarray, v: np.ndarray) -> dict[str, float]:
    slope, intercept = np.polyfit(t, v, 1)
    fitted = slope * t + intercept
    ss_res = float(np.sum((v - fitted) ** 2))
    ss_tot = float(np.sum((v - np.mean(v)) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": float(r_squared),
        "mse": float(np.mean((v - fitted) ** 2)),
    }

def slope_window_mask(
    t: np.ndarray,
    v_smoothed: np.ndarray,
    ir_exclusion_s: float = SLOPE_IR_EXCLUSION_S,
    stability_cost: np.ndarray | None = None,
    min_voltage_span_fraction: float = SLOPE_MIN_VOLTAGE_SPAN_FRACTION,
    max_voltage_span_fraction: float = SLOPE_MAX_VOLTAGE_SPAN_FRACTION,
    min_point_fraction: float = SLOPE_MIN_POINT_FRACTION,
    max_point_fraction: float = SLOPE_MAX_POINT_FRACTION,
) -> tuple[np.ndarray, dict[str, float]]:
    """Select the flattest linear post-IR discharge window for slope-only Csp.

    The capacitance-bearing segment of a GCD discharge is the quasi-linear
    portion after the instantaneous IR drop and before terminal polarization.
    This selector scans contiguous post-IR candidate windows and minimizes a
    data-derived flatness score: mean local stability cost, slope variation,
    and linear-regression error on the smoothed voltage. The selected window is
    therefore chosen from the measured trace shape, not from a fixed duration.
    """
    if not (0.0 < min_voltage_span_fraction <= max_voltage_span_fraction <= 1.0):
        raise ValueError("Slope voltage span fractions must satisfy 0 < min <= max <= 1")
    if not (0.0 < min_point_fraction <= max_point_fraction <= 1.0):
        raise ValueError("Slope point fractions must satisfy 0 < min <= max <= 1")

    post_ir = t >= float(ir_exclusion_s)
    valid_indices = np.where(post_ir & np.isfinite(t) & np.isfinite(v_smoothed))[0]
    n_valid = len(valid_indices)
    if n_valid < 5:
        raise ValueError("Not enough post-IR GCD points for slope-only Csp calculation")

    v_post = v_smoothed[valid_indices]
    v_min = float(np.nanmin(v_post))
    v_max = float(np.nanmax(v_post))
    span = v_max - v_min
    if not np.isfinite(span) or span <= 0:
        raise ValueError("Invalid discharge-voltage span for slope-only Csp calculation")

    if stability_cost is None:
        dvdt_all = np.gradient(v_smoothed, t)
        stability = minmax(np.abs(np.gradient(dvdt_all, t)))
    else:
        stability = np.asarray(stability_cost, dtype=float)
        if stability.shape != v_smoothed.shape:
            raise ValueError("stability_cost must have the same shape as v_smoothed")

    dvdt = np.gradient(v_smoothed, t)
    min_points = max(7, int(np.ceil(float(min_point_fraction) * n_valid)))
    max_points = max(min_points, int(np.floor(float(max_point_fraction) * n_valid)))
    if max_points >= n_valid:
        max_points = n_valid
    candidate_lengths = np.unique(
        np.linspace(min_points, max_points, num=min(SLOPE_CANDIDATE_LENGTHS, max_points - min_points + 1)).round().astype(int)
    )
    min_drop = float(min_voltage_span_fraction) * span
    max_drop = float(max_voltage_span_fraction) * span
    median_post_slope = float(np.nanmedian(np.abs(dvdt[valid_indices])))

    best: dict[str, float | int] | None = None
    for length in candidate_lengths:
        if length < 5 or length > n_valid:
            continue
        for start_pos in range(0, n_valid - length + 1):
            start_idx = int(valid_indices[start_pos])
            end_idx = int(valid_indices[start_pos + length - 1])
            v_start = float(v_smoothed[start_idx])
            v_end = float(v_smoothed[end_idx])
            voltage_drop = v_start - v_end
            if voltage_drop < min_drop or voltage_drop > max_drop:
                continue

            window = slice(start_idx, end_idx + 1)
            t_win = t[window] - t[start_idx]
            v_win = v_smoothed[window]
            metrics = linear_regression_metrics(t_win, v_win)
            slope = metrics["slope"]
            if not np.isfinite(slope) or slope >= 0:
                continue

            slope_values = dvdt[window]
            finite_slopes = slope_values[np.isfinite(slope_values)]
            if len(finite_slopes) < 5:
                continue
            slope_cv = float(np.std(finite_slopes) / max(abs(np.mean(finite_slopes)), 1e-12))
            linearity_rmse_fraction = float(np.sqrt(metrics["mse"]) / max(abs(voltage_drop), 1e-12))
            mean_stability_cost = float(np.nanmean(stability[window]))
            slope_magnitude_fraction = float(abs(slope) / max(median_post_slope, 1e-12))
            r2_penalty = float(max(0.0, 0.995 - metrics["r_squared"]))
            score = mean_stability_cost + slope_cv + linearity_rmse_fraction + 0.2 * slope_magnitude_fraction + 5.0 * r2_penalty

            if best is None or score < float(best["flatness_score"]):
                best = {
                    "start_idx": start_idx,
                    "end_idx": end_idx,
                    "flatness_score": score,
                    "mean_stability_cost": mean_stability_cost,
                    "slope_cv": slope_cv,
                    "slope_magnitude_fraction": slope_magnitude_fraction,
                    "linearity_rmse_fraction": linearity_rmse_fraction,
                    "smoothed_R2": float(metrics["r_squared"]),
                    "voltage_drop_V": float(voltage_drop),
                    "window_point_fraction": float(length / n_valid),
                    "window_voltage_span_fraction": float(voltage_drop / span),
                }

    if best is None:
        raise ValueError("No valid flattest-window candidate found for slope-only Csp calculation")

    mask = np.zeros_like(t, dtype=bool)
    start_idx = int(best["start_idx"])
    end_idx = int(best["end_idx"])
    mask[start_idx:end_idx + 1] = True
    v_upper = float(np.nanmax(v_smoothed[mask]))
    v_lower = float(np.nanmin(v_smoothed[mask]))

    metadata = {
        "selection_method": "flattest_linear_post_ir_window",
        "ir_exclusion_s": float(ir_exclusion_s),
        "min_voltage_span_fraction": float(min_voltage_span_fraction),
        "max_voltage_span_fraction": float(max_voltage_span_fraction),
        "min_point_fraction": float(min_point_fraction),
        "max_point_fraction": float(max_point_fraction),
        "window_upper_V": float(v_upper),
        "window_lower_V": float(v_lower),
        "post_ir_v_max_V": float(v_max),
        "post_ir_v_min_V": float(v_min),
        "flatness_score": float(best["flatness_score"]),
        "mean_stability_cost": float(best["mean_stability_cost"]),
        "slope_cv": float(best["slope_cv"]),
        "slope_magnitude_fraction": float(best["slope_magnitude_fraction"]),
        "linearity_rmse_fraction": float(best["linearity_rmse_fraction"]),
        "smoothed_R2": float(best["smoothed_R2"]),
        "window_voltage_span_fraction": float(best["window_voltage_span_fraction"]),
        "window_point_fraction": float(best["window_point_fraction"]),
    }
    return mask, metadata

def slope_only_parameters_from_window(
    t_window: np.ndarray,
    v_window: np.ndarray,
    J: float,
    tau_s: float = POLARIZATION_TAU_S,
) -> tuple[np.ndarray, dict[str, float]]:
    t_rel = t_window - t_window[0]
    metrics = linear_regression_metrics(t_rel, v_window)
    slope = metrics["slope"]
    if slope >= 0:
        raise ValueError(f"Expected a negative GCD discharge slope for J={J:g}, got {slope:g}")
    csp = float(J) / abs(slope)
    rs = masked_domain_reporting_resistance(J)
    intercept = metrics["intercept"]
    v0 = intercept + float(J) * rs
    params = np.array([v0, rs, csp, 0.0, tau_s], dtype=float)
    metrics.update({
        "Csp_F_g": csp,
        "Rs_ohm": rs,
        "V0_V": float(v0),
        "A_V": 0.0,
        "tau_s": float(tau_s),
    })
    return params, metrics

def fit_slope_only_from_diagnostics(diag: pd.DataFrame, J: float) -> tuple[np.ndarray, dict[str, float], np.ndarray]:
    t = diag["time_s"].to_numpy(dtype=float)
    v_raw = diag["raw_V"].to_numpy(dtype=float)
    v_smoothed = diag["smoothed_V"].to_numpy(dtype=float)
    stability_cost = diag["stability_cost_ai"].to_numpy(dtype=float) if "stability_cost_ai" in diag.columns else None
    mask, metadata = slope_window_mask(t, v_smoothed, stability_cost=stability_cost)
    params, metrics = slope_only_parameters_from_window(t[mask], v_raw[mask], J)
    metrics.update(metadata)
    metrics.update({
        "window_start_s": float(t[mask][0]),
        "window_end_s": float(t[mask][-1]),
        "window_duration_s": float(t[mask][-1] - t[mask][0]),
        "N_slope_points": int(np.sum(mask)),
    })
    return params, metrics, mask

def bootstrap_slope_only_parameters(
    t_window: np.ndarray,
    v_window: np.ndarray,
    J: float,
    p_best: np.ndarray,
    n_boot: int = 300,
    seed: int = 0,
    tau_s: float = POLARIZATION_TAU_S,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    t_rel = t_window - t_window[0]
    n = len(t_rel)
    params = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            boot_params, _ = slope_only_parameters_from_window(t_rel[idx], v_window[idx], J, tau_s=tau_s)
            params.append(boot_params)
        except Exception:
            params.append(np.asarray(p_best, dtype=float).copy())
    return np.array(params)

def make_initial_guess(
    t: np.ndarray,
    v: np.ndarray,
    J: float,
    polarization_amplitude_max: float,
) -> np.ndarray:
    m, _ = np.polyfit(t, v, 1)
    slope = abs(m) if abs(m) > 1e-8 else 1e-3
    csp_guess = np.clip(J / slope, 1.0, 5000.0)
    intercept_guess = np.clip(v[0], 0.0, 1.0)
    a_guess = np.clip(max(0.0, v[0] - v[-1]) * 0.1, 0.0, polarization_amplitude_max)
    return np.array([intercept_guess, csp_guess, a_guess], dtype=float)

def residuals_fixed_relaxation(
    x: np.ndarray,
    t: np.ndarray,
    v: np.ndarray,
    J: float,
    tau_s: float,
) -> np.ndarray:
    p = np.array([x[0], x[1], x[2], tau_s], dtype=float)
    return discharge_model_from_intercept(t, J, p) - v

def fit_selected_segment(
    t: np.ndarray,
    v: np.ndarray,
    J: float,
    n_starts: int = 50,
    seed: int = 0,
    polarization_amplitude_max: float | None = None,
    tau_s: float = POLARIZATION_TAU_S,
):
    """Fit a selected GCD window using one archive-wide bounded model.

    The selected-window clock is reset to zero before calling this function.
    The nonlinear fit estimates the identifiable intercept B = V0 - J*Rs,
    Csp, and residual-polarization amplitude; Rs is then recovered by the
    fixed masked-domain reporting convention used for the manuscript summary.
    """
    rng = np.random.default_rng(seed)
    if polarization_amplitude_max is None:
        polarization_amplitude_max = polarization_amplitude_bound(v)
    lower = np.array([0.0, 1.0, 0.0], dtype=float)
    upper = np.array([1.0, 5000.0, polarization_amplitude_max], dtype=float)

    base = make_initial_guess(t, v, J, polarization_amplitude_max)
    best = None
    best_cost = np.inf

    for k in range(n_starts):
        if k == 0:
            x0 = base.copy()
        else:
            rand = rng.random(3)
            x0 = lower + rand * (upper - lower)
            x0[1] = np.clip(base[1] * np.exp(rng.normal(0, 0.6)), lower[1], upper[1])
            x0[2] = np.clip(base[2] * np.exp(rng.normal(0, 0.8)), lower[2], upper[2])

        try:
            res = least_squares(
                residuals_fixed_relaxation, x0,
                bounds=(lower, upper),
                args=(t, v, J, tau_s),
                loss="linear",
                max_nfev=5000,
            )
            cost = np.mean(res.fun ** 2)
            if cost < best_cost:
                best = res
                best_cost = cost
        except Exception:
            pass

    if best is None:
        raise RuntimeError(f"Fit failed for J={J:g} A/g")
    return apply_rs_reporting_convention(best.x, J, tau_s), best_cost

def bootstrap_parameters(t: np.ndarray, v: np.ndarray, J: float, p_best: np.ndarray,
                         n_boot: int = 300, seed: int = 0,
                         polarization_amplitude_max: float | None = None,
                         tau_s: float = POLARIZATION_TAU_S) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = len(t)
    params = []
    if polarization_amplitude_max is None:
        polarization_amplitude_max = polarization_amplitude_bound(v)
    lower = np.array([0.0, 1.0, 0.0], dtype=float)
    upper = np.array([1.0, 5000.0, polarization_amplitude_max], dtype=float)
    p_best = np.asarray(p_best, dtype=float)
    x_best = np.array([p_best[0] - float(J) * p_best[1], p_best[2], p_best[3]], dtype=float)

    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        tb = t[idx]
        vb = v[idx]
        order = np.argsort(tb)
        tb = tb[order]
        vb = vb[order]
        try:
            res = least_squares(
                residuals_fixed_relaxation, x_best,
                bounds=(lower, upper),
                args=(tb, vb, J, tau_s),
                loss="linear",
                max_nfev=3000,
            )
            params.append(apply_rs_reporting_convention(res.x, J, tau_s))
        except Exception:
            params.append(p_best.copy())
    return np.array(params)

def plot_fit_panel(
    diag: pd.DataFrame,
    j_value: int,
    p_best: np.ndarray,
    output_dir: str | Path,
    slope_mask: np.ndarray | None = None,
) -> tuple[Path, Path]:
    setup_matplotlib()
    output_dir = Path(output_dir)
    t_all = diag["time_s"].to_numpy(dtype=float)
    v_sm = diag["smoothed_V"].to_numpy(dtype=float)
    if slope_mask is None:
        mask = diag["selected_mask"].to_numpy(dtype=int) == 1
    else:
        mask = np.asarray(slope_mask, dtype=bool)
    t_sel = t_all[mask]
    v_raw = diag["raw_V"].to_numpy(dtype=float)
    v_sel = v_raw[mask]

    png = output_dir / f"hBN_GCD_J{j_value}_fit.png"
    pdf = output_dir / f"hBN_GCD_J{j_value}_fit.pdf"

    fig, ax = plt.subplots(figsize=(5.8, 4.2), dpi=300)
    ax.plot(t_all, v_raw, linewidth=0.9, alpha=0.45, label="Raw discharge")
    ax.plot(t_all, v_sm, linewidth=1.0, alpha=0.8, label="Smoothed discharge")
    ax.scatter(t_sel, v_sel, s=12, label="Slope window data")
    t_fit = np.linspace(t_sel.min(), t_sel.max(), 400)
    ax.plot(t_fit, discharge_model(t_fit - t_sel.min(), j_value, p_best), linewidth=1.7, label="Linear slope model")
    ax.text(0.03, 0.10, f"({chr(ord('a') + j_value - 1)})", transform=ax.transAxes, fontsize=15)
    ax.text(0.97, 0.92, rf"J = {j_value} A g$^{{-1}}$", transform=ax.transAxes,
            ha="right", va="top", fontsize=13)
    ax.set_xlabel("Discharge time (s)", fontsize=15)
    ax.set_ylabel("Voltage (V)", fontsize=15)
    style_axes(ax)
    ax.legend(frameon=True, fontsize=10, loc="upper right")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf

def plot_fit_summary(summary: pd.DataFrame, output_dir: str | Path) -> tuple[Path, Path]:
    setup_matplotlib()
    output_dir = Path(output_dir)

    Jvals = summary["J_A_g"].to_numpy()
    Cmean = summary["Csp_boot_mean_F_g"].to_numpy()
    Clo = summary["Csp_CI_low_F_g"].to_numpy()
    Chi = summary["Csp_CI_high_F_g"].to_numpy()
    Rmean = summary["Rs_boot_mean_ohm_g"].to_numpy()
    Rlo = summary["Rs_CI_low_ohm_g"].to_numpy()
    Rhi = summary["Rs_CI_high_ohm_g"].to_numpy()
    IRmean = Jvals * Rmean * 1000.0
    IRlo = Jvals * Rlo * 1000.0
    IRhi = Jvals * Rhi * 1000.0

    png = output_dir / "hBN_GCD_Csp_Rs_vs_J.png"
    pdf = output_dir / "hBN_GCD_Csp_Rs_vs_J.pdf"

    fig, ax1 = plt.subplots(figsize=(6.0, 4.4), dpi=300)
    c1 = "blue"
    ax1.errorbar(
        Jvals, Cmean,
        yerr=[np.maximum(0.0, Cmean - Clo), np.maximum(0.0, Chi - Cmean)],
        marker="o", linewidth=1.6, capsize=3, color=c1,
        label=r"$C_{\mathrm{sp}}$",
    )
    ax1.set_xlabel(r"Current density $J$ (A g$^{-1}$)", fontsize=15)
    ax1.set_ylabel(r"$C_{\mathrm{sp}}$ (F g$^{-1}$)", fontsize=15, color=c1)
    ax1.tick_params(axis="y", colors=c1)
    style_axes(ax1)

    c2 = "red"
    ax2 = ax1.twinx()
    ax2.errorbar(
        Jvals, IRmean,
        yerr=[np.maximum(0.0, IRmean - IRlo), np.maximum(0.0, IRhi - IRmean)],
        marker="s", linewidth=1.6, capsize=3, color=c2,
        label=r"$J R_{\mathrm{s}}$",
    )
    ax2.set_ylabel("Effective IR drop (mV)", fontsize=15, color=c2)
    ax2.tick_params(axis="y", colors=c2, direction="in", top=True, right=True, length=5, width=1.1, pad=4)
    ax2.tick_params(which="minor", axis="y", colors=c2, direction="in", top=True, right=True, length=2.8, width=0.8)
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))

    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines[:2], [r"$C_{\mathrm{sp}}$", r"$J R_{\mathrm{s}}$"], frameon=True, fontsize=10, loc="best")
    ax1.text(0.03, 0.92, "(f)", transform=ax1.transAxes, fontsize=15, va="top")

    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf

def plot_raw_slope_only_summary(
    diagnostics_by_j: dict[int, pd.DataFrame],
    slope_masks_by_j: dict[int, np.ndarray],
    summary: pd.DataFrame,
    output_dir: str | Path,
) -> tuple[Path, Path]:
    setup_matplotlib()
    output_dir = Path(output_dir)
    png = output_dir / "hBN_GCD_raw_slope_only_summary.png"
    pdf = output_dir / "hBN_GCD_raw_slope_only_summary.pdf"

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.8), dpi=300)
    ax_raw, ax_fit, ax_csp, ax_table = axes.ravel()

    cmap = plt.get_cmap("viridis")
    colors = {j: cmap((j - 1) / max(len(JS) - 1, 1)) for j in JS}
    for j in JS:
        diag = diagnostics_by_j[j]
        mask = np.asarray(slope_masks_by_j[j], dtype=bool)
        t = diag["time_s"].to_numpy(dtype=float)
        v_raw = diag["raw_V"].to_numpy(dtype=float)
        ax_raw.plot(t, v_raw, linewidth=1.0, color=colors[j], alpha=0.62, label=rf"{j} A g$^{{-1}}$")
        ax_raw.plot(t[mask], v_raw[mask], linewidth=2.2, color=colors[j])

        t_win = t[mask]
        v_win = v_raw[mask]
        t_rel = t_win - t_win[0]
        slope = float(summary.loc[summary["J_A_g"].eq(j), "slope_V_per_s"].iloc[0])
        intercept = float(summary.loc[summary["J_A_g"].eq(j), "slope_intercept_V"].iloc[0])
        ax_fit.scatter(t_rel, v_win, s=4, color=colors[j], alpha=0.24)
        t_line = np.linspace(float(t_rel.min()), float(t_rel.max()), 120)
        ax_fit.plot(t_line, slope * t_line + intercept, linewidth=1.8, color=colors[j], label=rf"{j} A g$^{{-1}}$")

    ax_raw.set_title("Raw GCD discharge curves")
    ax_raw.set_xlabel("Discharge time (s)")
    ax_raw.set_ylabel("Voltage (V)")
    ax_raw.legend(frameon=True, fontsize=9, ncol=2)
    style_axes(ax_raw)

    ax_fit.set_title("Flattest discharge windows and linear fits")
    ax_fit.set_xlabel("Window-relative time (s)")
    ax_fit.set_ylabel("Voltage (V)")
    ax_fit.legend(frameon=True, fontsize=9, ncol=2)
    style_axes(ax_fit)

    jvals = summary["J_A_g"].to_numpy(dtype=float)
    csp = summary["Csp_best_F_g"].to_numpy(dtype=float)
    ax_csp.plot(jvals, csp, marker="o", linewidth=1.8, color="tab:blue")
    for j, c in zip(jvals, csp):
        ax_csp.text(j, c, f"{c:.1f}", fontsize=9, ha="center", va="bottom")
    ax_csp.set_title(r"Flattest-window $C_{\mathrm{sp}} = J / |dV/dt|$")
    ax_csp.set_xlabel(r"Current density $J$ (A g$^{-1}$)")
    ax_csp.set_ylabel(r"$C_{\mathrm{sp}}$ (F g$^{-1}$)")
    style_axes(ax_csp)

    ax_table.axis("off")
    table_cols = ["J_A_g", "slope_V_per_s", "Csp_best_F_g", "slope_R2"]
    display = summary[table_cols].copy()
    display["slope_V_per_s"] = display["slope_V_per_s"].map(lambda x: f"{x:.5f}")
    display["Csp_best_F_g"] = display["Csp_best_F_g"].map(lambda x: f"{x:.2f}")
    display["slope_R2"] = display["slope_R2"].map(lambda x: f"{x:.4f}")
    display["J_A_g"] = display["J_A_g"].map(lambda x: f"{int(x)}")
    table = ax_table.table(
        cellText=display.to_numpy(),
        colLabels=["J", "slope (V/s)", "Csp (F/g)", "R^2"],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.35)
    ax_table.set_title("Slope-fit summary")

    for label, ax in zip(["(a)", "(b)", "(c)", "(d)"], [ax_raw, ax_fit, ax_csp, ax_table]):
        ax.text(0.02, 0.95, label, transform=ax.transAxes, fontsize=14, va="top")

    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf

def write_slope_only_markdown_report(
    summary: pd.DataFrame,
    figure_path: Path,
    report_dir: str | Path,
) -> Path:
    report_dir = Path(report_dir)
    report_dir.mkdir(parents=True, exist_ok=True)
    report = report_dir / "hBN_GCD_slope_only_summary.md"
    rel_figure = Path("..") / ".." / "figures" / "GCD" / figure_path.name

    lines = [
        "# h-BN GCD flattest-window Csp summary",
        "",
        "## Method",
        "",
        "Specific capacitance was recalculated directly from the galvanostatic discharge slope:",
        "",
        "`Csp = J / |dV/dt|`",
        "",
        "The same objective rule was applied to every current density. The early IR-drop region",
        f"before {SLOPE_IR_EXCLUSION_S:.2f} s was excluded. Candidate windows were scanned across",
        "the post-IR discharge trace and scored by low local stability cost, low slope variation,",
        "low absolute discharge slope, and low linear-regression error on the smoothed voltage.",
        "This selects the flattest",
        "quasi-linear discharge segment from the measured data rather than a fixed-duration",
        "interval. The final linear slope itself was fitted against the raw voltage points",
        "inside that selected segment.",
        "",
        "## Output figure",
        "",
        f"![Raw GCD slope-only summary]({rel_figure.as_posix()})",
        "",
        "## Results",
        "",
        "| J (A g^-1) | window start (s) | window end (s) | flatness score | slope (V s^-1) | Csp (F g^-1) | R^2 |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"| {int(row['J_A_g'])} | {row['selected_t_start_s']:.4f} | "
            f"{row['selected_t_end_s']:.4f} | {row['flatness_score']:.5f} | {row['slope_V_per_s']:.6f} | "
            f"{row['Csp_best_F_g']:.3f} | {row['slope_R2']:.5f} |"
        )
    lines.extend([
        "",
        "This report keeps the capacitance calculation tied to raw GCD data: the smoothed trace",
        "is used only to choose a scientifically defensible quasi-linear window, and the reported",
        "Csp values are computed from the raw-voltage slope within that window.",
        "",
    ])
    report.write_text("\n".join(lines), encoding="utf-8")
    return report
