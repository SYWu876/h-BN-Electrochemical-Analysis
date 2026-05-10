
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.signal import savgol_filter
from scipy.optimize import least_squares

JS = [1, 2, 3, 4, 5]

MANUSCRIPT_TAGS = {
    1: "Figure5",
    2: "FigureS1",
    3: "FigureS2",
    4: "FigureS3",
    5: "FigureS4",
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
    tag = MANUSCRIPT_TAGS[j_value]
    base = output_dir / f"hBN_GCD_{tag}"
    return {
        "a_png": base.with_name(base.name + f"a_J{j_value}_raw_smoothed.png"),
        "b_png": base.with_name(base.name + f"b_J{j_value}_dVdt.png"),
        "c_png": base.with_name(base.name + f"c_J{j_value}_curvature_stability.png"),
        "d_png": base.with_name(base.name + f"d_J{j_value}_IRdrop_region.png"),
        "e_png": base.with_name(base.name + f"e_J{j_value}_selected_window.png"),
        "combined_png": base.with_name(base.name + f"a_e_J{j_value}_combined.png"),
        "a_pdf": base.with_name(base.name + f"a_J{j_value}_raw_smoothed.pdf"),
        "b_pdf": base.with_name(base.name + f"b_J{j_value}_dVdt.pdf"),
        "c_pdf": base.with_name(base.name + f"c_J{j_value}_curvature_stability.pdf"),
        "d_pdf": base.with_name(base.name + f"d_J{j_value}_IRdrop_region.pdf"),
        "e_pdf": base.with_name(base.name + f"e_J{j_value}_selected_window.pdf"),
        "combined_pdf": base.with_name(base.name + f"a_e_J{j_value}_combined.pdf"),
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

def residuals(p: np.ndarray, t: np.ndarray, v: np.ndarray, J: float) -> np.ndarray:
    return discharge_model(t, J, p) - v

def make_initial_guess(t: np.ndarray, v: np.ndarray, J: float) -> np.ndarray:
    m, _ = np.polyfit(t, v, 1)
    slope = abs(m) if abs(m) > 1e-8 else 1e-3
    csp_guess = np.clip(J / slope, 1.0, 5000.0)
    v0_guess = np.clip(v[0] + 0.01, 0.05, 1.0)
    rs_guess = np.clip(max(0.0, (np.max(v) - v[0]) / max(J, 1e-6) * 0.1), 0.0, 0.02)
    a_guess = np.clip(max(0.0, v[0] - v[-1]) * 0.1, 0.0, 0.3)
    tau_guess = np.clip(0.25 * (t[-1] - t[0]), 1e-4, 50.0)
    return np.array([v0_guess, rs_guess, csp_guess, a_guess, tau_guess], dtype=float)

def fit_selected_segment(t: np.ndarray, v: np.ndarray, J: float, n_starts: int = 50, seed: int = 0):
    rng = np.random.default_rng(seed)
    lower = np.array([0.05, 0.0, 1.0, 0.0, 1e-4], dtype=float)
    upper = np.array([1.00, 0.02, 5000.0, 0.30, 50.0], dtype=float)

    base = make_initial_guess(t, v, J)
    best = None
    best_cost = np.inf

    for k in range(n_starts):
        if k == 0:
            x0 = base.copy()
        else:
            rand = rng.random(5)
            x0 = lower + rand * (upper - lower)
            x0[2] = np.clip(base[2] * np.exp(rng.normal(0, 0.6)), lower[2], upper[2])
            x0[4] = np.clip(base[4] * np.exp(rng.normal(0, 0.8)), lower[4], upper[4])

        try:
            res = least_squares(
                residuals, x0,
                bounds=(lower, upper),
                args=(t, v, J),
                loss="soft_l1",
                f_scale=1e-3,
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
    return best.x, best_cost

def bootstrap_parameters(t: np.ndarray, v: np.ndarray, J: float, p_best: np.ndarray,
                         n_boot: int = 300, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = len(t)
    params = []
    lower = np.array([0.05, 0.0, 1.0, 0.0, 1e-4], dtype=float)
    upper = np.array([1.00, 0.02, 5000.0, 0.30, 50.0], dtype=float)

    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        tb = t[idx]
        vb = v[idx]
        order = np.argsort(tb)
        tb = tb[order]
        vb = vb[order]
        try:
            res = least_squares(
                residuals, p_best,
                bounds=(lower, upper),
                args=(tb, vb, J),
                loss="soft_l1",
                f_scale=1e-3,
                max_nfev=3000,
            )
            params.append(res.x)
        except Exception:
            params.append(p_best.copy())
    return np.array(params)

def plot_fit_panel(diag: pd.DataFrame, j_value: int, p_best: np.ndarray,
                   output_dir: str | Path) -> tuple[Path, Path]:
    setup_matplotlib()
    output_dir = Path(output_dir)
    t_all = diag["time_s"].to_numpy(dtype=float)
    v_sm = diag["smoothed_V"].to_numpy(dtype=float)
    mask = diag["selected_mask"].to_numpy(dtype=int) == 1
    t_sel = t_all[mask]
    v_sel = diag["raw_V"].to_numpy(dtype=float)[mask]

    panel_idx = {1: "6a", 2: "6b", 3: "6c", 4: "6d", 5: "6e"}[j_value]
    png = output_dir / f"hBN_Figure{panel_idx}_J{j_value}_fit.png"
    pdf = output_dir / f"hBN_Figure{panel_idx}_J{j_value}_fit.pdf"

    fig, ax = plt.subplots(figsize=(5.8, 4.2), dpi=300)
    ax.plot(t_all, v_sm, linewidth=1.0, alpha=0.8, label="Smoothed discharge")
    ax.scatter(t_sel, v_sel, s=12, label="Selected data")
    t_fit = np.linspace(t_sel.min(), t_sel.max(), 400)
    ax.plot(t_fit, discharge_model(t_fit, j_value, p_best), linewidth=1.7, label="Best-fit model")
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

    png = output_dir / "hBN_Figure6f_Csp_Rs_vs_J_replot.png"
    pdf = output_dir / "hBN_Figure6f_Csp_Rs_vs_J_replot.pdf"

    fig, ax1 = plt.subplots(figsize=(6.0, 4.4), dpi=300)
    c1 = "blue"
    ax1.errorbar(
        Jvals, Cmean,
        yerr=[Cmean - Clo, Chi - Cmean],
        marker="o", linewidth=1.6, capsize=3, color=c1,
        label=r"$C_{\mathrm{sp}}$",
    )
    ax1.set_xlabel(r"Current density $J$ (A g$^{-1}$)", fontsize=15)
    ax1.set_ylabel(r"$C_{\mathrm{sp}}$ (F g$^{-1}$)", fontsize=15, color=c1)
    ax1.tick_params(axis="y", colors=c1)
    ax1.set_ylim(60, 325)
    style_axes(ax1)

    c2 = "red"
    ax2 = ax1.twinx()
    ax2.errorbar(
        Jvals, Rmean,
        yerr=[Rmean - Rlo, Rhi - Rmean],
        marker="s", linewidth=1.6, capsize=3, color=c2,
        label=r"$R_{\mathrm{s}}$",
    )
    ax2.set_ylabel(r"$R_{\mathrm{s}}$ ($\Omega$ g)", fontsize=15, color=c2)
    ax2.tick_params(axis="y", colors=c2, direction="in", top=True, right=True, length=5, width=1.1, pad=4)
    ax2.tick_params(which="minor", axis="y", colors=c2, direction="in", top=True, right=True, length=2.8, width=0.8)
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.set_ylim(0.003, 0.014)

    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines[:2], [r"$C_{\mathrm{sp}}$", r"$R_{\mathrm{s}}$"], frameon=True, fontsize=10, loc="best")
    ax1.text(0.03, 0.10, "(f)", transform=ax1.transAxes, fontsize=15)

    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf
