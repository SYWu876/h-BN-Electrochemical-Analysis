#!/usr/bin/env python3
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import curve_fit
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

def baseline_als(y, lam=1e7, p=0.01, niter=20):
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2), dtype=float)
    w = np.ones(L, dtype=float)
    for _ in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.T)
        z = spsolve(Z.tocsc(), w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z

def pseudo_voigt(x, A, x0, sigma, gamma, eta):
    gauss = np.exp(-((x - x0) ** 2) / (2 * sigma ** 2))
    lorentz = gamma ** 2 / ((x - x0) ** 2 + gamma ** 2)
    return A * (eta * lorentz + (1 - eta) * gauss)

def save_fig_ns14a(x, y, baseline, x0, fwhm, y_model_full, out):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":18,"axes.labelsize":22,"xtick.labelsize":18,"ytick.labelsize":18,"legend.fontsize":15})
    fig, ax = plt.subplots(figsize=(8.2, 6.4), dpi=300)
    ax.plot(x, y, linewidth=2.6, label="Experimental")
    ax.plot(x, baseline, linewidth=2.2, label="ALS baseline")
    ax.plot(x, baseline + y_model_full, linewidth=2.6, label="Pseudo-Voigt fit")
    ax.axvline(x0, linestyle="--", linewidth=1.8)
    ax.set_xlabel("Raman shift (cm$^{-1}$)")
    ax.set_ylabel("Intensity (a.u.)")
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)
    ax.tick_params(direction="in", length=7, width=1.6, top=True, right=True, pad=8)
    ax.text(0.98, 0.96, f"center = {x0:.2f} cm$^{{-1}}$\nFWHM = {fwhm:.2f} cm$^{{-1}}$", transform=ax.transAxes, ha="right", va="top", fontsize=16)
    leg = ax.legend(loc="upper left", bbox_to_anchor=(0.08, 1.0), frameon=True)
    leg.get_frame().set_linewidth(1.2)
    ax.margins(x=0.02)
    fig.tight_layout()
    for ext in ("png","pdf","svg"):
        fig.savefig(out/f"Raman_raw_baseline_fit.{ext}", dpi=600 if ext=="png" else None, bbox_inches="tight")
    plt.close(fig)

def save_fig_ns14b(x, y_corr, x_fit, y_model_fit, x0, fwhm, x_left, x_right, out):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":18,"axes.labelsize":22,"xtick.labelsize":18,"ytick.labelsize":18,"legend.fontsize":15})
    fig, ax = plt.subplots(figsize=(8.2, 6.4), dpi=300)
    plot_mask = (x >= 1335) & (x <= 1395)
    ax.plot(x[plot_mask], y_corr[plot_mask], linewidth=2.6, label="Baseline-corrected")
    ax.plot(x_fit, y_model_fit, linewidth=2.6, label="Pseudo-Voigt fit")
    ax.axvline(x0, linestyle="--", linewidth=1.8, label="Peak center")
    half_max = y_model_fit.max() / 2.0
    ax.hlines(half_max, x_left, x_right, linewidth=2.0, label="FWHM")
    ax.vlines([x_left, x_right], half_max * 0.96, half_max * 1.04, linewidth=1.6)
    ax.set_xlabel("Raman shift (cm$^{-1}$)")
    ax.set_ylabel("Baseline-corrected intensity (a.u.)")
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)
    ax.tick_params(direction="in", length=7, width=1.6, top=True, right=True, pad=8)
    ax.text(0.98, 0.96, f"center = {x0:.2f} cm$^{{-1}}$\nFWHM = {fwhm:.2f} cm$^{{-1}}$", transform=ax.transAxes, ha="right", va="top", fontsize=16)
    leg = ax.legend(loc="upper left", frameon=True)
    leg.get_frame().set_linewidth(1.2)
    ax.margins(x=0.02, y=0.08)
    fig.tight_layout()
    for ext in ("png","pdf","svg"):
        fig.savefig(out/f"Raman_baseline_corrected_fit.{ext}", dpi=600 if ext=="png" else None, bbox_inches="tight")
    plt.close(fig)

def save_fig_ns14c(xp, yp, labels, name_map, low_lab, mid_lab, high_lab, out):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":18,"axes.labelsize":22,"xtick.labelsize":18,"ytick.labelsize":18,"legend.fontsize":15})
    fig, ax = plt.subplots(figsize=(8.2, 6.4), dpi=300)
    ax.plot(xp, yp, linewidth=1.8, alpha=0.9, label="Baseline-corrected profile")
    for lab in [low_lab, mid_lab, high_lab]:
        mask = labels == lab
        ax.scatter(xp[mask], yp[mask], s=24, label=name_map[lab])
    ax.text(1371.5, np.interp(1368.5, xp, yp) + 55, "Peak core", ha="left", va="bottom", fontsize=15)
    low_x = xp[labels == low_lab]
    if len(low_x) > 0:
        xl = float(np.min(low_x)) + 1.0
        ax.text(xl, np.interp(xl, xp, yp) + 45, "Flank I", ha="left", va="bottom", fontsize=14)
    mid_x = xp[labels == mid_lab]
    if len(mid_x) > 0:
        xr = float(np.max(mid_x)) - 1.0
        ax.text(xr, np.interp(xr, xp, yp) + 45, "Flank II", ha="right", va="bottom", fontsize=14)
    ax.set_xlabel("Raman shift (cm$^{-1}$)")
    ax.set_ylabel("Baseline-corrected intensity (a.u.)")
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)
    ax.tick_params(direction="in", length=7, width=1.6, top=True, right=True, pad=8)
    leg = ax.legend(loc="upper left", frameon=True)
    leg.get_frame().set_linewidth(1.2)
    ax.margins(x=0.02, y=0.10)
    fig.tight_layout()
    for ext in ("png","pdf","svg"):
        fig.savefig(out/f"Raman_segmentation.{ext}", dpi=600 if ext=="png" else None, bbox_inches="tight")
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/raw/Raman/hBN-3.txt")
    parser.add_argument("--processed-dir", default="data/processed/Raman")
    parser.add_argument("--figures-dir", default="outputs/figures/Raman")
    args = parser.parse_args()

    input_path = Path(args.input)
    processed_dir = Path(args.processed_dir)
    figures_dir = Path(args.figures_dir)
    processed_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, sep=r"\s+|\t+|,+", engine="python", header=None, names=["RamanShift","Intensity"])
    x = df["RamanShift"].to_numpy(float)
    y = df["Intensity"].to_numpy(float)

    baseline = baseline_als(y)
    y_corr = y - baseline
    fit_mask = (x > 1330) & (x < 1400)
    x_fit = x[fit_mask]
    y_fit = y_corr[fit_mask]
    p0 = [y_fit.max(), x_fit[np.argmax(y_fit)], 5.0, 8.0, 0.5]
    bounds = ([0, 1355, 0.5, 0.5, 0], [np.inf, 1375, 30, 30, 1])
    params, _ = curve_fit(pseudo_voigt, x_fit, y_fit, p0=p0, bounds=bounds, maxfev=50000)
    A, x0, sigma, gamma, eta = params
    y_model_fit = pseudo_voigt(x_fit, *params)
    y_model_full = pseudo_voigt(x, *params)

    half_max = y_model_fit.max()/2
    idx = np.where(y_model_fit >= half_max)[0]
    x_left, x_right = x_fit[idx[0]], x_fit[idx[-1]]
    fwhm = x_right - x_left
    peak_max_raw = x[np.argmax(y)]
    left_hw = x0 - x_left
    right_hw = x_right - x0
    asym_ratio = right_hw / left_hw if left_hw != 0 else np.nan

    seg_mask = (x >= 1335) & (x <= 1395)
    xp = x[seg_mask]
    yp = y_corr[seg_mask]
    dy = np.gradient(yp, xp)
    d2y = np.gradient(dy, xp)
    X = StandardScaler().fit_transform(np.column_stack([yp, dy, d2y]))
    labels = KMeans(n_clusters=3, random_state=42, n_init=20).fit_predict(X)
    cluster_means = [(lab, yp[labels == lab].mean()) for lab in range(3)]
    cluster_means_sorted = sorted(cluster_means, key=lambda t: t[1])
    low_lab = cluster_means_sorted[0][0]
    mid_lab = cluster_means_sorted[1][0]
    high_lab = cluster_means_sorted[2][0]
    name_map = {low_lab:"Flank region I", mid_lab:"Flank region II", high_lab:"Peak core"}

    pd.DataFrame({"RamanShift_cm-1":x, "Intensity_raw":y, "Baseline_ALS":baseline, "Intensity_corrected":y_corr}).to_csv(processed_dir/"raman_baseline_corrected.csv", index=False)
    pd.DataFrame({"RamanShift_cm-1":x_fit, "Intensity_corrected":y_fit, "PseudoVoigt_fit":y_model_fit}).to_csv(processed_dir/"raman_fit_window.csv", index=False)
    pd.DataFrame([{
        "raw_peak_max_cm-1":float(peak_max_raw), "fit_center_cm-1":float(x0), "fwhm_cm-1":float(fwhm),
        "left_half_width_cm-1":float(left_hw), "right_half_width_cm-1":float(right_hw),
        "asymmetry_ratio_right_over_left":float(asym_ratio), "A":float(A), "sigma":float(sigma), "gamma":float(gamma), "eta":float(eta),
        "als_lambda":1e7, "als_p":0.01, "als_niter":20, "kmeans_k":3, "segmentation_window_min_cm-1":1335.0, "segmentation_window_max_cm-1":1395.0
    }]).to_csv(processed_dir/"raman_fit_summary.csv", index=False)
    pd.DataFrame({
        "RamanShift_cm-1":xp, "Intensity_corrected":yp, "dI_domega":dy, "d2I_domega2":d2y,
        "cluster_id":labels, "region_label":[name_map[i] for i in labels]
    }).to_csv(processed_dir/"raman_local_descriptors_and_segmentation.csv", index=False)
    rows=[]
    for lab in [low_lab, mid_lab, high_lab]:
        xr = xp[labels==lab]
        yr = yp[labels==lab]
        rows.append({"region_label":name_map[lab], "cluster_id":int(lab), "x_min_cm-1":float(xr.min()), "x_max_cm-1":float(xr.max()), "n_points":int(len(xr)), "mean_corrected_intensity":float(np.mean(yr))})
    pd.DataFrame(rows).to_csv(processed_dir/"raman_segment_ranges.csv", index=False)

    save_fig_ns14a(x, y, baseline, x0, fwhm, y_model_full, figures_dir)
    save_fig_ns14b(x, y_corr, x_fit, y_model_fit, x0, fwhm, x_left, x_right, figures_dir)
    save_fig_ns14c(xp, yp, labels, name_map, low_lab, mid_lab, high_lab, figures_dir)
    print("Done.")

if __name__ == "__main__":
    main()
