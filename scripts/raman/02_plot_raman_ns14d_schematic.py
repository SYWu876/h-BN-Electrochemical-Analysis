#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path

def add_box(ax, x, y, w, h, title, lines, title_fs=16, text_fs=11):
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.02", linewidth=1.8, facecolor="white", edgecolor="black")
    ax.add_patch(box)
    ax.text(x + w/2, y + h - 0.035, title, ha="center", va="top", fontsize=title_fs, fontweight="bold")
    start = y + h - 0.11
    step = 0.047
    for i, line in enumerate(lines):
        ax.text(x + 0.02, start - i*step, line, ha="left", va="top", fontsize=text_fs)

def add_arrow(ax, x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=18, linewidth=1.8, color="black"))

def main():
    out = Path("outputs/figures/Raman")
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11})
    fig, ax = plt.subplots(figsize=(15, 7), dpi=300)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    add_box(ax, 0.03, 0.28, 0.25, 0.42, "Raman-derived\nheterogeneity", [
        "• E$_{2g}$ band retained near 1367 cm$^{-1}$",
        "• FWHM ≈ 12.9 cm$^{-1}$; mild asymmetry",
        "• Peak-core / flank segmentation",
        "• Consistent with edge-, strain-, defect-,",
        "  and stacking-perturbed h-BN domains"
    ], title_fs=17, text_fs=11.5)

    add_box(ax, 0.38, 0.41, 0.18, 0.18, "Interfacial\nconsequence", [
        "Distributed local environments",
        "for ion access, polarization,",
        "and charge transfer"
    ], title_fs=15.5, text_fs=11.3)

    add_box(ax, 0.68, 0.69, 0.27, 0.18, "CV response", [
        "• Mixed b-value behavior",
        "• Potential-dependent kinetic domains",
        "• Capacitive + diffusion-influenced current"
    ], title_fs=16.5, text_fs=11.5)

    add_box(ax, 0.68, 0.42, 0.27, 0.18, "GCD response", [
        "• Non-ideal discharge curvature",
        "• Stable-window fitting preferred",
        "• Current-dependent C$_{sp}$ and R$_s$"
    ], title_fs=16.5, text_fs=11.5)

    add_box(ax, 0.68, 0.15, 0.27, 0.18, "EIS response", [
        "• Depressed / dispersive impedance",
        "• CPE-governed relaxation",
        "• Broad time-constant distribution"
    ], title_fs=16.5, text_fs=11.5)

    add_arrow(ax, 0.28, 0.49, 0.38, 0.49)
    add_arrow(ax, 0.56, 0.49, 0.68, 0.78)
    add_arrow(ax, 0.56, 0.49, 0.68, 0.51)
    add_arrow(ax, 0.56, 0.49, 0.68, 0.24)

    ax.text(0.50, 0.06, "Raman-detected structural heterogeneity provides a consistent basis for the mixed CV kinetics, non-ideal GCD behavior, and distributed EIS relaxation of processed h-BN.", ha="center", va="center", fontsize=11.3)

    for ext in ("png","pdf","svg"):
        fig.savefig(out/f"Figure_NS14d_Raman.{ext}", dpi=600 if ext=="png" else None, bbox_inches="tight")
    plt.close(fig)
    print("Done.")

if __name__ == "__main__":
    main()
