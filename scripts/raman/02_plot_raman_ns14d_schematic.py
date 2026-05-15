#!/usr/bin/env python3
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


def wrapped_lines(lines, width):
    rendered = []
    for line in lines:
        parts = wrap(line, width=width, break_long_words=False, break_on_hyphens=False)
        if not parts:
            rendered.append("")
            continue
        rendered.append(f"- {parts[0]}")
        rendered.extend(f"  {part}" for part in parts[1:])
    return rendered


def add_box(ax, x, y, w, h, title, lines, title_fs=12.5, text_fs=8.6, wrap_width=34, bullet=True):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=1.6,
        facecolor="white",
        edgecolor="black",
    )
    ax.add_patch(box)
    ax.text(
        x + w / 2,
        y + h - 0.035,
        title,
        ha="center",
        va="top",
        fontsize=title_fs,
        fontweight="bold",
        linespacing=1.05,
    )

    title_lines = title.count("\n") + 1
    start = y + h - (0.085 + 0.045 * (title_lines - 1))
    step = 0.031
    text_lines = wrapped_lines(lines, wrap_width) if bullet else lines
    for i, line in enumerate(text_lines):
        ax.text(
            x + 0.02,
            start - i * step,
            line,
            ha="left",
            va="top",
            fontsize=text_fs,
            linespacing=1.05,
        )


def add_arrow(ax, x1, y1, x2, y2):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=18,
            linewidth=1.8,
            color="black",
        )
    )


def main():
    out = Path("outputs/figures/Raman")
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, ax = plt.subplots(figsize=(17.5, 7.4), dpi=300)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    add_box(
        ax,
        0.10,
        0.37,
        0.20,
        0.26,
        "Raman-derived\nheterogeneity",
        [
            "E$_{2g}$ near 1367 cm$^{-1}$",
            "FWHM $\\approx$ 12.9 cm$^{-1}$",
            "Core / flank heterogeneity",
        ],
        title_fs=12.3,
        text_fs=7.8,
        wrap_width=80,
    )

    add_box(
        ax,
        0.43,
        0.37,
        0.20,
        0.26,
        "Interfacial\nconsequence",
        [
            "Distributed local environments",
            "for ion access, polarization,",
            "and charge transfer",
        ],
        title_fs=12.3,
        text_fs=8.6,
        wrap_width=29,
        bullet=False,
    )

    add_box(
        ax,
        0.70,
        0.67,
        0.28,
        0.25,
        "CV response",
        [
            "Mixed b-value behavior",
            "Potential-dependent kinetic domains",
            "Capacitive + diffusion-influenced current",
        ],
        title_fs=12.8,
        text_fs=8.4,
        wrap_width=34,
    )

    add_box(
        ax,
        0.70,
        0.375,
        0.28,
        0.25,
        "GCD response",
        [
            "Non-ideal discharge curvature",
            "Stable-window fitting preferred",
            "Current-dependent C$_{sp}$ and R$_s$",
        ],
        title_fs=12.8,
        text_fs=8.4,
        wrap_width=34,
    )

    add_box(
        ax,
        0.70,
        0.08,
        0.28,
        0.25,
        "EIS response",
        [
            "Depressed / dispersive impedance",
            "CPE-governed relaxation",
            "Broad time-constant distribution",
        ],
        title_fs=12.8,
        text_fs=8.4,
        wrap_width=34,
    )

    add_arrow(ax, 0.30, 0.50, 0.43, 0.50)
    add_arrow(ax, 0.63, 0.50, 0.70, 0.79)
    add_arrow(ax, 0.63, 0.50, 0.70, 0.50)
    add_arrow(ax, 0.63, 0.50, 0.70, 0.205)

    ax.text(
        0.50,
        0.055,
        "Raman-detected structural heterogeneity provides a consistent basis for mixed CV kinetics,\n"
        "non-ideal GCD behavior, and distributed EIS relaxation in processed h-BN.",
        ha="center",
        va="center",
        fontsize=9.8,
        linespacing=1.2,
    )

    for ext in ("png", "pdf", "svg"):
        fig.savefig(out / f"Raman_heterogeneity_schematic.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
