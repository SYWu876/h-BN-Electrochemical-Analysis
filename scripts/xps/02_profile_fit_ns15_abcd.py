
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from utils_xps import (
    REGION_SPECS, table_x_dataframe, load_region, fit_fixed_centers,
    ensure_dir
)

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "XPS" / "HBN.xlsx"
if not RAW.exists():
    candidates = list((ROOT / "data" / "raw" / "XPS").glob("HBN*.xls*"))
    if not candidates:
        raise FileNotFoundError("No raw HBN .xls/.xlsx file found in data/raw/XPS/")
    RAW = candidates[0]

OUT_FIG = ROOT / "outputs" / "figures" / "XPS"
OUT_DATA = ROOT / "data" / "processed" / "XPS"
ensure_dir(OUT_FIG)
ensure_dir(OUT_DATA)

table_df = table_x_dataframe()

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.linewidth": 1.2,
    "axes.labelsize": 18,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 11,
})

panel_info = {
    "B1s": {"panel": "(a)", "title": "B 1s"},
    "N1s": {"panel": "(b)", "title": "N 1s"},
    "C1s": {"panel": "(c)", "title": "C 1s"},
    "O1s": {"panel": "(d)", "title": "O 1s"},
}

for region, spec in REGION_SPECS.items():
    df = load_region(RAW, spec["sheet"])
    x = df["BE"].to_numpy()
    y = df["counts"].to_numpy()
    bg = df["background"].to_numpy()
    net = df["net"].to_numpy()

    plot_window = spec["window"]
    fit_window = spec.get("fit_window", spec["window"])

    plot_mask = (x >= plot_window[0]) & (x <= plot_window[1])
    fit_mask = (x >= fit_window[0]) & (x <= fit_window[1])

    xp, yp, bgp = x[plot_mask], y[plot_mask], bg[plot_mask]
    xf, yf = x[fit_mask], net[fit_mask]

    comps_fit = fit_fixed_centers(xf, yf, spec["centers"], spec["amp0"], spec["sig0"], spec["eta0"], spec["ampmax"], spec["sigmax"])
    comps_plot = fit_fixed_centers(xp, net[plot_mask], spec["centers"], spec["amp0"], spec["sig0"], spec["eta0"], spec["ampmax"], spec["sigmax"])
    fit_total = bgp.copy()
    for comp in comps_plot:
        fit_total += comp["y"]

    # export summary
    rows = []
    for label, center, comp in zip(spec["labels"], spec["centers"], comps_fit):
        rows.append({
            "Region": region,
            "Label": label,
            "Center_eV": center,
            "FWHM_eV": float("nan"),
            "Area_fraction": float(table_df[(table_df.Region == region) & (table_df.Label == label)]["Area_fraction"].iloc[0]),
            "Eta": comp["eta"],
        })
    pd.DataFrame(rows).to_csv(OUT_DATA / f"{region}_profile_fit_summary.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.plot(xp, yp, linewidth=1.8, label="Raw data")
    ax.plot(xp, bgp, linewidth=1.6, label="Background")
    ax.plot(xp, fit_total, linewidth=1.8, label="Total fit")

    for label, comp in zip(spec["labels"], comps_plot):
        ax.plot(xp, bgp + comp["y"], linewidth=1.4, label=label)

    for xc in spec["centers"]:
        ax.axvline(xc, linestyle="--", linewidth=1.2)

    ax.invert_xaxis()
    ax.set_xlim(*spec["plot_xlim"])
    ax.set_xlabel("Binding energy (eV)")
    ax.set_ylabel("Intensity (a.u.)")
    ax.tick_params(which="major", direction="in", length=5, width=1.2, top=True, right=True)
    ax.tick_params(which="minor", direction="in", length=3, width=0.8, top=True, right=True)
    ax.minorticks_on()
    ax.text(-0.12, 1.03, panel_info[region]["panel"], transform=ax.transAxes, fontsize=24, ha="left", va="top")
    ax.text(0.83, 0.94, panel_info[region]["title"], transform=ax.transAxes, fontsize=18, ha="left", va="top")
    ax.legend(loc="upper left", frameon=True)

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    plt.tight_layout()
    png = OUT_FIG / f"Figure_NS15_{region}_overlay_on_background.png"
    pdf = OUT_FIG / f"Figure_NS15_{region}_overlay_on_background.pdf"
    plt.savefig(png, dpi=600, bbox_inches="tight", facecolor="white")
    plt.savefig(pdf, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Saved:", png, pdf)
