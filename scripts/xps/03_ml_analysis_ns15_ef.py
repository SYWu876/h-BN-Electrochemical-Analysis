
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, leaves_list
from utils_xps import build_descriptor_matrix, ensure_dir, MARKER_MAP

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

ml_df = build_descriptor_matrix(RAW)

# PCA
X = StandardScaler().fit_transform(ml_df[["DeltaE_eV", "FWHM_eV", "Area_fraction", "Eta"]].to_numpy())
pca = PCA(n_components=2)
scores = pca.fit_transform(X)
ml_df["PC1"] = scores[:, 0]
ml_df["PC2"] = scores[:, 1]
ml_df.to_csv(OUT_DATA / "Table_X1_corrected_13peak_descriptor_matrix.csv", index=False)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.linewidth": 2.0,
    "axes.labelsize": 22,
    "xtick.labelsize": 17,
    "ytick.labelsize": 17,
    "legend.fontsize": 14,
})

# NS15(e) symbols only
fig, ax = plt.subplots(figsize=(8.8, 6.8))
for region in ["B1s", "N1s", "C1s", "O1s"]:
    dfr = ml_df[ml_df["Region"] == region]
    ax.scatter(dfr["PC1"], dfr["PC2"], s=150, marker=MARKER_MAP[region], edgecolors="black", linewidths=1.4, zorder=3)
ax.axhline(0, linewidth=1.2, zorder=1)
ax.axvline(0, linewidth=1.2, zorder=1)
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
ax.tick_params(which="major", direction="in", length=6, width=1.8, top=True, right=True)
ax.tick_params(which="minor", direction="in", length=3.5, width=1.2, top=True, right=True)
ax.minorticks_on()
for spine in ax.spines.values():
    spine.set_linewidth(2.0)
ax.text(-0.14, 1.04, "(e)", transform=ax.transAxes, fontsize=30, ha="left", va="top")
plt.tight_layout()
plt.savefig(OUT_FIG / "Figure_NS15e_XPS_PCA_corrected_13peaks_symbols_only.png", dpi=600, bbox_inches="tight", facecolor="white")
plt.savefig(OUT_FIG / "Figure_NS15e_XPS_PCA_corrected_13peaks_symbols_only.pdf", dpi=600, bbox_inches="tight", facecolor="white")
plt.close(fig)

# NS15(e) labeled
fig, ax = plt.subplots(figsize=(9.0, 6.9))
for region in ["B1s", "N1s", "C1s", "O1s"]:
    dfr = ml_df[ml_df["Region"] == region]
    ax.scatter(dfr["PC1"], dfr["PC2"], s=160, marker=MARKER_MAP[region], edgecolors="black", linewidths=1.4, zorder=3)
offsets = {
    "B-def": (0.05, -0.12), "B–N": (0.05, 0.08), "B-ox": (0.05, -0.10),
    "N-def": (0.05, -0.12), "N–B": (0.05, 0.08), "N-ox": (0.05, 0.08),
    "C–C": (0.05, 0.08), "C–O": (0.05, -0.12), "C=O": (0.05, -0.10), "O–C=O": (0.05, 0.08),
    "O-lat": (0.05, -0.12), "O-ads": (0.05, 0.08), "H2O/O": (0.05, 0.08),
}
for _, row in ml_df.iterrows():
    dx, dy = offsets.get(row["Label"], (0.05, 0.05))
    ax.text(row["PC1"] + dx, row["PC2"] + dy, row["Label"], fontsize=13)
ax.axhline(0, linewidth=1.2, zorder=1)
ax.axvline(0, linewidth=1.2, zorder=1)
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
ax.tick_params(which="major", direction="in", length=6, width=1.8, top=True, right=True)
ax.tick_params(which="minor", direction="in", length=3.5, width=1.2, top=True, right=True)
ax.minorticks_on()
for spine in ax.spines.values():
    spine.set_linewidth(2.0)
ax.text(-0.14, 1.04, "(e)", transform=ax.transAxes, fontsize=30, ha="left", va="top")
plt.tight_layout()
plt.savefig(OUT_FIG / "Figure_NS15e_XPS_PCA_corrected_13peaks_with_labels.png", dpi=600, bbox_inches="tight", facecolor="white")
plt.savefig(OUT_FIG / "Figure_NS15e_XPS_PCA_corrected_13peaks_with_labels.pdf", dpi=600, bbox_inches="tight", facecolor="white")
plt.close(fig)

# NS15(f) clustered heat map
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.linewidth": 2.0,
    "axes.labelsize": 20,
    "xtick.labelsize": 14,
    "ytick.labelsize": 13,
})
row_labels = (ml_df["Region"] + " " + ml_df["Label"]).tolist()
feat_names = ["Area_fraction", "Eta", "DeltaE_eV", "FWHM_eV"]
Xz = StandardScaler().fit_transform(ml_df[feat_names].to_numpy())
row_order = leaves_list(linkage(Xz, method="ward"))
col_order = leaves_list(linkage(Xz.T, method="ward"))
Xz_ord = Xz[row_order][:, col_order]
row_labels_ord = [row_labels[i] for i in row_order]
feat_names_ord = [feat_names[i] for i in col_order]
fig, ax = plt.subplots(figsize=(8.8, 7.4))
im = ax.imshow(Xz_ord, aspect="auto", interpolation="nearest")
ax.set_xticks(np.arange(len(feat_names_ord)))
ax.set_xticklabels(feat_names_ord, rotation=25, ha="right")
ax.set_yticks(np.arange(len(row_labels_ord)))
ax.set_yticklabels(row_labels_ord)
ax.tick_params(which="major", direction="in", length=5, width=1.6, top=True, right=True)
for spine in ax.spines.values():
    spine.set_linewidth(2.0)
ax.text(-0.18, 1.04, "(f)", transform=ax.transAxes, fontsize=30, ha="left", va="top")
cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cbar.set_label("Standardized descriptor value", fontsize=16)
cbar.ax.tick_params(labelsize=12, width=1.2, length=4, direction="in")
plt.tight_layout()
plt.savefig(OUT_FIG / "Figure_NS15f_XPS_clustered_heatmap_corrected_13peaks.png", dpi=600, bbox_inches="tight", facecolor="white")
plt.savefig(OUT_FIG / "Figure_NS15f_XPS_clustered_heatmap_corrected_13peaks.pdf", dpi=600, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("Saved processed data and figures to:")
print(OUT_DATA)
print(OUT_FIG)
