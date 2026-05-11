
import matplotlib.pyplot as plt
from cv_utils import extract_anodic_branches, get_common_anodic_grid, ising_segmentation, qkpca_coordinates, processed_dir, set_nc_style, style_axes, ising_cmap

set_nc_style()
out = processed_dir()
branches = extract_anodic_branches()
E_grid, I_grid = get_common_anodic_grid(branches, n_points=500)
ising = ising_segmentation(E_grid, I_grid)
ising.to_csv(out / "CV_Ising_segmentation_data.csv", index=False)

# 4a
fig, ax = plt.subplots(figsize=(6.4, 2.6), dpi=300)
ax.imshow(
    ising["Ising_state"].to_numpy()[None, :],
    aspect="auto",
    interpolation="nearest",
    extent=[ising["Potential_V"].min(), ising["Potential_V"].max(), 0, 1],
    cmap=ising_cmap(),
    vmin=0, vmax=1
)
ax.set_xlabel("Potential V(V)", fontsize=16)
ax.set_ylabel("Ising state", fontsize=15)
ax.set_yticks([0.0, 1.0])
ax.set_yticklabels(["DD", "CD"])
ax.set_title("Ising kinetic segmentation (hBN)", fontsize=15, pad=10)
ax.tick_params(direction="in", top=True, right=True, length=5.5, width=1.2, pad=2)
ax.tick_params(which="minor", direction="in", top=True, right=True, length=3.0, width=0.8)
fig.tight_layout()
fig.savefig(out / "CV_Ising_segmentation.png", bbox_inches="tight")
fig.savefig(out / "CV_Ising_segmentation.pdf", bbox_inches="tight")
plt.close(fig)

# 4b
coords = qkpca_coordinates(branches)
coords.to_csv(out / "CV_QKPCA_coordinates.csv", index=False)
fig, ax = plt.subplots(figsize=(5.6, 4.6), dpi=300)
ax.scatter(coords["QKPCA1"], coords["QKPCA2"], marker="x", s=45, linewidths=1.2)
for _, row in coords.iterrows():
    ax.text(row["QKPCA1"] + 0.012, row["QKPCA2"] + 0.004, f'{int(row["scan_rate_mV_s"])}', fontsize=12)
ax.set_xlabel("Quantum KPCA 1", fontsize=16)
ax.set_ylabel("Quantum KPCA 2", fontsize=16)
ax.set_title("Quantum kernel PCA (hBN)", fontsize=15, pad=8)
style_axes(ax)
fig.tight_layout()
fig.savefig(out / "CV_QKPCA.png", bbox_inches="tight")
fig.savefig(out / "CV_QKPCA.pdf", bbox_inches="tight")
