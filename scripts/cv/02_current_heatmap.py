
import numpy as np
import matplotlib.pyplot as plt
from cv_utils import load_cv_dataframe, extract_all_curves, processed_dir, set_nc_style, style_axes

set_nc_style()
curves = extract_all_curves(load_cv_dataframe())
out = processed_dir()

# Use absolute current on a common descending grid over the overlap of full scans
E_min = max(c["E_V"].min() for c in curves)
E_max = min(c["E_V"].max() for c in curves)
E_grid = np.linspace(E_min, E_max, 600)
Iabs = np.vstack([np.interp(E_grid, c["E_V"][::-1], np.abs(c["I_mA"])[::-1]) if c["E_V"][0] > c["E_V"][-1]
                  else np.interp(E_grid, c["E_V"], np.abs(c["I_mA"])) for c in curves])

fig, ax = plt.subplots(figsize=(6.3, 4.8), dpi=300)
im = ax.imshow(
    Iabs,
    aspect="auto",
    origin="lower",
    extent=[E_grid.min(), E_grid.max(), 5, 50],
    interpolation="nearest"
)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r"$|I(E)|$ (mA)", fontsize=15)
ax.set_xlabel("Potential E(V)", fontsize=16)
ax.set_ylabel(r"Scan rate $\nu$ (mV/s)", fontsize=16)
style_axes(ax)
fig.tight_layout()
fig.savefig(out / "Figure2b_hBN_current_heatmap.png", bbox_inches="tight")
fig.savefig(out / "Figure2b_hBN_current_heatmap.pdf", bbox_inches="tight")
