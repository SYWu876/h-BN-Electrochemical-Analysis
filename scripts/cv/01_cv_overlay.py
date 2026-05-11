
from pathlib import Path
import matplotlib.pyplot as plt
from cv_utils import load_cv_dataframe, extract_all_curves, processed_dir, set_nc_style, style_axes

set_nc_style()
curves = extract_all_curves(load_cv_dataframe())
out = processed_dir()

fig, ax = plt.subplots(figsize=(6.3, 5.1), dpi=300)
for c in curves:
    ax.plot(c["E_V"], c["I_mA"], linewidth=1.5, label=fr'{int(c["scan_rate_mV_s"])} mV s$^{{-1}}$')

ax.set_xlabel("Potential E(V)", fontsize=16)
ax.set_ylabel("Current I (mA)", fontsize=16)
style_axes(ax)
leg = ax.legend(title="Scan rate", frameon=True, fontsize=11, title_fontsize=11, loc="upper left")
leg.get_frame().set_linewidth(1.0)
fig.tight_layout()
fig.savefig(out / "CV_overlay.png", bbox_inches="tight")
fig.savefig(out / "CV_overlay.pdf", bbox_inches="tight")
