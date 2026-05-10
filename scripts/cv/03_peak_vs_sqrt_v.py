
import numpy as np
import matplotlib.pyplot as plt
from cv_utils import peak_summary, processed_dir, set_nc_style, style_axes, linear_fit

set_nc_style()
out = processed_dir()
df = peak_summary()
df.to_csv(out / "Figure2c_peak_summary.csv", index=False)

x = df["sqrt_scan_rate_sqrt_mV_s"].to_numpy()
y = df["peak_current_mA"].to_numpy()
coef, y_fit, r2 = linear_fit(x, y)
xfit = np.linspace(x.min(), x.max(), 300)
yline = np.polyval(coef, xfit)

fig, ax = plt.subplots(figsize=(5.8, 4.7), dpi=300)
ax.plot(x, y, linestyle="none", marker="x", markersize=9, markeredgewidth=1.6, label="data")
ax.plot(xfit, yline, linewidth=1.6, label=fr"fit (R$^2$={r2:.3f})")
ax.set_xlabel(r"$\sqrt{\nu}\ \mathrm{(mV/s)^{1/2}}$", fontsize=16)
ax.set_ylabel(r"$I_{\mathrm{peak}}\ \mathrm{(mA)}$", fontsize=16)
style_axes(ax)
leg = ax.legend(frameon=True, fontsize=12, loc="upper left")
leg.get_frame().set_linewidth(1.0)
fig.tight_layout()
fig.savefig(out / "Figure2c_hBN_peak_vs_sqrtv.png", bbox_inches="tight")
fig.savefig(out / "Figure2c_hBN_peak_vs_sqrtv.pdf", bbox_inches="tight")
