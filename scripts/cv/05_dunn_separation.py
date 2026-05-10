
import matplotlib.pyplot as plt
from cv_utils import extract_anodic_branches, get_common_anodic_grid, dunn_partition, processed_dir, set_nc_style, style_axes

set_nc_style()
out = processed_dir()
branches = extract_anodic_branches()
E_grid, I_grid = get_common_anodic_grid(branches)
prof, frac = dunn_partition(E_grid, I_grid)
prof.to_csv(out / "Figure3_Dunn_partition_data.csv", index=False)
frac.to_csv(out / "Figure3c_fractional_capacitive.csv", index=False)

# 3a
fig, ax = plt.subplots(figsize=(5.6, 4.0), dpi=300)
ax.plot(prof["Potential_E_V"], prof["k1_mA_per_mV_per_s"], linewidth=1.5)
ax.set_xlabel("Potential E(V)", fontsize=16)
ax.set_ylabel(r"$k_1(E)$ (mA (mV s$^{-1}$)$^{-1}$)", fontsize=15)
style_axes(ax)
fig.tight_layout()
fig.savefig(out / "Figure3a_hBN_k1.png", bbox_inches="tight")
fig.savefig(out / "Figure3a_hBN_k1.pdf", bbox_inches="tight")
plt.close(fig)

# 3b
fig, ax = plt.subplots(figsize=(5.6, 4.0), dpi=300)
ax.plot(prof["Potential_E_V"], prof["k2_mA_per_sqrt_mV_per_s"], linewidth=1.5)
ax.set_xlabel("Potential E(V)", fontsize=16)
ax.set_ylabel(r"$k_2(E)$ (mA (mV s$^{-1}$)$^{-1/2}$)", fontsize=15)
style_axes(ax)
fig.tight_layout()
fig.savefig(out / "Figure3b_hBN_k2.png", bbox_inches="tight")
fig.savefig(out / "Figure3b_hBN_k2.pdf", bbox_inches="tight")
plt.close(fig)

# 3c
fig, ax = plt.subplots(figsize=(4.8, 4.0), dpi=300)
ax.plot(frac["scan_rate_mV_s"], frac["fractional_capacitive_fcap"], marker="o", linewidth=1.5)
ax.set_xlabel(r"Scan rate $\nu$ (mV s$^{-1}$)", fontsize=16)
ax.set_ylabel(r"Fractional capacitive $f_{\mathrm{cap}}$", fontsize=15)
style_axes(ax)
fig.tight_layout()
fig.savefig(out / "Figure3c_hBN_fcap.png", bbox_inches="tight")
fig.savefig(out / "Figure3c_hBN_fcap.pdf", bbox_inches="tight")
