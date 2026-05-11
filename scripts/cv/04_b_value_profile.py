
import matplotlib.pyplot as plt
from cv_utils import extract_anodic_branches, get_common_anodic_grid, b_value_profile, processed_dir, set_nc_style, style_axes

set_nc_style()
out = processed_dir()
branches = extract_anodic_branches()
E_grid, I_grid = get_common_anodic_grid(branches)
bdf = b_value_profile(E_grid, I_grid)
bdf.to_csv(out / "CV_b_values.csv", index=False)

fig, ax = plt.subplots(figsize=(6.0, 4.6), dpi=300)
ax.plot(bdf["Potential_E_V"], bdf["b_smooth"], linewidth=1.5, label="b(E)")
ax.axhline(1.0, linestyle="--", linewidth=1.2, label="b=1 (surface)")
ax.axhline(0.5, linestyle="--", linewidth=1.2, label="b=0.5 (diffusion)")
ax.set_xlabel("Potential E(V)", fontsize=16)
ax.set_ylabel("b value", fontsize=16)
style_axes(ax)
ax.set_ylim(min(-0.05, bdf["b_smooth"].min() - 0.05), max(1.1, bdf["b_smooth"].max() + 0.05))
leg = ax.legend(frameon=True, fontsize=12, loc="lower left")
leg.get_frame().set_linewidth(1.0)
fig.tight_layout()
fig.savefig(out / "CV_b_value_profile.png", bbox_inches="tight")
fig.savefig(out / "CV_b_value_profile.pdf", bbox_inches="tight")
