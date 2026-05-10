
from pathlib import Path
import pandas as pd
from gcd_utils import JS, fit_selected_segment, bootstrap_parameters, plot_fit_panel, plot_fit_summary

ROOT = Path(__file__).resolve().parents[2]
DIAG_DIR = ROOT / "data" / "processed" / "GCD" / "diagnostics"
FIG_DIR = ROOT / "outputs" / "figures" / "GCD"
TAB_DIR = ROOT / "data" / "processed" / "GCD" / "tables"

def main():
    summary_rows = []
    for j in JS:
        diag = pd.read_csv(DIAG_DIR / f"hBN_GCD_J{j}_processed_diagnostics.csv")
        mask = diag["selected_mask"].to_numpy(dtype=int) == 1
        t_sel = diag["time_s"].to_numpy(dtype=float)[mask]
        v_sel = diag["raw_V"].to_numpy(dtype=float)[mask]

        p_best, mse_best = fit_selected_segment(t_sel, v_sel, j, n_starts=12, seed=100 + j)
        boot = bootstrap_parameters(t_sel, v_sel, j, p_best, n_boot=50, seed=500 + j)

        import numpy as np
        ci_lo = np.percentile(boot, 2.5, axis=0)
        ci_hi = np.percentile(boot, 97.5, axis=0)
        mean_boot = np.mean(boot, axis=0)

        V0, Rs, Csp, A, tau = p_best
        summary_rows.append({
            "J_A_g": j,
            "V0_best_V": V0,
            "V0_boot_mean_V": mean_boot[0],
            "V0_CI_low_V": ci_lo[0],
            "V0_CI_high_V": ci_hi[0],
            "Rs_best_ohm_g": Rs,
            "Rs_boot_mean_ohm_g": mean_boot[1],
            "Rs_CI_low_ohm_g": ci_lo[1],
            "Rs_CI_high_ohm_g": ci_hi[1],
            "Csp_best_F_g": Csp,
            "Csp_boot_mean_F_g": mean_boot[2],
            "Csp_CI_low_F_g": ci_lo[2],
            "Csp_CI_high_F_g": ci_hi[2],
            "A_best_V": A,
            "A_boot_mean_V": mean_boot[3],
            "A_CI_low_V": ci_lo[3],
            "A_CI_high_V": ci_hi[3],
            "tau_best_s": tau,
            "tau_boot_mean_s": mean_boot[4],
            "tau_CI_low_s": ci_lo[4],
            "tau_CI_high_s": ci_hi[4],
            "MSE_best": mse_best,
            "selected_t_start_s": t_sel[0],
            "selected_t_end_s": t_sel[-1],
            "selected_duration_s": t_sel[-1] - t_sel[0],
            "N_selected": len(t_sel),
            "deltaV_IR_mV_from_bestRs": j * Rs * 1000.0,
        })

        plot_fit_panel(diag, j, p_best, FIG_DIR)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TAB_DIR / "hBN_GCD_fit_summary_J1_to_J5.csv", index=False)
    summary.to_excel(TAB_DIR / "hBN_GCD_fit_summary_J1_to_J5.xlsx", index=False)

    plot_fit_summary(summary, FIG_DIR)
    print("Finished bounded fitting, bootstrap CI, and Figure 6 export.")

if __name__ == "__main__":
    main()
