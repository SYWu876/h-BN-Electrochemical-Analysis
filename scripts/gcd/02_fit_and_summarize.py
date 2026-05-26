
from pathlib import Path
import numpy as np
import pandas as pd
from gcd_utils import (
    JS,
    bootstrap_slope_only_parameters,
    discharge_model,
    fit_slope_only_from_diagnostics,
    plot_fit_panel,
    plot_fit_summary,
    plot_raw_slope_only_summary,
    write_slope_only_markdown_report,
)

ROOT = Path(__file__).resolve().parents[2]
DIAG_DIR = ROOT / "data" / "processed" / "GCD" / "diagnostics"
FIG_DIR = ROOT / "outputs" / "figures" / "GCD"
REPORT_DIR = ROOT / "outputs" / "reports" / "GCD"
GCD_DIR = ROOT / "data" / "processed" / "GCD"
TAB_DIR = ROOT / "data" / "processed" / "GCD" / "tables"


def selected_mask_from_diag(diag: pd.DataFrame) -> np.ndarray:
    if "selected_mask" in diag.columns:
        return diag["selected_mask"].to_numpy(dtype=int) == 1

    mask_columns = [column for column in diag.columns if column.startswith("mask_")]
    if not mask_columns:
        raise ValueError("GCD diagnostics table does not include selected_mask or mask_* columns")
    return diag[mask_columns[0]].to_numpy(dtype=int) == 1


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    final_rows = []
    base_rows = []
    bootstrap_rows = []
    residual_rows = []
    residual_summary_rows = []
    diagnostics_by_j = {}
    slope_masks_by_j = {}
    for j in JS:
        diag = pd.read_csv(DIAG_DIR / f"hBN_GCD_J{j}_processed_diagnostics.csv")
        mask = selected_mask_from_diag(diag)
        diag = diag.copy()
        diag["selected_mask"] = mask.astype(int)
        t_abs = diag["time_s"].to_numpy(dtype=float)

        p_best, slope_metrics, slope_mask = fit_slope_only_from_diagnostics(diag, j)
        diagnostics_by_j[j] = diag
        slope_masks_by_j[j] = slope_mask
        t_slope_abs = t_abs[slope_mask]
        v_slope = diag["raw_V"].to_numpy(dtype=float)[slope_mask]
        boot = bootstrap_slope_only_parameters(
            t_slope_abs,
            v_slope,
            j,
            p_best,
            n_boot=100,
            seed=500 + j,
        )

        ci_lo = np.percentile(boot, 2.5, axis=0)
        ci_hi = np.percentile(boot, 97.5, axis=0)
        mean_boot = np.mean(boot, axis=0)
        median_boot = np.median(boot, axis=0)

        V0, Rs, Csp, A, tau = p_best
        t_slope_rel = t_slope_abs - t_slope_abs[0]
        v_fit = discharge_model(t_slope_rel, j, p_best)
        residual_mV = (v_slope - v_fit) * 1000.0
        for boot_id, boot_params in enumerate(boot):
            bootstrap_rows.append({
                "J_A_g^-1": j,
                "bootstrap_id": boot_id,
                "V0": boot_params[0],
                "Rs": boot_params[1],
                "Csp_F_g^-1": boot_params[2],
                "A": boot_params[3],
                "tau_s": boot_params[4],
                "effective_IR_drop_mV": j * boot_params[1] * 1000.0,
            })
        for time_abs, time_rel, v_exp, v_model, residual in zip(t_slope_abs, t_slope_rel, v_slope, v_fit, residual_mV):
            residual_rows.append({
                "J_A_g^-1": j,
                "time_abs_s": time_abs,
                "time_rel_s": time_rel,
                "V_exp_V": v_exp,
                "V_fit_V": v_model,
                "Residual_mV": residual,
            })
        residual_summary_rows.append({
            "J_A_g^-1": j,
            "Residual_mean_mV": float(np.mean(residual_mV)),
            "Residual_std_mV": float(np.std(residual_mV, ddof=1)),
            "Residual_max_abs_mV": float(np.max(np.abs(residual_mV))),
            "N_points": len(residual_mV),
            "window_start_s": slope_metrics["window_start_s"],
            "window_end_s": slope_metrics["window_end_s"],
        })
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
            "MSE_best": slope_metrics["mse"],
            "slope_V_per_s": slope_metrics["slope"],
            "slope_intercept_V": slope_metrics["intercept"],
            "slope_R2": slope_metrics["r_squared"],
            "slope_window_upper_V": slope_metrics["window_upper_V"],
            "slope_window_lower_V": slope_metrics["window_lower_V"],
            "window_selection_method": slope_metrics["selection_method"],
            "flatness_score": slope_metrics["flatness_score"],
            "mean_stability_cost": slope_metrics["mean_stability_cost"],
            "slope_cv": slope_metrics["slope_cv"],
            "slope_magnitude_fraction": slope_metrics["slope_magnitude_fraction"],
            "linearity_rmse_fraction": slope_metrics["linearity_rmse_fraction"],
            "smoothed_R2": slope_metrics["smoothed_R2"],
            "window_voltage_span_fraction": slope_metrics["window_voltage_span_fraction"],
            "window_point_fraction": slope_metrics["window_point_fraction"],
            "selected_t_start_s": slope_metrics["window_start_s"],
            "selected_t_end_s": slope_metrics["window_end_s"],
            "selected_duration_s": slope_metrics["window_duration_s"],
            "N_selected": slope_metrics["N_slope_points"],
            "deltaV_IR_mV_from_bestRs": j * Rs * 1000.0,
        })
        final_rows.append({
            "J_A_g^-1": j,
            "window_start_s": slope_metrics["window_start_s"],
            "window_end_s": slope_metrics["window_end_s"],
            "retained_window_length_s": slope_metrics["window_duration_s"],
            "mask_points": slope_metrics["N_slope_points"],
            "window_selection_method": slope_metrics["selection_method"],
            "flatness_score": slope_metrics["flatness_score"],
            "V0_base_V": V0,
            "Rs_base_ohm": Rs,
            "effective_IR_drop_base_mV": j * Rs * 1000.0,
            "Csp_base_F_g^-1": Csp,
            "A_base": A,
            "tau_base_s": tau,
            "Csp_CI2.5_F_g^-1": ci_lo[2],
            "Csp_median_F_g^-1": median_boot[2],
            "Csp_CI97.5_F_g^-1": ci_hi[2],
            "IRdrop_CI2.5_mV": j * ci_lo[1] * 1000.0,
            "IRdrop_median_mV": j * median_boot[1] * 1000.0,
            "IRdrop_CI97.5_mV": j * ci_hi[1] * 1000.0,
            "A_CI2.5": ci_lo[3],
            "A_median": median_boot[3],
            "A_CI97.5": ci_hi[3],
            "tau_CI2.5_s": ci_lo[4],
            "tau_median_s": median_boot[4],
            "tau_CI97.5_s": ci_hi[4],
        })
        base_rows.append({
            "J_A_g^-1": j,
            "V0": V0,
            "Rs": Rs,
            "Csp_F_g^-1": Csp,
            "A": A,
            "tau_s": tau,
            "effective_IR_drop_mV": j * Rs * 1000.0,
            "mask_points": slope_metrics["N_slope_points"],
            "window_start_s": slope_metrics["window_start_s"],
            "window_end_s": slope_metrics["window_end_s"],
            "window_selection_method": slope_metrics["selection_method"],
            "flatness_score": slope_metrics["flatness_score"],
        })

        plot_fit_panel(diag, j, p_best, FIG_DIR, slope_mask=slope_mask)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TAB_DIR / "hBN_GCD_fit_summary_refit_auxiliary.csv", index=False)
    summary.to_excel(TAB_DIR / "hBN_GCD_fit_summary_refit_auxiliary.xlsx", index=False)
    final_summary = pd.DataFrame(final_rows)
    final_summary.to_csv(TAB_DIR / "hBN_GCD_bounded_fit_summary_final.csv", index=False)
    final_summary.to_excel(TAB_DIR / "hBN_GCD_bounded_fit_summary_final.xlsx", index=False)
    final_summary.to_csv(TAB_DIR / "hBN_GCD_fit_summary_J1_to_J5.csv", index=False)
    final_summary.to_excel(TAB_DIR / "hBN_GCD_fit_summary_J1_to_J5.xlsx", index=False)
    final_summary.to_csv(GCD_DIR / "hBN_GCD_bounded_fit_summary.csv", index=False)
    pd.DataFrame(base_rows).to_csv(GCD_DIR / "hBN_GCD_base_fit_parameters.csv", index=False)
    pd.DataFrame(bootstrap_rows).to_csv(GCD_DIR / "hBN_GCD_bootstrap_samples.csv", index=False)
    pd.DataFrame({
        "J_A_g^-1": final_summary["J_A_g^-1"],
        "Csp_base_F_g^-1": final_summary["Csp_base_F_g^-1"],
        "Csp_median_F_g^-1": final_summary["Csp_median_F_g^-1"],
        "Csp_CI2.5_F_g^-1": final_summary["Csp_CI2.5_F_g^-1"],
        "Csp_CI97.5_F_g^-1": final_summary["Csp_CI97.5_F_g^-1"],
        "IRdrop_base_mV": final_summary["effective_IR_drop_base_mV"],
        "IRdrop_median_mV": final_summary["IRdrop_median_mV"],
        "IRdrop_CI2.5_mV": final_summary["IRdrop_CI2.5_mV"],
        "IRdrop_CI97.5_mV": final_summary["IRdrop_CI97.5_mV"],
        "A_median": final_summary["A_median"],
        "A_CI2.5": final_summary["A_CI2.5"],
        "A_CI97.5": final_summary["A_CI97.5"],
        "tau_median_s": final_summary["tau_median_s"],
        "tau_CI2.5_s": final_summary["tau_CI2.5_s"],
        "tau_CI97.5_s": final_summary["tau_CI97.5_s"],
    }).to_csv(GCD_DIR / "hBN_GCD_bootstrap_summary.csv", index=False)
    pd.DataFrame(residual_rows).to_csv(GCD_DIR / "hBN_GCD_fit_residuals.csv", index=False)
    pd.DataFrame(residual_summary_rows).to_csv(GCD_DIR / "hBN_GCD_residual_summary.csv", index=False)
    plot_fit_summary(summary, FIG_DIR)
    raw_summary_png, _ = plot_raw_slope_only_summary(diagnostics_by_j, slope_masks_by_j, summary, FIG_DIR)
    report = write_slope_only_markdown_report(summary, raw_summary_png, REPORT_DIR)

    print("Finished slope-only GCD fitting summary export.")
    print(f"Markdown report: {report}")

if __name__ == "__main__":
    main()
