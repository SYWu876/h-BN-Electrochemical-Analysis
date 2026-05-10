
from pathlib import Path
from gcd_utils import JS, load_gcd_pair, extract_discharge_branch, compute_diagnostics, select_stable_window, plot_diagnostic_panels

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "GCD" / "GCD_total.csv"
DIAG_DIR = ROOT / "data" / "processed" / "GCD" / "diagnostics"
FIG_DIR = ROOT / "outputs" / "figures" / "GCD"

DIAGNOSTIC_EXPORT_COLUMNS = {
    1: [
        "time_s", "raw_V", "smoothed_V", "dVdt_V_per_s", "abs_d2Vdt2_V_per_s2",
        "local_residual_V", "curvature_norm", "residual_norm", "stability_cost_ai",
        "IR_drop_candidate", "mask_qaoa_like",
    ],
    2: [
        "time_s", "raw_V", "smoothed_V", "dVdt_V_per_s", "abs_d2Vdt2_V_per_s2",
        "local_residual_V", "curvature_norm", "residual_norm", "stability_cost_ai",
        "IR_drop_candidate", "mask_qaoa_like",
    ],
    3: [
        "time_s", "raw_V", "smoothed_V", "dVdt_V_per_s", "abs_d2Vdt2_V_per_s2",
        "local_residual_V", "curvature_norm", "residual_norm", "stability_cost_ai",
        "IR_drop_candidate", "mask_qaoa_like",
    ],
    4: [
        "time_s", "raw_V", "smoothed_V", "dVdt_V_per_s", "abs_d2Vdt2_V_per_s2",
        "local_residual_V", "curvature_norm", "residual_norm", "stability_cost_ai",
        "mask_qaoa_like_tighter",
    ],
    5: [
        "time_s", "raw_V", "smoothed_V", "dVdt_V_per_s", "abs_d2Vdt2_V_per_s2",
        "local_residual_V", "stability_cost_ai", "IR_drop_candidate", "mask_qaoa_like",
    ],
}


def export_diagnostics(diag, j):
    columns = [column for column in DIAGNOSTIC_EXPORT_COLUMNS[j] if column in diag.columns]
    return diag.loc[:, columns]


def main():
    DIAG_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    for j in JS:
        t_full, v_full = load_gcd_pair(RAW, j)
        t_dis, v_dis = extract_discharge_branch(t_full, v_full)
        diag = compute_diagnostics(t_dis, v_dis)
        diag = select_stable_window(diag, j)
        export_diagnostics(diag, j).to_csv(DIAG_DIR / f"hBN_GCD_J{j}_processed_diagnostics.csv", index=False)
        plot_diagnostic_panels(diag, j, FIG_DIR)
    print("Finished preprocessing, stable-window selection, and diagnostic panel export.")

if __name__ == "__main__":
    main()
