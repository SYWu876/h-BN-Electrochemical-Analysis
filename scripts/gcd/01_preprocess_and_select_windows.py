
from pathlib import Path
from gcd_utils import JS, load_gcd_pair, extract_discharge_branch, compute_diagnostics, select_stable_window, plot_diagnostic_panels

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "GCD" / "GCD_total.csv"
DIAG_DIR = ROOT / "data" / "processed" / "GCD" / "diagnostics"
FIG_DIR = ROOT / "outputs" / "figures" / "GCD"

def main():
    DIAG_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    for j in JS:
        t_full, v_full = load_gcd_pair(RAW, j)
        t_dis, v_dis = extract_discharge_branch(t_full, v_full)
        diag = compute_diagnostics(t_dis, v_dis)
        diag = select_stable_window(diag, j)
        diag.to_csv(DIAG_DIR / f"hBN_GCD_J{j}_processed_diagnostics.csv", index=False)
        plot_diagnostic_panels(diag, j, FIG_DIR)
    print("Finished preprocessing, stable-window selection, and diagnostic panel export.")

if __name__ == "__main__":
    main()
