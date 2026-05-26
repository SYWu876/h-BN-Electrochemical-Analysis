from __future__ import annotations

import numpy as np
import pandas as pd

from cv_utils import (
    dunn_partition,
    extract_anodic_branches,
    get_common_anodic_grid,
    peak_summary,
    processed_dir,
)


def anodic_charge_mC(branch: dict[str, object]) -> float:
    """Return signed anodic charge from I(E), using mA, V, and V/s units."""
    scan_rate_v_s = float(branch["scan_rate_mV_s"]) / 1000.0
    potential_v = np.asarray(branch["E_V"], dtype=float)
    current_mA = np.asarray(branch["I_mA"], dtype=float)
    return float(np.trapezoid(current_mA, potential_v) / scan_rate_v_s)


def build_descriptor_summary() -> pd.DataFrame:
    branches = extract_anodic_branches()
    peak_df = peak_summary()
    _, frac_df = dunn_partition(*get_common_anodic_grid(branches))

    charge_by_rate = {
        float(branch["scan_rate_mV_s"]): anodic_charge_mC(branch)
        for branch in branches
    }
    frac_by_rate = dict(
        zip(
            frac_df["scan_rate_mV_s"].astype(float),
            frac_df["fractional_capacitive_fcap"].astype(float),
        )
    )

    rows = []
    for _, row in peak_df.iterrows():
        scan_rate = float(row["scan_rate_mV_s"])
        total_charge = charge_by_rate[scan_rate]
        capacitive_fraction = frac_by_rate[scan_rate]
        capacitive_charge = total_charge * capacitive_fraction
        diffusion_charge = total_charge - capacitive_charge
        rows.append(
            {
                "Scan_rate_mV_s": int(scan_rate) if scan_rate.is_integer() else scan_rate,
                "Anodic_peak_potential_V_vs_NHE": float(row["peak_potential_V"]),
                "Anodic_peak_current_mA": float(row["peak_current_mA"]),
                "Total_anodic_charge_Qtot_mC": total_charge,
                "Capacitive_charge_Qcap_mC": capacitive_charge,
                "Diffusion_influenced_charge_Qdiff_mC": diffusion_charge,
                "Capacitive_fraction_chi_cap": capacitive_fraction,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    out = processed_dir()
    descriptor_df = build_descriptor_summary()
    descriptor_df.to_csv(out / "hBN_CV_descriptors.csv", index=False)
    print("Saved CV descriptor summary:", out / "hBN_CV_descriptors.csv")


if __name__ == "__main__":
    main()
