from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import AutoMinorLocator
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


DOMAINS = ["TEM", "Raman", "XPS", "CV", "GCD", "EIS"]
EXPECTED_DESCRIPTORS_PER_DOMAIN = 4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build conservative cross-domain descriptor tables for h-BN."
    )
    parser.add_argument("--repo-root", type=Path, default=Path("."), help="Repository root.")
    parser.add_argument(
        "--make-figures",
        action="store_true",
        help="Also write local heatmap/PCA figure exports under outputs/figures/integrated/.",
    )
    return parser.parse_args()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def add_descriptor(
    rows: list[dict[str, object]],
    domain: str,
    descriptor: str,
    value: float,
    unit: str = "",
    source_path: str = "",
    note: str = "",
) -> None:
    rows.append(
        {
            "domain": domain,
            "descriptor": descriptor,
            "value": float(value),
            "unit": unit,
            "source_path": source_path,
            "note": note,
        }
    )


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    return float(np.average(values.to_numpy(dtype=float), weights=weights.to_numpy(dtype=float)))


def collect_tem(root: Path, rows: list[dict[str, object]]) -> None:
    rel = "data/processed/TEM/TEM_descriptors.csv"
    df = pd.read_csv(root / rel)
    weights = df["geometric_patch_weight_wi"]
    add_descriptor(rows, "TEM", "mean_LOI_0_100", df["LOI_0_100"].mean(), "0-100", rel)
    add_descriptor(rows, "TEM", "weighted_LOI_0_100", weighted_mean(df["LOI_0_100"], weights), "0-100", rel)
    add_descriptor(rows, "TEM", "max_patch_weight", weights.max(), "-", rel)
    add_descriptor(rows, "TEM", "weighted_Delta_q_FFT", weighted_mean(df["Delta_q_FFT"], weights), "a.u.", rel)


def collect_raman(root: Path, rows: list[dict[str, object]]) -> None:
    summary_rel = "data/processed/Raman/raman_fit_summary.csv"
    local_rel = "data/processed/Raman/raman_local_descriptors_and_segmentation.csv"
    summary = pd.read_csv(root / summary_rel).iloc[0]
    local = pd.read_csv(root / local_rel)
    add_descriptor(rows, "Raman", "fit_center_cm-1", summary["fit_center_cm-1"], "cm^-1", summary_rel)
    add_descriptor(rows, "Raman", "fwhm_cm-1", summary["fwhm_cm-1"], "cm^-1", summary_rel)
    add_descriptor(
        rows,
        "Raman",
        "asymmetry_ratio_right_over_left",
        summary["asymmetry_ratio_right_over_left"],
        "-",
        summary_rel,
    )
    add_descriptor(rows, "Raman", "mean_abs_d2I_domega2", local["d2I_domega2"].abs().mean(), "a.u.", local_rel)


def collect_xps(root: Path, rows: list[dict[str, object]]) -> None:
    matrix_rel = "data/processed/XPS/XPS_corrected_13peak_descriptor_matrix.csv"
    matrix = pd.read_csv(root / matrix_rel)
    add_descriptor(rows, "XPS", "defect_area_fraction_sum", matrix[matrix["Label"].str.contains("def", case=False)]["Area_fraction"].sum(), "-", matrix_rel)
    add_descriptor(rows, "XPS", "oxidized_area_fraction_sum", matrix[matrix["Label"].str.contains("ox|O-", case=False, regex=True)]["Area_fraction"].sum(), "-", matrix_rel)
    add_descriptor(rows, "XPS", "pc1_mean", matrix["PC1"].mean(), "-", matrix_rel)
    add_descriptor(rows, "XPS", "pc2_span", matrix["PC2"].max() - matrix["PC2"].min(), "-", matrix_rel)


def collect_cv(root: Path, rows: list[dict[str, object]]) -> None:
    peaks_rel = "data/processed/CV/CV_peak_summary.csv"
    b_rel = "data/processed/CV/CV_b_values.csv"
    frac_rel = "data/processed/CV/CV_fractional_capacitive.csv"
    qkpca_rel = "data/processed/CV/CV_QKPCA_coordinates.csv"
    peaks = pd.read_csv(root / peaks_rel)
    bvals = pd.read_csv(root / b_rel)
    frac = pd.read_csv(root / frac_rel)
    qkpca = pd.read_csv(root / qkpca_rel)
    slope, _ = np.polyfit(peaks["sqrt_scan_rate_sqrt_mV_s"], peaks["peak_current_mA"], 1)
    add_descriptor(rows, "CV", "peak_current_vs_sqrt_scan_slope", slope, "mA/(mV s^-1)^0.5", peaks_rel)
    add_descriptor(rows, "CV", "mean_b_smooth", bvals["b_smooth"].mean(), "-", b_rel)
    add_descriptor(rows, "CV", "max_capacitive_fraction", frac["fractional_capacitive_fcap"].max(), "-", frac_rel)
    add_descriptor(rows, "CV", "QKPCA1_span", qkpca["QKPCA1"].max() - qkpca["QKPCA1"].min(), "-", qkpca_rel)


def collect_gcd(root: Path, rows: list[dict[str, object]]) -> None:
    current_rel = "data/processed/GCD/tables/hBN_GCD_fit_summary_J1_to_J5.csv"
    reference_rel = "data/processed/GCD/tables/hBN_GCD_bounded_fit_summary_final.csv"
    rel = current_rel if (root / current_rel).exists() else reference_rel
    df = pd.read_csv(root / rel)
    csp = df["Csp_base_F_g^-1"]
    rs = df["Rs_base_ohm"]
    ir_drop = df["effective_IR_drop_base_mV"]
    add_descriptor(rows, "GCD", "Csp_mean_F_g-1", csp.mean(), "F g^-1", rel)
    add_descriptor(rows, "GCD", "Csp_max_F_g-1", csp.max(), "F g^-1", rel)
    add_descriptor(rows, "GCD", "Rs_mean_ohm", rs.mean(), "Ohm", rel)
    add_descriptor(rows, "GCD", "IRdrop_max_mV", ir_drop.max(), "mV", rel)


def collect_eis(root: Path, rows: list[dict[str, object]]) -> None:
    params_rel = "data/processed/EIS/classical_fit/hBN_EIS_final_fit_parameters.csv"
    deviations_rel = "data/processed/EIS/quantum_branches/hBN_parameter_deviations.csv"
    params = pd.read_csv(root / params_rel)
    deviations = pd.read_csv(root / deviations_rel)
    lookup = params.set_index("Parameter")["Value"]
    add_descriptor(rows, "EIS", "Rs_ohm", lookup["Rs"], "Ohm", params_rel)
    add_descriptor(rows, "EIS", "R1_ohm", lookup["R1"], "Ohm", params_rel)
    add_descriptor(rows, "EIS", "alpha1", lookup["alpha1"], "-", params_rel)
    add_descriptor(rows, "EIS", "max_discrete_relative_deviation_percent", deviations["Relative_deviation_discrete_percent"].max(), "%", deviations_rel)


def build_long_table(root: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    collect_tem(root, rows)
    collect_raman(root, rows)
    collect_xps(root, rows)
    collect_cv(root, rows)
    collect_gcd(root, rows)
    collect_eis(root, rows)
    df = pd.DataFrame(rows)
    df["domain"] = pd.Categorical(df["domain"], categories=DOMAINS, ordered=True)
    return df.sort_values(["domain", "descriptor"]).reset_index(drop=True)


def build_heatmap_matrix(long_df: pd.DataFrame) -> pd.DataFrame:
    domains = []
    values = []
    present_domains = set(long_df["domain"].astype(str))
    missing_domains = [domain for domain in DOMAINS if domain not in present_domains]
    if missing_domains:
        raise ValueError(
            "Missing cross-domain descriptor rows for domain(s): "
            + ", ".join(missing_domains)
            + ". Rebuild or inspect the corresponding processed input tables."
        )

    for domain in DOMAINS:
        sub = long_df[long_df["domain"] == domain].sort_values("descriptor")
        count = len(sub)
        if count != EXPECTED_DESCRIPTORS_PER_DOMAIN:
            source_paths = ", ".join(sorted(str(path) for path in sub["source_path"].dropna().unique()))
            source_hint = f" Source paths present: {source_paths}." if source_paths else ""
            raise ValueError(
                f"Expected {EXPECTED_DESCRIPTORS_PER_DOMAIN} descriptors for domain {domain}, "
                f"found {count}.{source_hint}"
            )
        domain_values = sub["value"].to_numpy(dtype=float)
        if not np.all(np.isfinite(domain_values)):
            raise ValueError(f"Non-finite descriptor value found for domain {domain}.")
        domains.append(domain)
        values.append(domain_values)
    matrix = np.vstack(values)
    scaled = StandardScaler().fit_transform(matrix)
    return pd.DataFrame(scaled, index=domains, columns=[f"descriptor_{i+1}" for i in range(matrix.shape[1])])


def build_pca_tables(heatmap: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    n_components = min(2, heatmap.shape[0], heatmap.shape[1])
    pca = PCA(n_components=n_components)
    projection = pca.fit_transform(heatmap.to_numpy(dtype=float))
    projection_columns = [f"PC{i+1}" for i in range(n_components)]
    projection_df = pd.DataFrame(projection, columns=projection_columns)
    projection_df.insert(0, "domain", heatmap.index)
    projection_df["explained_variance_ratio_PC1"] = pca.explained_variance_ratio_[0]
    projection_df["explained_variance_ratio_PC2"] = pca.explained_variance_ratio_[1] if n_components > 1 else np.nan
    loadings_df = pd.DataFrame(
        pca.components_.T,
        columns=projection_columns,
    )
    loadings_df.insert(0, "descriptor_axis", heatmap.columns)
    return projection_df, loadings_df


def style_axes(ax: plt.Axes) -> None:
    ax.tick_params(direction="in", top=True, right=True, length=5, width=1.0)
    ax.tick_params(which="minor", direction="in", top=True, right=True, length=2.5, width=0.8)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def make_figures(heatmap: pd.DataFrame, pca_projection: pd.DataFrame, out_dir: Path) -> None:
    ensure_dir(out_dir)
    plt.rcParams.update({"font.family": "serif", "font.size": 11, "axes.linewidth": 1.1})

    fig, ax = plt.subplots(figsize=(6.6, 4.2), dpi=300)
    im = ax.imshow(heatmap.to_numpy(dtype=float), aspect="auto", cmap="coolwarm", vmin=-2.0, vmax=2.0)
    ax.set_yticks(np.arange(len(heatmap.index)))
    ax.set_yticklabels(heatmap.index)
    ax.set_xticks(np.arange(len(heatmap.columns)))
    ax.set_xticklabels(heatmap.columns, rotation=45, ha="right")
    ax.set_title("Cross-domain descriptor matrix (z-scored)")
    fig.colorbar(im, ax=ax, label="z-score")
    fig.tight_layout()
    fig.savefig(out_dir / "hBN_cross_domain_heatmap.png", bbox_inches="tight")
    fig.savefig(out_dir / "hBN_cross_domain_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    has_pc2 = "PC2" in pca_projection.columns and pca_projection["PC2"].notna().any()
    y_values = pca_projection["PC2"] if has_pc2 else np.zeros(len(pca_projection))

    fig, ax = plt.subplots(figsize=(5.2, 4.4), dpi=300)
    ax.scatter(pca_projection["PC1"], y_values, s=70)
    for _, row in pca_projection.iterrows():
        y = row["PC2"] if has_pc2 else 0.0
        ax.text(row["PC1"], y, f"  {row['domain']}", va="center")
    ax.set_xlabel("PC1 (exploratory)")
    ax.set_ylabel("PC2 (exploratory)" if has_pc2 else "PC2 unavailable")
    ax.set_title("Cross-domain exploratory PCA")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(out_dir / "hBN_cross_domain_pca.png", bbox_inches="tight")
    fig.savefig(out_dir / "hBN_cross_domain_pca.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    root = args.repo_root.resolve()
    out_dir = root / "data" / "processed" / "integrated"
    ensure_dir(out_dir)

    long_df = build_long_table(root)
    heatmap = build_heatmap_matrix(long_df)
    pca_projection, pca_loadings = build_pca_tables(heatmap)

    long_df.to_csv(out_dir / "hBN_cross_domain_descriptor_long.csv", index=False, lineterminator="\n")
    heatmap.to_csv(out_dir / "hBN_cross_domain_heatmap_matrix.csv", index_label="domain", lineterminator="\n")
    pca_projection.to_csv(out_dir / "hBN_cross_domain_pca_projection.csv", index=False, lineterminator="\n")
    pca_loadings.to_csv(out_dir / "hBN_cross_domain_pca_loadings.csv", index=False, lineterminator="\n")

    if args.make_figures:
        make_figures(heatmap, pca_projection, root / "outputs" / "figures" / "integrated")

    print(f"Wrote integrated descriptor outputs to {out_dir}")


if __name__ == "__main__":
    main()
