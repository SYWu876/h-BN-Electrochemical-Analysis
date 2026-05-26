from __future__ import annotations

import csv
import importlib.util
import json
import os
import py_compile
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_TIMEOUT_SECONDS = 120
FULL_PIPELINE_TIMEOUT_SECONDS = 360


def load_integrated_module():
    script = ROOT / "scripts" / "integrated" / "00_build_cross_domain_evidence.py"
    spec = importlib.util.spec_from_file_location("build_cross_domain_evidence", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_cv_current_heatmap_module():
    script_dir = ROOT / "scripts" / "cv"
    script = script_dir / "02_current_heatmap.py"
    sys.path.insert(0, str(script_dir))
    try:
        spec = importlib.util.spec_from_file_location("cv_current_heatmap", script)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


def test_readme_and_citation_contact_email_match() -> None:
    readme_text = (ROOT / "README.md").read_text(encoding="utf-8")
    citation_text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")

    readme_match = re.search(r"^Email:\s*([^\s]+)\s*$", readme_text, flags=re.MULTILINE)
    citation_emails = re.findall(r"^\s*email:\s*([^\s]+)\s*$", citation_text, flags=re.MULTILINE)

    assert readme_match is not None, "README correspondence email was not found"
    assert citation_emails, "No author email entries were found in CITATION.cff"
    assert readme_match.group(1) in citation_emails, (
        f"README email {readme_match.group(1)!r} was not found in CITATION emails {citation_emails!r}"
    )
    assert "scripts/eis/01_eis_classical_anchor_fit.py" in readme_text
    assert "scripts/eis/02_eis_quantum_comparison_from_anchor.py" in readme_text
    assert "scripts/eis/03_eis_surrogate_qaoa_landscape.py" in readme_text
    assert "scripts/eis/04_eis_shared_objective_full_pipeline.py" in readme_text
    assert "scripts/cv/" in readme_text
    assert "scripts/gcd/" in readme_text
    assert "scripts/raman/" in readme_text
    assert "scripts/tem/" in readme_text
    assert "scripts/xps/" in readme_text
    assert "scripts/integrated/" in readme_text
    assert "All EIS surrogate slices used in the manuscript can be regenerated" in readme_text
    assert "R1_Q1_slice.csv" in readme_text
    assert "Rs_alpha1_slice.csv" in readme_text


def test_analysis_scripts_compile(tmp_path: Path) -> None:
    for script in sorted((ROOT / "scripts").rglob("*.py")):
        pyc_name = f"{script.relative_to(ROOT).as_posix().replace('/', '__')}.pyc"
        py_compile.compile(str(script), cfile=str(tmp_path / pyc_name), doraise=True)


def test_integrated_domain_archive_files_are_present() -> None:
    expected_paths = [
        "scripts/qc_circuit/README.md",
        "scripts/qc_circuit/generate_gcd_eis_qc_circuits.ipynb",
        "scripts/cv/00_run_all_cv_analysis.py",
        "scripts/eis/01_eis_classical_anchor_fit.py",
        "scripts/eis/02_eis_quantum_comparison_from_anchor.py",
        "scripts/eis/03_eis_surrogate_qaoa_landscape.py",
        "scripts/eis/04_eis_shared_objective_full_pipeline.py",
        "scripts/gcd/00_run_all_gcd_analysis.py",
        "scripts/integrated/00_build_cross_domain_evidence.py",
        "scripts/raman/01_raman_analysis_pipeline.py",
        "scripts/tem/00_tem_patch_ensemble_analysis.py",
        "scripts/xps/04_run_xps_pipeline.py",
        "data/processed/CV/CV_peak_summary.csv",
        "data/processed/GCD/diagnostics/hBN_GCD_J1_processed_diagnostics.csv",
        "data/processed/GCD/tables/hBN_GCD_fit_summary_J1_to_J5.csv",
        "data/processed/GCD/tables/hBN_GCD_bounded_fit_summary_final.csv",
        "data/processed/integrated/hBN_cross_domain_descriptor_long.csv",
        "data/processed/integrated/hBN_cross_domain_heatmap_matrix.csv",
        "data/processed/integrated/hBN_cross_domain_pca_projection.csv",
        "data/processed/integrated/hBN_cross_domain_pca_loadings.csv",
        "data/processed/EIS/classical_fit/hBN_EIS_final_fit_parameters.csv",
        "data/processed/EIS/quantum_branches/hBN_parameter_deviations.csv",
        "data/processed/CV/CV_b_values.csv",
        "data/processed/CV/CV_fractional_capacitive.csv",
        "data/processed/CV/CV_QKPCA_coordinates.csv",
        "data/processed/Raman/raman_fit_summary.csv",
        "data/processed/Raman/raman_local_descriptors_and_segmentation.csv",
        "data/processed/TEM/TEM_descriptors.csv",
        "data/processed/XPS/XPS_13_peak_assignments.csv",
        "data/processed/XPS/XPS_corrected_13peak_descriptor_matrix.csv",
        "data/raw/TEM/roi_boxes_template.csv",
        "data/raw/XPS/HBN.xlsx",
        "docs/Note_S7_Raman_GitHub.md",
        "docs/README_TEM_reproducibility.md",
        "docs/XPS_REPRODUCIBILITY_NOTES.md",
    ]
    for path in expected_paths:
        artifact = ROOT / path
        assert artifact.exists(), path
        assert artifact.stat().st_size > 0, path

    with (ROOT / "data" / "raw" / "TEM" / "roi_boxes_template.csv").open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        assert reader.fieldnames == ["roi", "x", "y", "width", "height"]
        roi_rows = list(reader)
    assert [row["roi"] for row in roi_rows] == ["ROI-1", "ROI-2", "ROI-3", "ROI-4", "ROI-5"]
    roi_boxes = [(int(row["x"]), int(row["y"]), int(row["width"]), int(row["height"])) for row in roi_rows]
    assert len(set(roi_boxes)) == 5
    assert all(x > 0 and y > 0 and width > 0 and height > 0 for x, y, width, height in roi_boxes)

    with (ROOT / "data" / "processed" / "TEM" / "TEM_descriptors.csv").open(newline="", encoding="utf-8") as f:
        assert next(csv.reader(f)) == [
            "ROI",
            "P_FFT",
            "A_off",
            "R_bp",
            "Delta_q_FFT",
            "LOI_0_100",
            "geometric_patch_weight_wi",
        ]

    with (ROOT / "data" / "processed" / "TEM" / "TEM_descriptors_normalized.csv").open(
        newline="", encoding="utf-8"
    ) as f:
        assert next(csv.reader(f)) == [
            "ROI",
            "P_FFT_norm",
            "A_off_norm",
            "R_bp_norm",
            "Delta_q_FFT_norm_raw",
            "Delta_q_FFT_ordering_tendency_inverted",
            "LOI_0_100",
            "geometric_patch_weight_wi",
        ]

    with (ROOT / "data" / "processed" / "integrated" / "hBN_cross_domain_descriptor_long.csv").open(
        newline="", encoding="utf-8"
    ) as f:
        assert next(csv.reader(f)) == ["domain", "descriptor", "value", "unit", "source_path", "note"]

    with (ROOT / "data" / "processed" / "integrated" / "hBN_cross_domain_heatmap_matrix.csv").open(
        newline="", encoding="utf-8"
    ) as f:
        assert next(csv.reader(f)) == ["domain", "descriptor_1", "descriptor_2", "descriptor_3", "descriptor_4"]

    with (ROOT / "data" / "processed" / "integrated" / "hBN_cross_domain_pca_projection.csv").open(
        newline="", encoding="utf-8"
    ) as f:
        assert next(csv.reader(f)) == [
            "domain",
            "PC1",
            "PC2",
            "explained_variance_ratio_PC1",
            "explained_variance_ratio_PC2",
        ]

    assert not (ROOT / "data" / "raw" / "Raman" / "hBN-3(1).txt").exists()
    assert not (ROOT / "data" / "raw" / "XPS" / "HBN(1).xls").exists()
    tem_image = ROOT / "data" / "raw" / "TEM" / "OneView 200kV 800kX 39972.jpg"
    assert tem_image.exists()
    assert tem_image.stat().st_size > 0


def test_gcd_outputs_use_unified_bounded_fit_contract() -> None:
    final_summary = pd.read_csv(
        ROOT / "data" / "processed" / "GCD" / "tables" / "hBN_GCD_bounded_fit_summary_final.csv"
    )
    public_summary = pd.read_csv(
        ROOT / "data" / "processed" / "GCD" / "tables" / "hBN_GCD_fit_summary_J1_to_J5.csv"
    )
    pd.testing.assert_frame_equal(final_summary, public_summary)

    expected_csp = [
        99.0876457235766,
        90.02908328210498,
        82.62651540927513,
        70.89006794626557,
        61.27702420765139,
    ]
    assert final_summary["Csp_base_F_g^-1"].to_list() == pytest.approx(expected_csp, rel=1e-8)
    assert final_summary["Csp_base_F_g^-1"].is_monotonic_decreasing
    assert set(final_summary["window_selection_method"]) == {"flattest_linear_post_ir_window"}
    assert (final_summary["flatness_score"] > 0).all()
    assert final_summary["tau_base_s"].to_list() == pytest.approx([9.0] * 5)
    assert final_summary["A_base"].to_list() == pytest.approx([0.0] * 5)
    assert final_summary["Rs_base_ohm"].to_list() == pytest.approx(
        [
            0.00986700312913102,
            0.009732716467688456,
            0.009627444607500105,
            0.009465826381994863,
            0.009339393602113992,
        ],
        rel=1e-8,
    )
    assert final_summary["effective_IR_drop_base_mV"].to_list() == pytest.approx(
        [
            9.86700312913102,
            19.465432935376914,
            28.882333822500314,
            37.86330552797945,
            46.69696801056996,
        ],
        rel=1e-8,
    )

    gcd_text = "\n".join(
        path.read_text(encoding="utf-8")
        for path in sorted((ROOT / "scripts" / "gcd").glob("*.py"))
    )
    assert "--use-reference-final-table" not in gcd_text
    assert "table_s3_path" not in gcd_text
    assert "slope_only_parameters_from_window" in gcd_text
    assert "flattest_linear_post_ir_window" in gcd_text
    assert "SLOPE_MIN_VOLTAGE_SPAN_FRACTION" in gcd_text
    assert "\"raw_V\"" in gcd_text
    assert "masked_domain_reporting_resistance" in gcd_text

    readme_text = (ROOT / "README.md").read_text(encoding="utf-8")
    note_s9_text = (ROOT / "docs" / "Note_S9_GitHub_v4.md").read_text(encoding="utf-8")
    assert "--use-reference-final-table" not in readme_text + note_s9_text
    assert "flattest quasi-linear segment" in readme_text
    assert "one shared algorithm" in note_s9_text


def test_xps_outputs_use_si_aligned_c1s_reporting_values() -> None:
    expected_labels = ["C–C", "C–O", "C=O", "O–C=O"]
    expected_centers = [284.59, 285.50, 286.35, 288.43]
    expected_fwhm = [1.50, 0.85, 1.45, 6.00]
    expected_area = [0.7161, 0.0967, 0.0923, 0.0949]

    descriptor = pd.read_csv(ROOT / "data" / "processed" / "XPS" / "XPS_corrected_13peak_descriptor_matrix.csv")
    c1s_descriptor = descriptor[descriptor["Region"].eq("C1s")].reset_index(drop=True)
    assert c1s_descriptor["Label"].to_list() == expected_labels
    assert c1s_descriptor["Center_eV"].to_list() == pytest.approx(expected_centers)
    assert c1s_descriptor["FWHM_eV"].to_list() == pytest.approx(expected_fwhm)
    assert c1s_descriptor["Area_fraction"].to_list() == pytest.approx(expected_area)

    profile = pd.read_csv(ROOT / "data" / "processed" / "XPS" / "C1s_profile_fit_summary.csv")
    assert profile["Label"].to_list() == expected_labels
    assert profile["Center_eV"].to_list() == pytest.approx(expected_centers)
    assert profile["FWHM_eV"].to_list() == pytest.approx(expected_fwhm)
    assert profile["Area_fraction"].to_list() == pytest.approx(expected_area)

    fit_summary = pd.read_csv(ROOT / "data" / "processed" / "XPS" / "XPS_C1s_fit_summary.csv")
    assert fit_summary["Center (eV)"].to_list() == pytest.approx(expected_centers)
    assert fit_summary["Approx. FWHM (eV)"].to_list() == pytest.approx(expected_fwhm)
    assert fit_summary["Relative area fraction"].to_list() == pytest.approx(expected_area)

    assignments = pd.read_csv(ROOT / "data" / "processed" / "XPS" / "XPS_13_peak_assignments.csv")
    c1s_assignments = assignments[assignments["Region"].eq("C1s")].reset_index(drop=True)
    assert c1s_assignments["Label"].to_list() == expected_labels
    assert c1s_assignments["Area_fraction"].to_list() == pytest.approx(expected_area)

    ml_script = (ROOT / "scripts" / "xps" / "03_ml_analysis_ns15_ef.py").read_text(encoding="utf-8")
    assert "component_area_fractions" not in ml_script
    assert "fit_fixed_centers" not in ml_script


def test_qc_circuit_notebook_is_manuscript_facing() -> None:
    notebook = json.loads(
        (ROOT / "scripts" / "qc_circuit" / "generate_gcd_eis_qc_circuits.ipynb").read_text(
            encoding="utf-8"
        )
    )
    text = "\n".join("".join(cell.get("source", [])) for cell in notebook["cells"])

    assert "GCD and EIS quantum-circuit figure notebook" in text
    assert "GCD discharge-window QAOA circuit" in text
    assert "EIS 21-qubit QUBO/QAOA circuit" in text
    assert "從 PPT" not in text
    assert "投影片" not in text
    assert "placeholder" not in text.lower()
    assert "%pip install" not in text
    assert "scripts/" + "QC Circuit" not in text

    code_cells = [cell for cell in notebook["cells"] if cell.get("cell_type") == "code"]
    assert code_cells
    assert all(not cell.get("outputs") for cell in code_cells)
    assert all("execution" not in cell.get("metadata", {}) for cell in code_cells)


def test_integrated_heatmap_validation_errors_are_clear(tmp_path: Path) -> None:
    module = load_integrated_module()
    rows = []
    for domain in module.DOMAINS:
        for idx in range(module.EXPECTED_DESCRIPTORS_PER_DOMAIN):
            rows.append(
                {
                    "domain": domain,
                    "descriptor": f"descriptor_{idx + 1}",
                    "value": float(idx + 1),
                    "unit": "-",
                    "source_path": f"data/processed/{domain}/source.csv",
                    "note": "",
                }
            )

    valid_df = module.pd.DataFrame(rows)
    heatmap = module.build_heatmap_matrix(valid_df)
    assert heatmap.shape == (len(module.DOMAINS), module.EXPECTED_DESCRIPTORS_PER_DOMAIN)

    missing_domain_df = valid_df[valid_df["domain"] != "EIS"]
    with pytest.raises(ValueError, match="Missing cross-domain descriptor rows.*EIS"):
        module.build_heatmap_matrix(missing_domain_df)

    uneven_df = valid_df[
        ~((valid_df["domain"] == "CV") & (valid_df["descriptor"] == "descriptor_4"))
    ]
    with pytest.raises(ValueError, match="Expected 4 descriptors for domain CV, found 3"):
        module.build_heatmap_matrix(uneven_df)

    one_component_projection = module.pd.DataFrame(
        {
            "domain": ["TEM", "Raman"],
            "PC1": [0.25, -0.25],
            "explained_variance_ratio_PC1": [1.0, 1.0],
            "explained_variance_ratio_PC2": [module.np.nan, module.np.nan],
        }
    )
    one_axis_heatmap = module.pd.DataFrame({"descriptor_1": [-1.0, 1.0]}, index=["TEM", "Raman"])
    module.make_figures(one_axis_heatmap, one_component_projection, tmp_path)
    assert (tmp_path / "hBN_cross_domain_pca.png").exists()


def test_cv_current_heatmap_uses_split_cv_branches() -> None:
    module = load_cv_current_heatmap_module()
    curves = module.extract_all_curves(module.load_cv_dataframe())
    E_grid = module.common_potential_grid(curves, "anodic")
    cathodic_matrix = module.interpolate_branch_matrix(curves, "cathodic", E_grid)
    anodic_matrix = module.interpolate_branch_matrix(curves, "anodic", E_grid)
    intensity_matrix = module.branch_rms_intensity_matrix(cathodic_matrix, anodic_matrix)

    assert cathodic_matrix.shape == anodic_matrix.shape == (len(module.SCAN_RATES), E_grid.size)
    assert intensity_matrix.shape == cathodic_matrix.shape
    assert module.np.all(module.np.diff(E_grid) > 0)

    # The raw CV loops contain a large negative cathodic feature and a large positive
    # anodic feature. A collapsed full-loop interpolation produces only a smooth
    # scan-rate gradient and fails these range checks.
    assert cathodic_matrix.min() < -4.5
    assert anodic_matrix.max() > 8.5
    assert cathodic_matrix.std(axis=1).max() > 1.0
    assert anodic_matrix.std(axis=1).max() > 1.5
    assert 5.5 < intensity_matrix.max() < 6.8
    assert intensity_matrix[:, (0.28 < E_grid) & (E_grid < 0.31)].mean() > 2.0


def test_package_excludes_cache_and_os_metadata() -> None:
    result = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=SCRIPT_TIMEOUT_SECONDS,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    paths = [Path(line) for line in result.stdout.splitlines() if line.strip()]
    forbidden_names = {".DS_Store", "__MACOSX", "__pycache__"}
    for path in paths:
        assert path.name not in forbidden_names, path.as_posix()
        assert not path.name.startswith("._"), path.as_posix()
        assert path.suffix != ".pyc", path.as_posix()
        assert not path.as_posix().startswith("outputs/figures/integrated/"), path.as_posix()

    for artifact in (ROOT / "data" / "processed").rglob("*"):
        if artifact.is_file():
            assert not artifact.name.startswith(("Figure", "Table")), artifact.relative_to(ROOT).as_posix()


def test_tem_pipeline_rejects_invalid_geometric_eps(tmp_path: Path) -> None:
    for invalid_eps in ("0", "-0.1", "nan", "inf"):
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts" / "tem" / "00_tem_patch_ensemble_analysis.py"),
                "--geometric-eps",
                invalid_eps,
                "--output-dir",
                str(tmp_path / f"tem-{invalid_eps}"),
            ],
            capture_output=True,
            text=True,
            timeout=SCRIPT_TIMEOUT_SECONDS,
        )

        assert result.returncode != 0
        assert "--geometric-eps" in result.stderr
        assert "must be a positive finite float" in result.stderr


def test_tem_pipeline_uses_reference_descriptors_and_real_roi_boxes(tmp_path: Path) -> None:
    output_dir = tmp_path / "tem"
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "tem" / "00_tem_patch_ensemble_analysis.py"),
            "--repo-root",
            str(ROOT),
            "--source-image",
            "data/raw/TEM/OneView 200kV 800kX 39972.jpg",
            "--roi-csv",
            "data/raw/TEM/roi_boxes_template.csv",
            "--output-dir",
            str(output_dir),
        ],
        capture_output=True,
        text=True,
        timeout=SCRIPT_TIMEOUT_SECONDS,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "Descriptor source:" in result.stdout
    assert "data\\processed\\TEM\\TEM_descriptors.csv" in result.stdout or "data/processed/TEM/TEM_descriptors.csv" in result.stdout
    assert (output_dir / "TEM_roi_selection_overlay.png").exists()
    assert (output_dir / "TEM_descriptors.csv").read_text(encoding="utf-8").strip() == (
        ROOT / "data" / "processed" / "TEM" / "TEM_descriptors.csv"
    ).read_text(encoding="utf-8").strip()


def test_eis_scripts_write_outputs_in_temporary_project(tmp_path: Path) -> None:
    project = tmp_path / "project"
    shutil.copytree(ROOT / "scripts", project / "scripts")

    raw_eis = project / "data" / "raw" / "EIS"
    raw_eis.mkdir(parents=True)
    shutil.copy2(ROOT / "data" / "raw" / "EIS" / "hbn_EIS_1.csv", raw_eis / "hbn_EIS_1.csv")

    expected_outputs = {
        "eis/01_eis_classical_anchor_fit.py": [
            "data/processed/EIS/classical_fit/hBN_EIS_final_fitting_curve.csv",
            "data/processed/EIS/classical_fit/hBN_EIS_final_fit_parameters.csv",
            "data/processed/EIS/classical_fit/hBN_EIS_residuals_final_model.csv",
        ],
        "eis/02_eis_quantum_comparison_from_anchor.py": [
            "data/processed/EIS/quantum_branches/hBN_EIS_quantum_branch_parameters.csv",
            "data/processed/EIS/quantum_branches/hBN_EIS_quantum_branch_metrics.csv",
            "data/processed/EIS/quantum_branches/hBN_EIS_quantum_branch_overlays.csv",
            "data/processed/EIS/quantum_branches/hBN_parameter_deviations.csv",
        ],
        "eis/03_eis_surrogate_qaoa_landscape.py": [
            "data/processed/EIS/qaoa_landscapes/hBN_surrogate_slice.csv",
            "data/processed/EIS/qaoa_landscapes/hBN_qaoa_coarse_landscape.csv",
            "data/processed/EIS/qaoa_landscapes/hBN_qaoa_refined_landscape.csv",
            "data/processed/EIS/qaoa_landscapes/hBN_surrogate_fidelity.csv",
        ],
    }

    for script_name, output_paths in expected_outputs.items():
        try:
            result = subprocess.run(
                [sys.executable, str(project / "scripts" / script_name)],
                cwd=project,
                capture_output=True,
                text=True,
                timeout=SCRIPT_TIMEOUT_SECONDS,
            )
        except subprocess.TimeoutExpired as exc:
            stdout = exc.stdout or ""
            stderr = exc.stderr or ""
            raise AssertionError(
                f"{script_name} timed out after {SCRIPT_TIMEOUT_SECONDS} seconds\n"
                f"stdout:\n{stdout}\n\nstderr:\n{stderr}"
            ) from exc
        assert result.returncode == 0, result.stdout + result.stderr

        for output_path in output_paths:
            output = project / output_path
            assert output.exists(), output_path
            assert output.stat().st_size > 0, output_path

    qaoa_output = project / "data" / "processed" / "EIS" / "qaoa_landscapes"
    assert not list(qaoa_output.glob("Figure_12*.csv"))
    assert not (qaoa_output / "hBN_qaoa_summary.csv").exists()

    classical_parameters = project / "data" / "processed" / "EIS" / "classical_fit" / "hBN_EIS_final_fit_parameters.csv"
    with classical_parameters.open(newline="", encoding="utf-8") as f:
        assert next(csv.reader(f)) == ["Parameter", "Value", "Unit", "Meaning"]

    quantum_deviations = project / "data" / "processed" / "EIS" / "quantum_branches" / "hBN_parameter_deviations.csv"
    with quantum_deviations.open(newline="", encoding="utf-8") as f:
        assert next(csv.reader(f)) == [
            "Parameter",
            "Classical",
            "Continuous_branch",
            "Discrete_branch",
            "Absolute_deviation_continuous",
            "Absolute_deviation_discrete",
            "Relative_deviation_continuous_percent",
            "Relative_deviation_discrete_percent",
            "Unit",
        ]


def test_full_shared_objective_pipeline_writes_manuscript_slices(tmp_path: Path) -> None:
    project = tmp_path / "project"
    shutil.copytree(ROOT / "scripts", project / "scripts")

    raw_eis = project / "data" / "raw" / "EIS"
    raw_eis.mkdir(parents=True)
    shutil.copy2(ROOT / "data" / "raw" / "EIS" / "hbn_EIS_1.csv", raw_eis / "hbn_EIS_1.csv")

    output_dir = project / "outputs" / "eis_shared_objective"
    env = {**os.environ, "MPLBACKEND": "Agg"}
    try:
        result = subprocess.run(
            [
                sys.executable,
                str(project / "scripts" / "eis" / "04_eis_shared_objective_full_pipeline.py"),
                "--input",
                "data/raw/EIS/hbn_EIS_1.csv",
                "--output",
                "outputs/eis_shared_objective",
                "--qaoa-grid",
                "3",
                "--shots",
                "256",
            ],
            cwd=project,
            capture_output=True,
            text=True,
            timeout=FULL_PIPELINE_TIMEOUT_SECONDS,
            env=env,
        )
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        raise AssertionError(
            f"04_eis_shared_objective_full_pipeline.py timed out after "
            f"{FULL_PIPELINE_TIMEOUT_SECONDS} seconds\nstdout:\n{stdout}\n\nstderr:\n{stderr}"
        ) from exc
    assert result.returncode == 0, result.stdout + result.stderr

    expected_outputs = [
        "R1_Q1_slice.csv",
        "Rs_alpha1_slice.csv",
        "surrogate_2d_slices.csv",
        "qaoa_p1_landscape.csv",
        "decoded_bitstring_summary.json",
        "nyquist_classical_discrete.png",
        "nyquist_classical_continuous_discrete.png",
        "qaoa_p1_landscape.png",
        "eis_fitted_spectra_classical_continuous_discrete.csv",
    ]
    for output_name in expected_outputs:
        output = output_dir / output_name
        assert output.exists(), output_name
        assert output.stat().st_size > 0, output_name

    qaoa_landscape = pd.read_csv(output_dir / "qaoa_p1_landscape.csv")
    assert set(qaoa_landscape["mode"]) == {"exact_statevector_p1"}
    anchor_parameters = pd.read_csv(output_dir / "classical_anchor_parameters.csv")
    rct = float(anchor_parameters.loc[anchor_parameters["parameter"] == "Rct", "classical_anchor_value"].iloc[0])
    assert 20.0 < rct < 50.0
    branch_loss = pd.read_csv(output_dir / "branch_loss_comparison.csv")
    assert {"unweighted_true_SSE", "relative_weighted_SSE", "surrogate_energy"}.issubset(branch_loss.columns)
    assert "continuous_surrogate_minimum" in set(branch_loss["branch"])
    spectra = pd.read_csv(output_dir / "eis_fitted_spectra_classical_continuous_discrete.csv")
    assert {"Zreal_continuous_ohm", "minus_Zimag_continuous_ohm", "phase_continuous_deg"}.issubset(spectra.columns)
    with (output_dir / "decoded_bitstring_summary.json").open(encoding="utf-8") as f:
        decoded_summary = json.load(f)
    assert decoded_summary["qaoa_summary"]["qaoa_landscape_mode"] == "exact_statevector_p1"
    assert decoded_summary["continuous_branch"]["source"] == "continuous_surrogate_minimum"
    assert decoded_summary["loss_definition"] == "unweighted complex SSE over real and imaginary residuals"

    with (output_dir / "R1_Q1_slice.csv").open(newline="", encoding="utf-8") as f:
        assert next(csv.reader(f)) == [
            "slice",
            "param_x",
            "param_y",
            "internal_param_x",
            "internal_param_y",
            "normalized_x",
            "normalized_y",
            "surrogate_energy",
        ]


def test_full_shared_objective_pipeline_rejects_invalid_cli_values(tmp_path: Path) -> None:
    env = {**os.environ, "MPLBACKEND": "Agg"}
    try:
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts" / "eis" / "04_eis_shared_objective_full_pipeline.py"),
                "--input",
                str(ROOT / "data" / "raw" / "EIS" / "hbn_EIS_1.csv"),
                "--output",
                str(tmp_path / "outputs"),
                "--bits-per-parameter",
                "0",
            ],
            capture_output=True,
            text=True,
            timeout=SCRIPT_TIMEOUT_SECONDS,
            env=env,
        )
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        raise AssertionError(
            f"04_eis_shared_objective_full_pipeline.py invalid CLI check timed out after "
            f"{SCRIPT_TIMEOUT_SECONDS} seconds\nstdout:\n{stdout}\n\nstderr:\n{stderr}"
        ) from exc

    assert result.returncode != 0
    assert "--bits-per-parameter" in result.stderr
    assert "must be greater than 0" in result.stderr
