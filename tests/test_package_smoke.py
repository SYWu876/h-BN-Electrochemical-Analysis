from __future__ import annotations

import py_compile
import re
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_readme_and_citation_contact_email_match() -> None:
    readme_text = (ROOT / "README.md").read_text(encoding="utf-8")
    citation_text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")

    readme_match = re.search(r"^Email:\s*([^\s]+)\s*$", readme_text, flags=re.MULTILINE)
    citation_emails = re.findall(r"^\s*email:\s*([^\s]+)\s*$", citation_text, flags=re.MULTILINE)

    assert readme_match is not None
    assert citation_emails
    assert readme_match.group(1) in citation_emails


def test_analysis_scripts_compile(tmp_path: Path) -> None:
    for script in sorted((ROOT / "scripts").glob("*.py")):
        py_compile.compile(str(script), cfile=str(tmp_path / f"{script.stem}.pyc"), doraise=True)


def test_eis_scripts_write_outputs_in_temporary_project(tmp_path: Path) -> None:
    project = tmp_path / "project"
    shutil.copytree(ROOT / "scripts", project / "scripts")

    raw_eis = project / "data" / "raw" / "EIS"
    raw_eis.mkdir(parents=True)
    shutil.copy2(ROOT / "data" / "raw" / "EIS" / "hbn_EIS_1.csv", raw_eis / "hbn_EIS_1.csv")

    expected_outputs = {
        "01_classical_eis_fit.py": [
            "data/processed/EIS/classical_fit/hBN_EIS_final_fitting_curve.csv",
            "data/processed/EIS/classical_fit/hBN_EIS_final_fit_parameters.csv",
            "data/processed/EIS/classical_fit/hBN_EIS_residuals_final_model.csv",
        ],
        "02_quantum_branch_comparison.py": [
            "data/processed/EIS/quantum_branches/hBN_EIS_quantum_branch_parameters.csv",
            "data/processed/EIS/quantum_branches/hBN_EIS_quantum_branch_metrics.csv",
            "data/processed/EIS/quantum_branches/hBN_EIS_quantum_branch_overlays.csv",
        ],
        "03_surrogate_qaoa_landscape.py": [
            "data/processed/EIS/qaoa_landscapes/hBN_surrogate_slice.csv",
            "data/processed/EIS/qaoa_landscapes/hBN_qaoa_coarse_landscape.csv",
            "data/processed/EIS/qaoa_landscapes/hBN_qaoa_refined_landscape.csv",
            "data/processed/EIS/qaoa_landscapes/hBN_surrogate_fidelity.csv",
        ],
    }

    for script_name, output_paths in expected_outputs.items():
        result = subprocess.run(
            [sys.executable, str(project / "scripts" / script_name)],
            cwd=project,
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert result.returncode == 0, result.stdout + result.stderr

        for output_path in output_paths:
            output = project / output_path
            assert output.exists(), output_path
            assert output.stat().st_size > 0, output_path

    qaoa_output = project / "data" / "processed" / "EIS" / "qaoa_landscapes"
    assert not list(qaoa_output.glob("Figure_12*.csv"))
    assert not (qaoa_output / "hBN_qaoa_summary.csv").exists()
